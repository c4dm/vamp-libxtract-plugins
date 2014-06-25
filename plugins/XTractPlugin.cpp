/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugins using Jamie Bullock's
    libxtract audio feature extraction library.

    Centre for Digital Music, Queen Mary, University of London.
    This file copyright 2006-2012 Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "XTractPlugin.h"

#include <cassert>
#include <cstdio>
#include <math.h>
#include <stdio.h>

using std::cerr;
using std::endl;
using std::string;

xtract_function_descriptor_t *
XTractPlugin::m_xtDescriptors = 0;

int
XTractPlugin::m_xtDescRefCount = 0;

XTractPlugin::XTractPlugin(unsigned int xtFeature, float inputSampleRate) :
    Plugin(inputSampleRate),
    m_xtFeature(xtFeature),
    m_channels(0),
    m_stepSize(0),
    m_blockSize(0),
    m_resultBuffer(0),
    m_peakThreshold(10),
    m_rolloffThreshold(90),
    m_harmonicThreshold(.1),
    m_minFreq(80),
    m_maxFreq(18000),
    m_coeffCount(40),
    m_highestCoef(20),
    m_lowestCoef(0),
    m_mfccFilters(0),
    m_mfccStyle((int)XTRACT_EQUAL_GAIN),
    m_spectrumType((int)XTRACT_MAGNITUDE_SPECTRUM),
    m_dc(0),
    m_normalise(0),
    m_barkBandLimits(0),
    m_outputBinCount(0),
    m_initialised(false)
{
    if (m_xtDescRefCount++ == 0) {
        m_xtDescriptors =
            (xtract_function_descriptor_t *)xtract_make_descriptors();
    }
}

XTractPlugin::~XTractPlugin()
{
    if (m_mfccFilters) {
        for (size_t i = 0; i < m_coeffCount; ++i) {
            delete[] m_mfccFilters[i];
        } 
        delete[] m_mfccFilters;
    }
    if (m_barkBandLimits) {
        delete[] m_barkBandLimits;
    }
    if (m_resultBuffer) {
        delete[] m_resultBuffer;
    }

    if (--m_xtDescRefCount == 0) {
        xtract_free_descriptors(m_xtDescriptors);
    }
}

string
XTractPlugin::getIdentifier() const
{
    return xtDescriptor()->algo.name;
}

string
XTractPlugin::getName() const
{
    return xtDescriptor()->algo.p_name;
}

string
XTractPlugin::getDescription() const
{
    return xtDescriptor()->algo.p_desc;
}
    

string
XTractPlugin::getMaker() const
{
    return "libxtract by Jamie Bullock (plugin by Chris Cannam)";
}

int
XTractPlugin::getPluginVersion() const
{
    return 4;
}

string
XTractPlugin::getCopyright() const
{
    string text = "Copyright 2006-2012 Jamie Bullock, plugin Copyright 2006-2012 Queen Mary, University of London. ";

    string method = "";

    method += xtDescriptor()->algo.author;

    if (method != "") {
        int year = xtDescriptor()->algo.year;
        if (year != 0) {
            char yearstr[12];
            sprintf(yearstr, " (%d)", year);
            method += yearstr;
        }
        text += "Method from " + method + ". ";
    }

    text += "Distributed under the GNU General Public License";
    return text;
}

XTractPlugin::InputDomain
XTractPlugin::getInputDomain() const
{
	
    if (xtDescriptor()->data.format == XTRACT_AUDIO_SAMPLES)
	return TimeDomain;
    else
	return FrequencyDomain;
}
   

bool XTractPlugin::m_anyInitialised = false;

bool
XTractPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{

    int donor = *(xtDescriptor()->argv.donor),
	data_format = xtDescriptor()->data.format;

    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    if (blockSize != getPreferredBlockSize()) {
        cerr << "XTractPlugin::initialise: ERROR: "
             << "Only the standard block size of " << getPreferredBlockSize()
             << " is supported (owing to global FFT initialisation requirements)" << endl;
        return false;
    }

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    if (!m_anyInitialised) {
        m_anyInitialised = true;
        // initialise libxtract
        xtract_init_fft(m_blockSize, XTRACT_SPECTRUM);
        xtract_init_fft(m_blockSize, XTRACT_AUTOCORRELATION_FFT);
        xtract_init_fft(m_blockSize, XTRACT_DCT);
        xtract_init_fft(m_blockSize, XTRACT_MFCC);
    }        

    if (donor == XTRACT_INIT_MFCC) {

        m_mfccFilters = new double *[m_coeffCount];
        for (size_t i = 0; i < m_coeffCount; ++i) {
            m_mfccFilters[i] = new double[m_blockSize];
        }

        int error = (int)xtract_init_mfcc(m_blockSize, m_inputSampleRate/2,
                                          m_mfccStyle, m_minFreq, m_maxFreq,
                                          m_coeffCount, m_mfccFilters);
        if (error != XTRACT_SUCCESS) {
            cerr << "XTractPlugin::initialise: ERROR: "
                 << "xtract_init_mfcc returned error code " << error << endl;
            return false;
        }

    } else if (donor == XTRACT_BARK_COEFFICIENTS ||
               donor == XTRACT_INIT_BARK ||
               data_format == XTRACT_BARK_COEFFS) {

        m_barkBandLimits = new int[XTRACT_BARK_BANDS];

        /*int error = *(int)*/xtract_init_bark(m_blockSize, m_inputSampleRate,
                                          m_barkBandLimits);
//        if (error != SUCCESS) {
//            cerr << "XTractPlugin::initialise: ERROR: "
//                 << "xtract_init_bark returned error code " << error << endl;
//            return false;
//        }
    }

    switch (m_xtFeature) {
	case XTRACT_SPECTRUM:	     
	    m_outputBinCount = m_blockSize / 2 + (m_dc ? 1 : 0); break;
	case XTRACT_HARMONIC_SPECTRUM:	     
	case XTRACT_PEAK_SPECTRUM:           
	    m_outputBinCount = m_blockSize / 2; break;
	case XTRACT_DCT:                 
	case XTRACT_AUTOCORRELATION_FFT: 
	case XTRACT_AUTOCORRELATION:     
	case XTRACT_AMDF:                
	case XTRACT_ASDF:                
	    m_outputBinCount = m_blockSize; break;
	case XTRACT_MFCC:                
	    m_outputBinCount = (m_highestCoef - m_lowestCoef)+1; break;
	case XTRACT_BARK_COEFFICIENTS:   
	    m_outputBinCount = XTRACT_BARK_BANDS; break;
	default:			     
	    m_outputBinCount = 1; break;
    }

    m_outputDescriptors.clear();
    setupOutputDescriptors();

    m_initialised = true;

    return true;
}

void
XTractPlugin::reset()
{
}

size_t
XTractPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
XTractPlugin::getMaxChannelCount() const
{
    return 1;
}

size_t
XTractPlugin::getPreferredStepSize() const
{
    if (getInputDomain() == FrequencyDomain) {
        return getPreferredBlockSize();
    } else {
        return getPreferredBlockSize() / 2;
    }
}

size_t
XTractPlugin::getPreferredBlockSize() const
{
    return 1024;
}

XTractPlugin::ParameterList
XTractPlugin::getParameterDescriptors() const
{
    ParameterList list;
    ParameterDescriptor desc;

    if (m_xtFeature == XTRACT_MFCC) {

        desc.identifier = "minfreq";
        desc.name = "Minimum Frequency";
        desc.minValue = 0;
        desc.maxValue = m_inputSampleRate / 2;
        desc.defaultValue = 80;
        desc.isQuantized = false;
        desc.unit = "Hz";
        list.push_back(desc);

        desc.identifier = "maxfreq";
        desc.name = "Maximum Frequency";
        desc.defaultValue = 18000;
        if (desc.defaultValue > m_inputSampleRate * 0.875) {
            desc.defaultValue = m_inputSampleRate * 0.875;
        }
        list.push_back(desc);

        desc.identifier = "bands";
        desc.name = "# Mel Frequency Bands";
        desc.minValue = 10;
        desc.maxValue = 80;
        desc.defaultValue = 40;
        desc.unit = "";
        desc.isQuantized = true;
        desc.quantizeStep = 1;
        list.push_back(desc);

        desc.identifier = "lowestcoef";
        desc.name = "Lowest Coefficient Returned";
        desc.minValue = 0;
        desc.maxValue = 80;
        desc.defaultValue = 0;
        desc.unit = "";
        desc.isQuantized = true;
        desc.quantizeStep = 1;
        list.push_back(desc);

        desc.identifier = "highestcoef";
        desc.name = "Highest Coefficient Returned";
        desc.minValue = 0;
        desc.maxValue = 80;
        desc.defaultValue = 20;
        desc.unit = "";
        desc.isQuantized = true;
        desc.quantizeStep = 1;
        list.push_back(desc);

        desc.identifier = "style";
        desc.name = "MFCC Type";
        desc.minValue = 0;
        desc.maxValue = 1;
        desc.defaultValue = 0;
        desc.valueNames.push_back("Equal Gain");
        desc.valueNames.push_back("Equal Area");
        list.push_back(desc);
    }

    if (m_xtFeature == XTRACT_SPECTRUM) {

        desc.identifier = "spectrumtype";
        desc.name = "Type";
        desc.minValue = 0;
        desc.maxValue = 3;
        desc.defaultValue = int(XTRACT_MAGNITUDE_SPECTRUM);
        desc.isQuantized = true;
        desc.quantizeStep = 1;
        desc.valueNames.push_back("Magnitude Spectrum");
        desc.valueNames.push_back("Log Magnitude Spectrum");
        desc.valueNames.push_back("Power Spectrum");
        desc.valueNames.push_back("Log Power Spectrum");
        list.push_back(desc);

        desc.identifier = "dc";
        desc.name = "Include DC";
        desc.maxValue = 1;
        desc.defaultValue = 0;
        desc.valueNames.clear();
        list.push_back(desc);

        desc.identifier = "normalise";
        desc.name = "Normalise";
        list.push_back(desc);
    }

    if (needPeakThreshold()) {
        
        desc.identifier = "peak-threshold";
        desc.name = "Peak Threshold";
        desc.minValue = 0;
        desc.maxValue = 100;
        desc.defaultValue = 10; /* Threshold as % of maximum peak found */
        desc.isQuantized = false;
        desc.valueNames.clear();
        desc.unit = "%";
        list.push_back(desc);

    } 
    
    if (needRolloffThreshold()) {

        desc.identifier = "rolloff-threshold";
        desc.name = "Rolloff Threshold";
        desc.minValue = 0;
        desc.maxValue = 100;
        desc.defaultValue = 90; /* Freq below which 90% of energy is */
        desc.isQuantized = false;
        desc.valueNames.clear();
        desc.unit = "%";
        list.push_back(desc);

    }

    if (needHarmonicThreshold()) {

	desc.identifier = "harmonic-threshold";
        desc.name = "Harmonic Threshold";
        desc.minValue = 0;
        desc.maxValue = 1.0;
        desc.defaultValue = .1; /* Distance from nearesst harmonic number */
        desc.isQuantized = false;
        desc.valueNames.clear();
        desc.unit = "";
        list.push_back(desc);
    }

    return list;
}

float
XTractPlugin::getParameter(string param) const
{
    if (m_xtFeature == XTRACT_MFCC) {
        if (param == "minfreq") return m_minFreq;
        if (param == "maxfreq") return m_maxFreq;
        if (param == "bands") return m_coeffCount;
        if (param == "lowestcoef") return m_lowestCoef;
        if (param == "highestcoef") return m_highestCoef;
        if (param == "style") return m_mfccStyle;
    }

    if (m_xtFeature == XTRACT_SPECTRUM) {
        if (param == "spectrumtype") return m_spectrumType;
        if (param == "dc") return m_dc;
        if (param == "normalise") return m_normalise;
    }

    if (param == "peak-threshold") return m_peakThreshold;
    if (param == "rolloff-threshold") return m_rolloffThreshold;
    if (param == "harmonic-threshold") return m_harmonicThreshold;

    return 0.f;
}

void
XTractPlugin::setParameter(string param, float value)
{
    if (m_xtFeature == XTRACT_MFCC) {
        if (param == "minfreq") m_minFreq = value;
        else if (param == "maxfreq") m_maxFreq = value;
        else if (param == "bands") m_coeffCount = int(value + .1);
        else if (param == "lowestcoef"){
        	m_lowestCoef  = int(value + .1);
        	if(m_lowestCoef >= m_coeffCount) m_lowestCoef = m_coeffCount - 1;
        	if(m_lowestCoef > m_highestCoef) m_lowestCoef = m_highestCoef;
        }
        else if (param == "highestcoef"){
        	m_highestCoef  = int(value + .1);
        	if(m_highestCoef >= m_coeffCount) m_highestCoef = m_coeffCount - 1;
        	if(m_highestCoef < m_lowestCoef) m_highestCoef = m_lowestCoef;
        }
        else if (param == "style") m_mfccStyle = int(value + .1);
    }

    if (m_xtFeature == XTRACT_SPECTRUM) {
        if (param == "spectrumtype") m_spectrumType = int(value + .1);
        if (param == "dc") m_dc = int(value + .1);
        if (param == "normalise") m_normalise = int(value + .1);
    }

    if (param == "peak-threshold") m_peakThreshold = value;
    if (param == "rolloff-threshold") m_rolloffThreshold = value;
    if (param == "harmonic-threshold") m_harmonicThreshold = value;
}

XTractPlugin::OutputList
XTractPlugin::getOutputDescriptors() const
{
    if (m_outputDescriptors.empty()) {
        setupOutputDescriptors();
    }
    return m_outputDescriptors;
}

void
XTractPlugin::setupOutputDescriptors() const
{
    OutputDescriptor d;
    const xtract_function_descriptor_t *xtFd = xtDescriptor();
    d.identifier = getIdentifier();  
    d.name = getName();
    d.description = getDescription();
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = m_outputBinCount;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep;

    if (xtFd->is_scalar){
	switch(xtFd->result.scalar.unit){
	    case XTRACT_HERTZ:	    d.unit = "Hz"; break;
	    case XTRACT_DBFS:	    d.unit = "dB"; break;
	    default:	    d.unit = ""; break;
	}
    }	
    else {
	if (xtFd->result.vector.format == XTRACT_SPECTRAL){

	    d.binCount /= 2;
	    d.identifier = "amplitudes";
	    d.name = "Peak Amplitudes";
            d.description = "";
	}
    } 

    m_outputDescriptors.push_back(d);
}

bool
XTractPlugin::needPeakThreshold() const
{
    const xtract_function_descriptor_t *xtFd = xtDescriptor();

    if(m_xtFeature == XTRACT_PEAK_SPECTRUM || 
	    xtFd->data.format == XTRACT_SPECTRAL_PEAKS ||
	    xtFd->data.format == XTRACT_SPECTRAL_PEAKS_MAGNITUDES ||
	    needHarmonicThreshold()) 
	return true;
    else return false;
}

bool
XTractPlugin::needHarmonicThreshold() const
{
    const xtract_function_descriptor_t *xtFd = xtDescriptor();

    if(m_xtFeature == XTRACT_HARMONIC_SPECTRUM || 
	    xtFd->data.format == XTRACT_SPECTRAL_HARMONICS_FREQUENCIES ||
	    m_xtFeature == XTRACT_NOISINESS ||
	    xtFd->data.format == XTRACT_SPECTRAL_HARMONICS_MAGNITUDES) 
	return true;
    else return false;
}

bool
XTractPlugin::needRolloffThreshold() const
{
    if(m_xtFeature == XTRACT_ROLLOFF) 
	return true;
    else
	return false;
}

XTractPlugin::FeatureSet
XTractPlugin::process(const float *const *inputBuffers,
                      Vamp::RealTime timestamp)
{
    if (m_outputDescriptors.empty()) {
        setupOutputDescriptors();
    }

    int rbs =
        // Add 2 here to accommodate extra data for spectrum with DC
        2 + (m_outputBinCount > m_blockSize ? m_outputBinCount : m_blockSize);
    if (!m_resultBuffer) {
        m_resultBuffer = new double[rbs];
    }

    int i;

    for (i = 0; i < rbs; ++i) m_resultBuffer[i] = 0.f;

    int N = m_blockSize, M = N >> 1;

    const double *data = 0;
    double *input_d = new double[N];
    for (int i = 0; i < N; ++i) {
        input_d[i] = inputBuffers[0][i];
    }

    double *fft_temp = 0, *data_temp = 0;
    void *argv = 0;
    bool isSpectral = false;
    xtract_function_descriptor_t *xtFd = xtDescriptor();

    FeatureSet fs;

    switch (xtFd->data.format) {
	case XTRACT_AUDIO_SAMPLES:
	    data = input_d;
	    break;
	case XTRACT_SPECTRAL:
	default:
	    // All the rest are derived from the  spectrum
	    // Need same format as would be output by xtract_spectrum
	    double q = m_inputSampleRate / N;
	    fft_temp = new double[N];
	    for (int n = 1; n < N/2; ++n) {
		fft_temp[n] = sqrt(input_d[n*2]   * 
			input_d[n*2] + input_d[n*2+1] * 
			input_d[n*2+1]) / N;
		fft_temp[N-n] = (N/2 - n) * q;
	    }
	    fft_temp[0]   = fabs(input_d[0]) / N;
	    fft_temp[N/2] = fabs(input_d[N]) / N;
	    data = &fft_temp[0];
	    isSpectral = true;
	    break;
    }

    assert(m_outputBinCount > 0);

    double *result = m_resultBuffer;

    double argf[XTRACT_MAXARGS];
    argv = &argf[0];
    argf[0] = 0.f; // handy for some, e.g. lowest_value which has a threshold

    double mean, variance, sd, npartials, nharmonics;

    bool needSD, needVariance, needMean, needPeaks, 
	 needBarkCoefficients, needHarmonics, needF0, needSFM, needMax, 
	 needNumPartials, needNumHarmonics;

    int donor;

    needSD = needVariance = needMean = needPeaks =  
	needBarkCoefficients = needF0 = needHarmonics = needSFM = needMax = 
	needNumPartials = needNumHarmonics = 0;

    mean = variance = sd = npartials = nharmonics = 0.f;

    i = xtFd->argc;

    while(i--){
        if (m_xtFeature == XTRACT_BARK_COEFFICIENTS) {
            /* "BARK_COEFFICIENTS is special because argc = BARK_BANDS" */
            break;
        }
	donor = xtFd->argv.donor[i];
	switch(donor){
	    case XTRACT_STANDARD_DEVIATION: 
	    case XTRACT_SPECTRAL_STANDARD_DEVIATION:
		needSD = 1;	
		break;
	    case XTRACT_VARIANCE:	    
	    case XTRACT_SPECTRAL_VARIANCE:
		needVariance = 1; 
		break;	
	    case XTRACT_MEAN:		  
	    case XTRACT_SPECTRAL_MEAN:	
		needMean = 1;		
		break;
	    case XTRACT_F0:
	    case XTRACT_FAILSAFE_F0:
		needF0 = 1;
		break;
	    case XTRACT_FLATNESS:
		needSFM = 1;
	    case XTRACT_HIGHEST_VALUE:
		needMax = 1;
		break;
	}
    }

    if(needHarmonicThreshold() && m_xtFeature != XTRACT_HARMONIC_SPECTRUM)
	needHarmonics = needF0  = 1;

    if(needPeakThreshold()  && m_xtFeature != XTRACT_PEAK_SPECTRUM) 
	needPeaks = 1; 

    if(xtFd->data.format == XTRACT_BARK_COEFFS && 
	    m_xtFeature != XTRACT_BARK_COEFFICIENTS){
	needBarkCoefficients = 1;
    }

    if (needMean) {
	if(isSpectral)
	    xtract_spectral_mean(data, N, 0, result);
	else
	    xtract_mean(data, M, 0, result);
        mean = *result;
        *result = 0.f;
    }

    if (needVariance || needSD) {
        argf[0] = mean;
	if(isSpectral)
	    xtract_spectral_variance(data, N, argv, result);
	else
	    xtract_variance(data, M, argv, result);
        variance = *result;
        *result = 0.f;
    } 

    if (needSD) {
        argf[0] = variance;
	if(isSpectral)
	    xtract_spectral_standard_deviation(data, N, argv, result);
	else
	    xtract_standard_deviation(data, M, argv, result);
        sd = *result;
        *result = 0.f;
    }

    if (needMax) {
	xtract_highest_value(data, M, argv, result);
	argf[1] = *result;
	*result = 0.f;
    }

    if (needSD) {
        argf[0] = mean;
        argf[1] = sd;
    } else if (needVariance) {
        argf[0] = variance;
    } else if (needMean) {
        argf[0] = mean;
    }
    
    // data should be now correct for all except:
    // XTRACT_SPECTRAL_CENTROID -- N/2 magnitude peaks and N/2 frequencies
    // TONALITY -- SFM
    // TRISTIMULUS_1/2/3 -- harmonic spectrum
    // ODD_EVEN_RATIO -- harmonic spectrum
    // LOUDNESS -- Bark coefficients
    // XTRACT_HARMONIC_SPECTRUM -- peak spectrum

    // argv should be now correct for all except:
    //
    // XTRACT_ROLLOFF -- (sr/N), threshold (%)
    // XTRACT_PEAK_SPECTRUM -- (sr / N), peak threshold (%)
    // XTRACT_HARMONIC_SPECTRUM -- f0, harmonic threshold
    // XTRACT_F0 -- samplerate
    // XTRACT_MFCC -- Mel filter coefficients
    // XTRACT_BARK_COEFFICIENTS -- Bark band limits
    // XTRACT_NOISINESS -- npartials, nharmonics.
    // XTRACT_SPECTRUM -- q, spectrum type, dc, normalise

    data_temp = new double[N];

    if (m_xtFeature == XTRACT_ROLLOFF || 
        m_xtFeature == XTRACT_PEAK_SPECTRUM || needPeaks) {
        argf[0] = m_inputSampleRate / N;
	if(m_xtFeature == XTRACT_ROLLOFF) 
	    argf[1] = m_rolloffThreshold;
	else 
	    argf[1] = m_peakThreshold;
        argv = &argf[0];
    }

    if (m_xtFeature == XTRACT_SPECTRUM) {
        argf[0] = 0; // xtract_spectrum will calculate this for us
        argf[1] = m_spectrumType;
        argf[2] = m_dc;
        argf[3] = m_normalise;
        argv = &argf[0];
    }

    if (needPeaks) {
	//We only read in the magnitudes (M)
        /*int rv = */ xtract_peak_spectrum(data, M, argv, result);
        for (int n = 0; n < N; ++n) {
            data_temp[n] = result[n];
            result[n] = 0.f;
        }
        // rv not trustworthy
//        if (rv != SUCCESS) {
//            cerr << "ERROR: XTractPlugin::process: xtract_peaks failed (error code = " << rv << ")" << endl;
//            goto done;
//        }
    }

    if (needNumPartials) {
	xtract_nonzero_count(data_temp, M, NULL, &npartials);
    }

    if (needF0 || m_xtFeature == XTRACT_FAILSAFE_F0 || 
	    m_xtFeature == XTRACT_F0) {
        argf[0] = m_inputSampleRate;
        argv = &argf[0];
    }

    if (needF0) {
	xtract_failsafe_f0(&input_d[0], N, (void *)&m_inputSampleRate, result);
	argf[0] = *result;
	argv = &argf[0];
    }

    if (needSFM) {
	xtract_flatness(data, N >> 1, 0, &argf[0]);
	argv = &argf[0];
    }

    if (needHarmonics || m_xtFeature == XTRACT_HARMONIC_SPECTRUM){
       argf[1] = m_harmonicThreshold;
    }       

    if (needHarmonics){
	xtract_harmonic_spectrum(data_temp, N, argv, result);
        for (int n = 0; n < N; ++n) {
            data_temp[n] = result[n];
            result[n] = 0.f;
        }
    }

    if (needNumHarmonics) {
	xtract_nonzero_count(data_temp, M, NULL, &nharmonics);
    }

    if (m_xtFeature == XTRACT_NOISINESS) {

	argf[0] = nharmonics;
	argf[1] = npartials;
	argv = &argf[0];

    }

    if (needBarkCoefficients || m_xtFeature == XTRACT_BARK_COEFFICIENTS) {
        argv = &m_barkBandLimits[0];
    }

    xtract_mel_filter mfccFilterBank;
    if (m_xtFeature == XTRACT_MFCC) {
        mfccFilterBank.n_filters = m_coeffCount;
        mfccFilterBank.filters = m_mfccFilters;
        argv = &mfccFilterBank;
    }

    if (needBarkCoefficients) {
	
        /*int rv = */ xtract_bark_coefficients(data, 0, argv, data_temp);
//        if (rv != SUCCESS) {
//            cerr << "ERROR: XTractPlugin::process: xtract_bark_coefficients failed (error code = " << rv << ")" << endl;
//            goto done;
//        }
        data = &data_temp[0]; 
        argv = 0;
    }
    
    if (xtFd->data.format == XTRACT_SPECTRAL_HARMONICS_FREQUENCIES) {

	N = M;
	data = &data_temp[N];

    } else if (xtFd->data.format == XTRACT_SPECTRAL_HARMONICS_MAGNITUDES) {

	N = M;
	data = &data_temp[0];
   
    } 

    // If we only want spectral magnitudes, use first half of the input array
    else if(xtFd->data.format == XTRACT_SPECTRAL_MAGNITUDES ||
	    xtFd->data.format == XTRACT_SPECTRAL_PEAKS_MAGNITUDES ||
	    xtFd->data.format == XTRACT_ARBITRARY_SERIES) {
	N = M;
    }

    else if(xtFd->data.format == XTRACT_BARK_COEFFS) {

        N = XTRACT_BARK_BANDS - 1; /* Because our SR is 44100 (< 54000)*/
    }

    if (needPeaks && !needHarmonics) {

	data = &data_temp[0];

    }

    // now the main result
    xtract[m_xtFeature](data, N, argv, result);

//haveResult:
//    {
        int index = 0;

        for (size_t output = 0; output < m_outputDescriptors.size(); ++output) {

            Feature feature;
            feature.hasTimestamp = false;
            bool good = true;

            for (size_t n = 0; n < m_outputDescriptors[output].binCount; ++n) {
                double value = m_resultBuffer[index + m_lowestCoef];
                if (isnan(value) || isinf(value)) {
                    good = false;
                    index += (m_outputDescriptors[output].binCount - n);
                    break;
                }
                feature.values.push_back(value);
                ++index;
            }

            if (good) fs[output].push_back(feature);
        }
//    }
   
//done:
    delete[] fft_temp;
    delete[] data_temp;
    delete[] input_d;

//    cerr << "XTractPlugin::process returning" << endl;

    return fs;
}

XTractPlugin::FeatureSet
XTractPlugin::getRemainingFeatures()
{
    return FeatureSet();
}

