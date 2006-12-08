/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugins using Jamie Bullock's
    libxtract audio feature extraction library.

    Centre for Digital Music, Queen Mary, University of London.
    This file copyright 2006 Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "XTractPlugin.h"

#include <cassert>

#define XTRACT
#include <xtract/libxtract.h>

using std::cerr;
using std::endl;
using std::string;

XTractPlugin::XTractPlugin(unsigned int xtFeature, float inputSampleRate) :
    Plugin(inputSampleRate),
    m_xtFeature(xtFeature),
    m_channels(0),
    m_stepSize(0),
    m_blockSize(0),
    m_resultBuffer(0),
    m_threshold(10),
    m_minFreq(80),
    m_maxFreq(18000),
    m_coeffCount(20),
    m_mfccFilters(0),
    m_mfccStyle((int)EQUAL_GAIN),
    m_barkBandLimits(0),
    m_outputBinCount(0),
    m_initialised(false)
{
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
}

string
XTractPlugin::getName() const
{
    return xtract_help_strings[m_xtFeature];
}

string
XTractPlugin::getDescription() const
{
    switch (m_xtFeature) {
    case MEAN:                return "Spectral Magnitude Mean";
    case VARIANCE:            return "Spectral Magnitude Variance";
    case STANDARD_DEVIATION:  return "Spectral Magnitude Standard Deviation";
    case AVERAGE_DEVIATION:   return "Spectral Magnitude Average Deviation";
    case SKEWNESS:            return "Spectral Skewness";
    case KURTOSIS:            return "Spectral Kurtosis";
    case CENTROID:            return "Spectral Centroid";
    case IRREGULARITY_K:      return "Spectral Irregularity (Krimphoff)";
    case IRREGULARITY_J:      return "Spectral Irregularity (Jensen)";
    case TRISTIMULUS_1:       return "Tristimulus (1st order)";
    case TRISTIMULUS_2:       return "Tristimulus (2nd order)";
    case TRISTIMULUS_3:       return "Tristimulus (3rd order)";
    case SMOOTHNESS:          return "Spectral Smoothness";
    case SPREAD:              return "Spectral Spread";
    case ZCR:                 return "Zero Crossing Rate";
    case ROLLOFF:             return "Spectral Rolloff";
    case LOUDNESS:            return "Loudness";
    case FLATNESS:            return "Spectral Flatness";
    case TONALITY:            return "Tonality";
    case CREST:               return "Spectral Crest";
    case NOISINESS:           return "Noisiness";
    case RMS_AMPLITUDE:       return "RMS Amplitude";
    case INHARMONICITY:       return "Inharmonicity";
    case POWER:               return "Spectral Power";
    case ODD_EVEN_RATIO:      return "Odd/Even Harmonic Ratio";
    case SHARPNESS:           return "Spectral Sharpness";
    case SLOPE:               return "Spectral Slope";
    case LOWEST_MATCH:        return "Lowest Match";
    case HPS:                 return "Harmonic Product Spectrum";
    case F0:                  return "Fundamental Frequency";
    case FLUX:                return "Spectral Flux";
    case ATTACK_TIME:         return "Attack Time";
    case DECAY_TIME:          return "Decay Time";
    case DELTA_FEATURE:       return "Delta of Feature";
    case AUTOCORRELATION:     return "Autocorrelation (Time Domain)";
    case AMDF:                return "Average Magnitude Difference Function";
    case ASDF:                return "Average Squared Difference Function";
    case BARK_COEFFICIENTS:   return "Bark Band Coefficients";
    case PEAKS:               return "Spectral Peaks";
    case MAGNITUDE_SPECTRUM:  return "Magnitude Spectrum";
    case AUTOCORRELATION_FFT: return "Autocorrelation (FFT)";
    case MFCC:                return "Mel-Frequency Cepstral Coefficients";
    case DCT:                 return "Discrete Cosine Transform";
    default:
        cerr << "ERROR: XTractPlugin::getDescription: Unexpected feature index " << m_xtFeature << endl;
        return "<unknown>";
    }
}

string
XTractPlugin::getMaker() const
{
    return "libxtract by Jamie Bullock (plugin by Chris Cannam)";
}

int
XTractPlugin::getPluginVersion() const
{
    return 1;
}

string
XTractPlugin::getCopyright() const
{
    string text = "Copyright 2006 Jamie Bullock, plugin Copyright 2006 Queen Mary, University of London. ";

    string method = "Method from ";

    switch (m_xtFeature) {
    case IRREGULARITY_K:      method += "Krimphoff (1994)"; break;
    case IRREGULARITY_J:      method += "Jensen (1999)"; break;
    case TRISTIMULUS_1:
    case TRISTIMULUS_2:
    case TRISTIMULUS_3:       method += "Pollard and Jansson (1982)"; break;
    case SMOOTHNESS:          method += "McAdams (1999)"; break;
    case SPREAD:              method += "Casagrande (2005)"; break;
    case ROLLOFF:             method += "Bee Suan Ong (2005)"; break;
    case LOUDNESS:            method += "Moore, Glasberg et al (2005)"; break;
    case FLATNESS:            method += "Tristan Jehan (2005)"; break;
    case TONALITY:            method += "Tristan Jehan (2005)"; break;
    case NOISINESS:           method += "Tae Hong Park (2000)"; break;
    case CREST:               method += "Peeters (2003)"; break;
    case POWER:               method += "Bee Suan Ong (2005)"; break;
    case MFCC:                method += "Rabiner"; break;
    default:                  method = "";
    }

    if (method != "") text += method + ". ";
    text += "Distributed under the GNU General Public License";
    return text;
}


XTractPlugin::InputDomain
XTractPlugin::getInputDomain() const
{
    if (getXTractInputType() == TimeDomainAudio ||
        getXTractInputType() == SpectralPeaksPlusF0) return TimeDomain;
    else return FrequencyDomain;
}

bool
XTractPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    if (m_xtFeature == MFCC) {

        m_mfccFilters = new float *[m_coeffCount];
        for (size_t i = 0; i < m_coeffCount; ++i) {
            m_mfccFilters[i] = new float[m_blockSize];
        }

        int error = (int)xtract_init_mfcc(m_blockSize, m_inputSampleRate/2,
                                          m_mfccStyle, m_minFreq, m_maxFreq,
                                          m_coeffCount, m_mfccFilters);
        if (error != SUCCESS) {
            cerr << "XTractPlugin::initialise: ERROR: "
                 << "xtract_init_mfcc returned error code " << error << endl;
            return false;
        }

    } else if (m_xtFeature == BARK_COEFFICIENTS ||
               getXTractInputType() == BarkCoefficients) {

        m_barkBandLimits = new int[BARK_BANDS];

        int error = (int)xtract_init_bark(m_blockSize, m_inputSampleRate/2,
                                          m_barkBandLimits);
//        if (error != SUCCESS) {
//            cerr << "XTractPlugin::initialise: ERROR: "
//                 << "xtract_init_bark returned error code " << error << endl;
//            return false;
//        }
    }

    switch (m_xtFeature) {
    case MAGNITUDE_SPECTRUM:  m_outputBinCount = m_blockSize; break;
    case AUTOCORRELATION_FFT: m_outputBinCount = m_blockSize; break;
    case MFCC:                m_outputBinCount = m_coeffCount; break;
    case DCT:                 m_outputBinCount = m_blockSize; break;
    case AUTOCORRELATION:     m_outputBinCount = m_blockSize; break;
    case AMDF:                m_outputBinCount = m_blockSize; break;
    case ASDF:                m_outputBinCount = m_blockSize; break;
    case BARK_COEFFICIENTS:   m_outputBinCount = BARK_BANDS; break;
    case PEAKS:               m_outputBinCount = m_blockSize; break;
    default:                  m_outputBinCount = 1; break;
    }

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
        return getPreferredBlockSize() / 2;
    } else {
        return getPreferredBlockSize();
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

    if (m_xtFeature == MFCC) {

        desc.name = "minfreq";
        desc.description = "Minimum Frequency";
        desc.minValue = 0;
        desc.maxValue = m_inputSampleRate / 2;
        desc.defaultValue = 80;
        desc.isQuantized = false;
        desc.unit = "Hz";
        list.push_back(desc);

        desc.name = "maxfreq";
        desc.description = "Maximum Frequency";
        desc.defaultValue = 18000;
        if (desc.defaultValue > m_inputSampleRate * 0.875) {
            desc.defaultValue = m_inputSampleRate * 0.875;
        }
        list.push_back(desc);

        desc.name = "bands";
        desc.description = "Mel Frequency Bands";
        desc.minValue = 10;
        desc.maxValue = 30;
        desc.defaultValue = 20;
        desc.unit = "";
        desc.isQuantized = true;
        desc.quantizeStep = 1;
        list.push_back(desc);

        desc.name = "style";
        desc.description = "MFCC Type";
        desc.minValue = 0;
        desc.maxValue = 1;
        desc.defaultValue = 0;
        desc.valueNames.push_back("Equal Gain");
        desc.valueNames.push_back("Equal Area");
        list.push_back(desc);
    }

    if (needPeakThreshold()) {
        
        desc.name = "threshold";
        desc.description = "Peak Threshold";
        desc.minValue = 0;
        desc.maxValue = 100;
        desc.defaultValue = 10; //!!!???
        desc.isQuantized = false;
        desc.valueNames.clear();
        desc.unit = "%";
        list.push_back(desc);

    } else if (m_xtFeature == ROLLOFF) {

        desc.name = "threshold";
        desc.description = "Rolloff Threshold";
        desc.minValue = 0;
        desc.maxValue = 100;
        desc.defaultValue = 10; //!!!???
        desc.isQuantized = false;
        desc.valueNames.clear();
        desc.unit = "%";
        list.push_back(desc);
    }

    return list;
}

float
XTractPlugin::getParameter(string param) const
{
    if (m_xtFeature == MFCC) {
        if (param == "minfreq") return m_minFreq;
        if (param == "maxfreq") return m_maxFreq;
        if (param == "bands") return m_coeffCount;
        if (param == "style") return m_mfccStyle;
    }

    if (param == "threshold") return m_threshold;

    return 0.f;
}

void
XTractPlugin::setParameter(string param, float value)
{
    if (m_xtFeature == MFCC) {
        if (param == "minfreq") m_minFreq = value;
        else if (param == "maxfreq") m_maxFreq = value;
        else if (param == "bands") m_coeffCount = lrintf(value + .1);
        else if (param == "style") m_mfccStyle = lrintf(value + .1);
    }

    if (param == "threshold") m_threshold = value;
}

XTractPlugin::OutputList
XTractPlugin::getOutputDescriptors() const
{
    if (m_outputDescriptors.empty()) setupOutputDescriptors();
    return m_outputDescriptors;
}

void
XTractPlugin::setupOutputDescriptors() const
{
    OutputDescriptor d;
    d.name = getName();
    d.unit = "";
    d.description = getDescription();
    d.hasFixedBinCount = true;
    d.binCount = m_outputBinCount;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep;

    switch (m_xtFeature) {
    case CENTROID:            d.unit = "Hz"; break;
    case F0:                  d.unit = "Hz"; break;
    }

    if (m_xtFeature == PEAKS) {

        d.binCount /= 2;
        d.name = "frequencies";
        d.description = "Peak Frequencies";
        m_outputDescriptors.push_back(d);

        d.name = "amplitudes";
        d.description = "Peak Amplitudes";
        d.unit = "dB";

    } else if (m_xtFeature == MAGNITUDE_SPECTRUM) {

        d.binCount = d.binCount / 2 + 1;
    }

    m_outputDescriptors.push_back(d);
}

XTractPlugin::XTractInputType
XTractPlugin::getXTractInputType() const
{
    //!!! I don't like this -- split into needXXX functions
    // (needTimeDomainAudio, needBarkCoefficients etc)

    switch (m_xtFeature) {

    case MEAN:
    case VARIANCE:
    case STANDARD_DEVIATION:
    case AVERAGE_DEVIATION:
    case SKEWNESS:
    case KURTOSIS:
    case IRREGULARITY_K:
    case IRREGULARITY_J:
    case SMOOTHNESS:
    case SPREAD:
    case FLATNESS:
    case CREST:
    case NOISINESS:
    case POWER:
    case SHARPNESS:
    case SLOPE:
    case HPS:
    case BARK_COEFFICIENTS:
    case PEAKS:
    case FLUX:
    case ROLLOFF:
    case MFCC:
        return MagnitudeSpectrum;

    case CENTROID:
    case TONALITY:
        return SpectralPeaks;

    case TRISTIMULUS_1:
    case TRISTIMULUS_2:
    case TRISTIMULUS_3:
    case ODD_EVEN_RATIO:
        return HarmonicSpectrum;

    case ZCR:
    case RMS_AMPLITUDE:
    case F0:
    case AUTOCORRELATION:
    case AMDF:
    case ASDF:
    case MAGNITUDE_SPECTRUM:
    case AUTOCORRELATION_FFT:
    case DCT:
    case ATTACK_TIME:
    case DECAY_TIME:
    case DELTA_FEATURE:
        return TimeDomainAudio;

    case LOUDNESS:
        return BarkCoefficients;

    case INHARMONICITY:
    case LOWEST_MATCH:
        return SpectralPeaksPlusF0;

    default:
        cerr << "ERROR: XTractPlugin::getXTractInputType: Unexpected feature index " << m_xtFeature << endl;
        return TimeDomainAudio;
    }
}

bool
XTractPlugin::needPeakThreshold() const
{
    switch (m_xtFeature) {
    case PEAKS:
    case CENTROID:
    case TONALITY:
    case INHARMONICITY:
    case LOWEST_MATCH:
        return true;

    default:
        return false;
    }
}

XTractPlugin::FeatureSet
XTractPlugin::process(const float *const *inputBuffers,
                      Vamp::RealTime timestamp)
{
    if (m_outputDescriptors.empty()) setupOutputDescriptors();

    int rbs = m_outputBinCount > m_blockSize ? m_outputBinCount : m_blockSize;
    if (!m_resultBuffer) {
        m_resultBuffer = new float[rbs];
    }
    for (int i = 0; i < rbs; ++i) m_resultBuffer[i] = 0.f;

    float *data = 0;
    int N = m_blockSize;
    void *argv = 0;

    FeatureSet fs;

    XTractInputType inputType = getXTractInputType();

    switch (inputType) {

    case TimeDomainAudio:
    case SpectralPeaksPlusF0:
        assert(getInputDomain() == TimeDomain);
        // libxtract functions may modify their input data, so we must copy
        data = new float[N];
        for (int n = 1; n < N; ++n) {
            data[n] = inputBuffers[0][n];
        }
        break;

    default:
        // All the rest are derived from the magnitude spectrum
        assert(getInputDomain() == FrequencyDomain);
        // Need same format as would be output by xtract_magnitude_spectrum
        data = new float[N];
        for (int n = 1; n < N/2; ++n) {
            data[n]   = sqrt(inputBuffers[0][n*2]   * inputBuffers[0][n*2] +
                             inputBuffers[0][n*2+1] * inputBuffers[0][n*2+1]) / N;
            data[N-n] = 0.f;
        }
        data[0]   = fabs(inputBuffers[0][0]) / N;
        data[N/2] = fabs(inputBuffers[0][N]) / N;
        break;
    }

    assert(m_outputBinCount > 0);

    float *result = m_resultBuffer;

    float argf[2];
    argv = &argf[0];

    float mean = 0.f, variance = 0.f, sd = 0.f;

    bool needSD = (m_xtFeature == SKEWNESS ||
                   m_xtFeature == KURTOSIS);

    bool needVariance = (needSD ||
                         m_xtFeature == STANDARD_DEVIATION);

    bool needMean = (needVariance ||
                     m_xtFeature == VARIANCE ||
                     m_xtFeature == AVERAGE_DEVIATION);

    bool needPeaks = (needPeakThreshold() && m_xtFeature != PEAKS);
    bool needBarkCoefficients = (m_xtFeature == LOUDNESS);

    float *barkCoefficients = 0;

    if (inputType == SpectralPeaksPlusF0) {
        if (processSPF0(data)) goto haveResult;
        else goto done;
    }

    if (needMean) {
        xtract_mean(data, N, 0, result);
        mean = *result;
        *result = 0.f;
    }

    if (needVariance) {
        argf[0] = mean;
        xtract_variance(data, N, argv, result);
        variance = *result;
        *result = 0.f;
    } 

    if (needSD) {
        argf[0] = variance;
        xtract_standard_deviation(data, N, argv, result);
        sd = *result;
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
    // CENTROID -- N/2 frequency and N/2 magnitude peaks
    // TONALITY -- peaks (different from those for CENTROID?)
    // TRISTIMULUS_1/2/3 -- harmonic spectrum
    // ODD_EVEN_RATIO -- harmonic spectrum
    // LOUDNESS -- Bark coefficients

    // argv should be now correct for all except:
    //
    // ROLLOFF -- threshold
    // F0 -- samplerate
    // MFCC -- Mel filter coefficients
    // BARK_COEFFICIENTS -- Bark band limits
    // PEAKS -- peak threshold and samplerate

    if (m_xtFeature == ROLLOFF) {
        argf[0] = m_threshold / 100.f;
        argv = &argf[0];
    }

    if (m_xtFeature == F0) {
        argf[0] = m_inputSampleRate;
        argv = &argf[0];
    }
    
    if (m_xtFeature == PEAKS || needPeaks) {
        argf[0] = m_threshold;
        argf[1] = m_inputSampleRate;
        argv = &argf[0];
    }

    if (m_xtFeature == BARK_COEFFICIENTS || needBarkCoefficients) {
        argv = &m_barkBandLimits[0];
    }

    xtract_mel_filter mfccFilterBank;
    if (m_xtFeature == MFCC) {
        mfccFilterBank.n_filters = m_coeffCount;
        mfccFilterBank.filters = m_mfccFilters;
        argv = &mfccFilterBank;
    }

    if (needPeaks) {
        int rv = xtract_peaks(data, N, argv, result);
        for (int n = 0; n < N; ++n) {
            data[n] = result[n];
            result[n] = 0.f;
        }
        // rv not trustworthy
//        if (rv != SUCCESS) {
//            cerr << "ERROR: XTractPlugin::process: xtract_peaks failed (error code = " << rv << ")" << endl;
//            goto done;
//        }
    }

    if (needBarkCoefficients) {
        barkCoefficients = new float[BARK_BANDS];
        int rv = xtract_bark_coefficients(data, N, argv, barkCoefficients);
//        if (rv != SUCCESS) {
//            cerr << "ERROR: XTractPlugin::process: xtract_bark_coefficients failed (error code = " << rv << ")" << endl;
//            goto done;
//        }
        data = barkCoefficients;
        N = BARK_BANDS;
        argv = 0;
    }

    if (inputType == HarmonicSpectrum) {

        //!!! Not yet implemented

    } else if (inputType == SpectralPeaks) {

        assert(needPeaks);
    }

    // now the main result

    xtract[m_xtFeature](data, N, argv, result);

haveResult:
    {
        int index = 0;

        for (size_t output = 0; output < m_outputDescriptors.size(); ++output) {

            Feature feature;
            feature.hasTimestamp = false;
            bool good = true;

            for (size_t n = 0; n < m_outputDescriptors[output].binCount; ++n) {
                float value = m_resultBuffer[index];
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
    }
   
done:
    if (barkCoefficients) delete[] barkCoefficients;
    delete[] data;

    cerr << "XTractPlugin::process returning" << endl;

    return fs;
}

bool
XTractPlugin::processSPF0(float *data)
{
    int N = m_blockSize;

    xtract_magnitude_spectrum(data, N, 0, m_resultBuffer);
    
    float *peaks = new float[N];
    float argf[2];
    argf[0] = m_threshold;
    argf[1] = m_inputSampleRate;
    int rv = xtract_peaks(m_resultBuffer, N, &argf[0], peaks);
    // rv not trustworthy
//    if (rv != SUCCESS) {
//        cerr << "ERROR: XTractPlugin::processSPF0: xtract_peaks failed (error code = " << rv << ")" << endl;
//        return false;
//    }

    float f0 = 0.f;
    rv = xtract_f0(data, N, &argf[1], &f0);
    if (rv != SUCCESS) {
        cerr << "ERROR: XTractPlugin::processSPF0: xtract_f0 failed (error code = " << rv << ")" << endl;
        return false;
    }

    float *mda[2];
    mda[0] = &f0;
    mda[1] = &(peaks[N/2]);

    for (int i = 0; i < N/2; ++i) m_resultBuffer[i] = 0.f;

    xtract[m_xtFeature](peaks, N/2, mda, m_resultBuffer);

    return true;
}


XTractPlugin::FeatureSet
XTractPlugin::getRemainingFeatures()
{
    return FeatureSet();
}

