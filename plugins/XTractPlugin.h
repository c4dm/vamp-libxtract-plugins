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

#ifndef _XTRACT_PLUGIN_H_
#define _XTRACT_PLUGIN_H_

#include <vamp-sdk/Plugin.h>

class XTractPlugin : public Vamp::Plugin
{
public:
    XTractPlugin(unsigned int xtFeature, float inputSampleRate);
    virtual ~XTractPlugin();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const;

    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string) const;
    void setParameter(std::string, float);

    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    enum XTractInputType {
        TimeDomainAudio,
        MagnitudeSpectrum,
        SpectralPeaks,
        SpectralPeaksPlusF0,
        HarmonicSpectrum,
        BarkCoefficients,
        CustomInputType
    };
    XTractInputType getXTractInputType() const;
    bool needPeakThreshold() const;

    mutable OutputList m_outputDescriptors;
    void setupOutputDescriptors() const;

    bool processSPF0(float *data);

    const unsigned int m_xtFeature;
    size_t m_channels;
    size_t m_stepSize;
    size_t m_blockSize;

    float *m_resultBuffer;

    float m_threshold;

    float m_minFreq;
    float m_maxFreq;

    size_t m_coeffCount;
    float **m_mfccFilters;
    int m_mfccStyle;

    int *m_barkBandLimits;

    size_t m_outputBinCount;
    bool m_initialised;
};


#endif
