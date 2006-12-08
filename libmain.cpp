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

#include <vamp/vamp.h>
#include <vamp-sdk/PluginAdapter.h>

#include "plugins/XTractPlugin.h"

#include <xtract/libxtract.h>

#include <map>
#include <cassert>

class XTractPluginAdapter : public Vamp::PluginAdapterBase
{
public:
    XTractPluginAdapter(unsigned int feature) :
        PluginAdapterBase(),
        m_xtFeature(feature) 
    { }

    virtual ~XTractPluginAdapter() { }

protected:
    Vamp::Plugin *createPlugin(float inputSampleRate) {
        XTractPlugin *xtp = new XTractPlugin(m_xtFeature, inputSampleRate);
        return xtp;
    }

    unsigned int m_xtFeature;
};

static std::map<unsigned int, XTractPluginAdapter *> pluginAdapterMap;

// Define this if libxtract has been compiled with XTRACT_FFT and
// linked with fftw3
#define HAVE_XTRACT_FFT 1

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int index)
{
    // The missingFeatures array contains the libxtract enum values
    // for features that we are not providing in this plugin, either
    // because they aren't very meaningful in isolation or because the
    // version of libxtract we're using doesn't support them yet.
    // 
    // These must appear in this array in the same order as in the
    // original enum in libxtract.h.

    const unsigned int missingFeatures[] = {

        TRISTIMULUS_1,   // not equipped to provide harmonic spectrum input
        TRISTIMULUS_2,   // not equipped to provide harmonic spectrum input
        TRISTIMULUS_3,   // not equipped to provide harmonic spectrum input
        TONALITY,        // nonsensical definition in scalar.c?
        CREST,           // not implemented
        NOISINESS,       // not implemented
        POWER,           // not implemented
        ODD_EVEN_RATIO,  // not equipped to provide harmonic spectrum input
        SHARPNESS,       // not implemented
        SLOPE,           // not implemented
        LOWEST_MATCH,    // not meaningful?
        HPS,             // "this function doesn't work properly"
        FLUX,            // not implemented
        ATTACK_TIME,     // not implemented
        DECAY_TIME,      // not implemented
        DELTA_FEATURE,   // not implemented (and not meaningful?)

#ifndef HAVE_XTRACT_FFT
        MAGNITUDE_SPECTRUM,  // requires fftw
        AUTOCORRELATION_FFT, // requires fftw
        MFCC,                // requires fftw
        DCT,                 // requires fftw
#endif
    };

    for (unsigned int i = 0;
         i < sizeof(missingFeatures)/sizeof(missingFeatures[0]);
         ++i) {
        if (i > 0) assert(missingFeatures[i] > missingFeatures[i-1]);
        if (index >= missingFeatures[i]) ++index;
    }

    if (index >= XTRACT_FEATURES) return 0;

    if (pluginAdapterMap.find(index) == pluginAdapterMap.end()) {
        pluginAdapterMap[index] = new XTractPluginAdapter(index);
    }

    return pluginAdapterMap[index]->getDescriptor();
}

