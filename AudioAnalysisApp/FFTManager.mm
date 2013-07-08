//
//  FFTManager.cpp
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#include "FFTManager.h"

const int LOG_N = 10;

@implementation FFTManager

- (void)performFFT:(SInt16 *)packets
{
    FFTSetup fftSetup = vDSP_create_fftsetup(LOG_N, kFFTRadix2);
    
    const int originalSize = sizeof(packets) / sizeof(SInt16) << LOG_N;
    
    // Buffer for real (time-domain) input and output signals
    float *x = new float[originalSize];
    float *y = new float[originalSize];
    
    // Initialize the input buffer with a sinusoid
    for (int k = 0; k < originalSize; k++) {
        x[k] = (float)packets[k];
    }
    
    // We need complex buffers for real and imaginary parts
    DSPComplex *tempComplex = new DSPComplex[originalSize / 2];
    
    DSPSplitComplex tempSplitComplex;
    tempSplitComplex.realp = new float[originalSize / 2];
    tempSplitComplex.imagp = new float[originalSize / 2];
    
    // For polar coordinates
    float *mag = new float[originalSize / 2];
    float *mag_Db = new float[originalSize / 2];
    float *phase = new float[originalSize / 2];
    
    // Scramble-pack the real data into the complex buffer in the way
    // that it's required by the real-to-complex FFT function.
    vDSP_ctoz((DSPComplex *)x, 2, &tempSplitComplex, 1, originalSize / 2);
    
    // real->complex forward FFT
    vDSP_fft_zrip(fftSetup, &tempSplitComplex, 1, LOG_N, kFFTDirection_Forward);
    
    printf("\nSpectrum:\n");
    for (int k = 0; k < originalSize / 2; k++) {
        printf("%3d\t%6.2f\t%6.2f\n", k, tempSplitComplex.realp[k], tempSplitComplex.imagp[k]);
    }
    
    // Convert from complex/rectangular (real & imaginary) coordinates
    // to polar (magnitude & phase) coordinates
    float mic_scale = (originalSize) * 2;
    vDSP_zvabs(&tempSplitComplex, 1, mag, 1, originalSize / 2);
    vDSP_zvphas(&tempSplitComplex, 1, phase, 1, originalSize / 2);
    vDSP_vdbcon(mag, 1, &mic_scale, mag_Db, 1, originalSize/2, 0);
    
    printf("\nMag / Phase:\n");
    for (int k = 0; k < originalSize / 2; k++) {
        printf("%3d\t%6.2f\t%6.2f\n", k, mag[k], phase[k]);
    }
    float freqStep =  22000/ (originalSize / 2);
    printf("\nDecibels:\n");
    for (int k = 0; k < originalSize / 2; k++) {
        printf("%f\t%f\n", (k * freqStep), mag_Db[k]);
    }
    
    // Convert from polar coords back to rectangular coords
    tempSplitComplex.realp = mag;
    tempSplitComplex.imagp = phase;
    
    vDSP_ztoc(&tempSplitComplex, 1, tempComplex, 2, originalSize / 2);
    vDSP_rect((float *)tempComplex, 2, (float *)tempComplex, 2, originalSize / 2);
    vDSP_ctoz(tempComplex, 2, &tempSplitComplex, 1, originalSize / 2);
    
    // Do Inverse FFT
    
    // complex->real inverse FFT
    vDSP_fft_zrip(fftSetup, &tempSplitComplex, 1, LOG_N, kFFTDirection_Inverse);
    
    // Unpack the result into a real vector
    vDSP_ztoc(&tempSplitComplex, 1, (DSPComplex *)y, 2, originalSize / 2);
    
    // Compensate for scaling for both FFTs
    float scale = 1.0 / originalSize;
    vDSP_vsmul(y, 1, &scale, y, 1, originalSize);
    
    // For a sinusodial wave, input x and input y vectors will have identical values
    printf("\nInput & output:\n");
    for (int k = 0; k < originalSize / 2; k++) {
        printf("%3d\t%6.2f\t%6.2f\n", k, x[k], y[k]);
    }
}

@end