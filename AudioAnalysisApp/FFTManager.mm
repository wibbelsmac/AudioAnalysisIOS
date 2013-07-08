//
//  FFTManager.cpp
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "FFTManager.h"

#import "GSGraph.h"

@implementation FFTManager

- (void)performFFT:(SInt16 *)packets
{    
    const int log2n = 12;
    const int n = (sizeof(packets) / sizeof(SInt16)) << log2n; // multiply size by 2^log2n
    const int nOver2 = n / 2;
    NSLog(@"%ld %ld %d", sizeof(packets), sizeof(SInt16), n);
    
    // An opaque type that contains setup information for a given double-precision FFT transform.
    FFTSetupD fftSetup = vDSP_create_fftsetupD(log2n, kFFTRadix2);

    // Buffers for input signals before processing and output signals after processing
    double *input = new double[n];
    double *output = new double[n];
    
    // Initialize the input buffer with packets from the audio file
    for (int k = 0; k < n; k++) {
        input[k] = (double)packets[k];
    }
    
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    DSPDoubleComplex *complex = new DSPDoubleComplex[nOver2];
    DSPDoubleSplitComplex splitComplex;
    splitComplex.realp = new double[nOver2];
    splitComplex.imagp = new double[nOver2];
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)input, 2, &splitComplex, 1, nOver2);

    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(fftSetup, &splitComplex, 1, log2n, kFFTDirection_Forward);
    
<<<<<<< HEAD
    // For polar coordinates
    float *mag = new float[originalSize / 2];
    float *mag_Db = new float[originalSize / 2];
    float *phase = new float[originalSize / 2];
    
    // Scramble-pack the real data into the complex buffer in the way
    // that it's required by the real-to-complex FFT function.
    vDSP_ctoz((DSPComplex *)x, 2, &tempSplitComplex, 1, originalSize / 2);
=======
    double window[nOver2];
    vDSP_hann_windowD(window, nOver2, vDSP_HANN_NORM);
    
    for (int i = 0; i < nOver2; i++) {
        splitComplex.realp[i] *= window[i];
    }
>>>>>>> a40a4d3bcfe58dedd0fbf56ddb09beca808ed436
    
    // For polar coordinates
    double *mag = new double[nOver2];
    double *phase = new double[nOver2];
    
    double tmpData[nOver2];
    double mAdjust0BD = 1.5849e-13;
    double one = 1;
    
    // Convert from complex/rectangular (real & imaginary) coordinates
    // to polar (magnitude & phase) coordinates
<<<<<<< HEAD
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
    
=======
    vDSP_zvmagsD(&splitComplex, 1, tmpData, 1, nOver2);
    vDSP_vsaddD(tmpData, 1, &mAdjust0BD, tmpData, 1, nOver2);
    vDSP_vdbconD(tmpData, 1, &one, tmpData, 1, nOver2, 0);
    
    // Get complex vector absolute values
    vDSP_zvabsD(&splitComplex, 1, mag, 1, nOver2);
    // Get complex vector phase
    vDSP_zvphasD(&splitComplex, 1, phase, 1, nOver2);
    
    splitComplex.realp = mag;
    splitComplex.imagp = phase;
    
    // For graphing values
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    for (int i = 0; i < nOver2 / 2; i++) {
        xValues[i] = [NSNumber numberWithInt:i];
        yValues[i] = [NSNumber numberWithDouble:tmpData[i]];
//        printf("%d %f\n", i, tmpData[i]);
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];

>>>>>>> a40a4d3bcfe58dedd0fbf56ddb09beca808ed436
    // Convert from polar coords back to rectangular coords
    vDSP_ztocD(&splitComplex, 1, complex, 2, nOver2);
    vDSP_rectD((double *)complex, 2, (double *)complex, 2, nOver2); // Polar to rectangular conversion
    vDSP_ctozD(complex, 2, &splitComplex, 1, nOver2);

    // Computes an in-place single-precision real discrete Fourier transform,
    // from the frequency domain to the time domain (inverse).
    vDSP_fft_zripD(fftSetup, &splitComplex, 1, log2n, kFFTDirection_Inverse);

    // Unpack the result into a real vector
<<<<<<< HEAD
    vDSP_ztoc(&tempSplitComplex, 1, (DSPComplex *)y, 2, originalSize / 2);
    
    // Compensate for scaling for both FFTs
    float scale = 1.0 / originalSize;
    vDSP_vsmul(y, 1, &scale, y, 1, originalSize);
    
    // For a sinusodial wave, input x and input y vectors will have identical values
    printf("\nInput & output:\n");
    for (int k = 0; k < originalSize / 2; k++) {
        printf("%3d\t%6.2f\t%6.2f\n", k, x[k], y[k]);
    }
=======
    vDSP_ztocD(&splitComplex, 1, (DSPDoubleComplex *)output, 2, nOver2);

    // Compensate for scaling for both FFTs. See Apple's vDSP documentation
    // to see how much scaling is required.
    // NOTE: Apple scales values to increase the efficiency of their functions.
    // As such, we have to scale them back.
    double scale = 1.0 / n;
    vDSP_vsmulD(output, 1, &scale, output, 1, nOver2);
    
//    // For graphing values
//    NSMutableArray *xValues = [[NSMutableArray alloc] init];
//    NSMutableArray *yValues = [[NSMutableArray alloc] init];
//    for (int i = 0; i < nOver2; i++) {
//        xValues[i] = [NSNumber numberWithInt:i];
//        yValues[i] = [NSNumber numberWithDouble:output[i]];
//        printf("%d %f\n", i, output[i]);
//    }
//    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
//    [graph drawAxisLines];
//    [graph plotXAndYValues];
    
    // Destroy the fft setup to avoid memory leakage
    vDSP_destroy_fftsetupD(fftSetup);
>>>>>>> a40a4d3bcfe58dedd0fbf56ddb09beca808ed436
}

@end