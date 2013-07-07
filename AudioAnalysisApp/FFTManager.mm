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
    
    double window[nOver2];
    vDSP_hann_windowD(window, nOver2, vDSP_HANN_NORM);
//    vDSP_hamm_windowD(window, nOver2, 0);
    
    for (int i = 0; i < nOver2; i++) {
        splitComplex.realp[i] *= window[i];
//        splitComplex.imagp[i] *= window[i];
    }
    
    // For polar coordinates
    double *mag = new double[nOver2];
    double *phase = new double[nOver2];
    
    // Convert from complex/rectangular (real & imaginary) coordinates
    // to polar (magnitude & phase) coordinates
    
    double tmpData[nOver2];
    double mAdjust0BD = 1.5849e-13;
    double one = 1;
    
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

    // Convert from polar coords back to rectangular coords
    vDSP_ztocD(&splitComplex, 1, complex, 2, nOver2);
    vDSP_rectD((double *)complex, 2, (double *)complex, 2, nOver2); // Polar to rectangular conversion
    vDSP_ctozD(complex, 2, &splitComplex, 1, nOver2);

    // Computes an in-place single-precision real discrete Fourier transform,
    // from the frequency domain to the time domain (inverse).
    vDSP_fft_zripD(fftSetup, &splitComplex, 1, log2n, kFFTDirection_Inverse);

    // Unpack the result into a real vector
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
    
    [self detechPitchFor:tmpData withNumberOfFrames:nOver2 sampleRate:44100.0 lowBoundFreq:50.0 highBoundFreq:150.0];
    
    // Destroy the fft setup to avoid memory leakage
    vDSP_destroy_fftsetupD(fftSetup);
}

- (void)detechPitchFor:(double *)samples withNumberOfFrames:(int)n sampleRate:(double)sr lowBoundFreq:(double)lbf highBoundFreq:(double)hbf
{
    double freq = 0;
    double *hann = (double *)malloc(sizeof(double) * n);
    double *result = (double *)malloc(sizeof(double) * n);
    
    vDSP_hann_windowD(hann, n, vDSP_HANN_NORM);
    
    int returnIndex = 0;
    float sum;
    bool goingUp = false;
    float normalize = 0;
    
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < n; j++) {
            sum += (samples[j] * samples[j+i]) * hann[j];
        }
        if (i == 0) {
            normalize = sum;
        }
        result[i] = sum / normalize;
    }
    
    for (int i = 0; i < (n - 8); i++) {
        if (result[i] < 0) {
            i += 2; // no peaks below 0, skip forward at a faster rate
        } else {
            if(result[i] > result[i-1] && goingUp == false && i > 1) {
                //local min at i - 1
                goingUp = true;
            } else if(goingUp == true && result[i] < result[i-1]) {
                //local max at i-1
                if(returnIndex == 0 && result[i-1] > (result[0] * 0.95)) {
                    returnIndex = i - 1;
                    break;
                    //############### NOTE ##################################
                    // My implemenation breaks out of this loop when it finds the first peak.
                    // This is (probably) the greatest source of error, so if you would like to
                    // improve this algorithm, start here. the next else if() will trigger on
                    // future local maxima (if you first take out the break; above this paragraph)
                    //#######################################################
                } else if(result[i-1]>result[0]*0.85) {
                    // Do nothing
                }
                goingUp = false;
            }
        }
    }
    
    freq = sr / interpolate(result[returnIndex-1], result[returnIndex], result[returnIndex+1], returnIndex);
    NSLog(@"Pitch is %f", freq);
}

float interpolate(double y1, double y2, double y3, int k) {
    
    float d, kp;
    d = (y3 - y1) / (2 * (2 * y2 - y1 - y3));
    //printf("%f = %d + %f\n", k+d, k, d);
    kp  =  k + d;
    return kp;
}

@end