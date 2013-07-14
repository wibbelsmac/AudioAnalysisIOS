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
#include "AAAButterLow.h"

#import "GSGraph.h"

@implementation FFTManager

const int sampleFreq = 44100;
const float sampleFrameTime = .1f;
const int numSamplesPerFrame = sampleFrameTime * sampleFreq;
const float subFrameFract = .60f;
const int numSamplesPerSubFrame = numSamplesPerFrame * subFrameFract;
const int diffSubandFrame = numSamplesPerFrame - numSamplesPerSubFrame;

- (void)performFFT:(SInt16 *)packets ByteCount:(int) byteCount
{
    
    // initialize butter worth filter with 1k cuttoff
    AAAButterLow *bLow = [[AAAButterLow alloc] initWithSampleFreq:44100.0f CutoffFreq:1000.0f];
    int dataLength = byteCount/ sizeof(packets);
    
    float* floatData = (float*) malloc(dataLength * sizeof(float));
    vDSP_vflt16(packets, 1, floatData, 1, dataLength);
    float* filteredResult = (float*) malloc(dataLength * sizeof(float));

    [bLow filterArray:floatData DataLength:dataLength ResultArray:filteredResult ResultLength:dataLength];
    [self performCrossCorrelation:filteredResult NumSamples: dataLength];
    const int log2n = 13;
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
    
    for (int i = 0; i < nOver2; i++) {
        splitComplex.realp[i] *= window[i];
    }

    
    // For polar coordinates
    double *mag = new double[nOver2];
    double *phase = new double[nOver2];
    
    double tmpData[nOver2];
    double mAdjust0BD = 1.5849e-13;
    double one = 1;
    
    // Convert from complex/rectangular (real & imaginary) coordinates
    // to polar (magnitude & phase) coordinates
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
//    NSMutableArray *xValues = [[NSMutableArray alloc] init];
//    NSMutableArray *yValues = [[NSMutableArray alloc] init];
//    for (int i = 0; i < nOver2 / 2; i++) {
//        xValues[i] = [NSNumber numberWithInt:i];
//        yValues[i] = [NSNumber numberWithDouble:tmpData[i]];
////        printf("%d %f\n", i, tmpData[i]);
//    }
//    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
//    [graph drawAxisLines];
//    [graph plotXAndYValues];

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
    
    // Destroy the fft setup to avoid memory leakage
    vDSP_destroy_fftsetupD(fftSetup);
}
- (void)performCrossCorrelation:(float *)packets NumSamples:(int) numSamples {
    
    
    const int sampleCount = numSamples;
    // Populate *window with the values for a hamming window function
    float *hannWindow = (float *)malloc(sizeof(float) * numSamplesPerFrame);
    float *tempSamples = (float *)malloc(sizeof(float) * numSamplesPerFrame);
    float *subFrameSamples = (float *)malloc(sizeof(float) * numSamplesPerSubFrame);
    float *resultWindow = (float *)malloc(sizeof(float) * diffSubandFrame);
    
    vDSP_hamm_window(hannWindow, numSamplesPerFrame, 0);
    
    NSMutableArray *peakValues = [[NSMutableArray alloc] init];
    
    for(int i = 0; i < sampleCount - numSamplesPerFrame; i += numSamplesPerFrame) {
        // Window the samples
        int j = 0;
        for(; j < numSamplesPerSubFrame; j++) {
            tempSamples[j] = packets[i + j] * hannWindow[j];
            subFrameSamples[j] = packets[i + j] * hannWindow[j];
        }
        for(; j < numSamplesPerFrame; j++) {
            tempSamples[j] = packets[i + j] * hannWindow[j];
        }
        vDSP_conv(tempSamples,1, subFrameSamples,1,resultWindow, 1, diffSubandFrame, numSamplesPerSubFrame);
        [peakValues addObject:[NSNumber numberWithInt:[self findPitchPeakFreq:resultWindow NumFrames:diffSubandFrame]]];
    }
    
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    for (int i = 0; i < peakValues.count; i++) {
        xValues[i] = [NSNumber numberWithInt:i];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:peakValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
   
}

-(int) findPitchPeakFreq:(float*)result NumFrames:(int) numFrames {
    int peakArray[2] = {-1,-1};
    int i = 0;
    bool goingUp = true;
    for(int peak = 0; peak < 2; peak++) {
        for(; i < numFrames; i++)
        {
            if(result[i]<0) {
                //i+=2; // no peaks below 0, skip forward at a faster rate
            } else {
                if(result[i]>result[i-1] && goingUp == false && i >1) {
                
                    //local min at i-1
                    goingUp = true;
                
                } else if(goingUp == true && result[i]<result[i-1]) {
                
                    //local max at i-1
                    if(peakArray[peak] == 0 && result[i-1]>result[i]*0.95) {
                        
                    } else if(result[i-1]>result[0]*0.85) {
                    }
                    goingUp = false;
                    peakArray[peak] = i-1;
                    break;
                }
            }
        }
    }
    //NSLog(@"%d", sampleFreq/(peakArray[1] - peakArray[0]));
    if(peakArray[0] > -1 && peakArray[1] > -1) {
    printf("p0: %d, p1:%d, %d Hz \n", peakArray[0], peakArray[1], sampleFreq/(peakArray[1] - peakArray[0]));
    return sampleFreq / (peakArray[1] - peakArray[0]);
    } else {
        if((peakArray[0] > -1)) {
            printf("Could not find any peak\n");
        } else {
            printf("Could not find second peak\n");
        }
        return 0;
    }
}



@end