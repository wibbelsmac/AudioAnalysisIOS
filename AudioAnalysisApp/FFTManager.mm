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

#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

@interface FFTManager(Private)

- (void)plotRange:(int)range forPackets:(double *)packets;
- (void)plotAllValuesStartingPacket:(int)sp EndWithPacket:(int)ep forPackets:(double *)packets;
- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev;

@end

@implementation FFTManager

const int sampleFreq = 44100;
const float sampleFrameTime = .1f;
const int numSamplesPerFrame = sampleFrameTime * sampleFreq;
const float subFrameFract = .60f;
const int numSamplesPerSubFrame = numSamplesPerFrame * subFrameFract;
const int diffSubandFrame = numSamplesPerFrame - numSamplesPerSubFrame;

- (id)initWithPackets:(SInt16 *)packets PacketCount:(int)packetCount
{
    self = [super init];
    if (self) {
        _mLog2N = ceil(log2(packetCount));
        _mNumberOfPackets = (int) pow(2, _mLog2N);
        _mFFTLength = _mNumberOfPackets / 2;
        NSLog(@"n = %d", _mNumberOfPackets);
        
        // An opaque type that contains setup information for a given double-precision FFT transform.
        _fftSetup = vDSP_create_fftsetupD(_mLog2N, kFFTRadix2);
        
        // Initialize the input buffer with packets from the audio file
        _mAudioBuffer = (double *) calloc(_mNumberOfPackets, sizeof(double));
        for (int k = 0; k < _mNumberOfPackets; k++) {
            _mAudioBuffer[k] = (double) packets[k];
        }
        free(packets);
    }
    return self;
}

- (void)performFFTOnRangeStartingAt:(int)sp EndingAt:(int)ep
{
    _mLog2N = ceil(log2(ep - sp));
    _mNumberOfPackets = (int) pow(2, _mLog2N);
    _mFFTLength = _mNumberOfPackets / 2;
    
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    _mSplitComplex.realp = (double *) calloc(_mFFTLength, sizeof(double));
    _mSplitComplex.imagp = (double *) calloc(_mFFTLength, sizeof(double));
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)_mAudioBuffer, 2, &_mSplitComplex, 1, _mFFTLength);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(self.fftSetup, &_mSplitComplex, 1, _mLog2N, kFFTDirection_Forward);
    
    [self plotRange:ep - sp forPackets:_mAudioBuffer];
    
    double scale = 1.0 / (2 * _mNumberOfPackets);
    vDSP_vsmulD(_mSplitComplex.realp, 1, &scale, _mSplitComplex.realp, 1, _mFFTLength);
    vDSP_vsmulD(_mSplitComplex.imagp, 1, &scale, _mSplitComplex.imagp, 1, _mFFTLength);
    
    _mSplitComplex.realp[0] = 0; // zeroing out the DC value
    _mSplitComplex.imagp[0] = 0; // zeroing out the Nyquist value
    
    double tmpData[_mFFTLength];
    double mAdjust0BD = 1.5849e-13;
    double one = 1.0;
    
    // Convert from complex/rectangular (real & imaginary) coordinates
    // to polar (magnitude & phase) coordinates
    vDSP_zvmagsD(&_mSplitComplex, 1, tmpData, 1, _mFFTLength);
    vDSP_vsaddD(tmpData, 1, &mAdjust0BD, tmpData, 1, _mFFTLength);
    vDSP_vdbconD(tmpData, 1, &one, tmpData, 1, _mFFTLength, 0);
    
    for (int i = 0; i < _mFFTLength; i++) {
        if (tmpData[i] < 0.0) {
            tmpData[i] = 0;
        }
    }
    
//    [self plotRange:_mFFTLength forPackets:tmpData];
    
    DSPDoubleComplex *complex = new DSPDoubleComplex[_mFFTLength];
    
    // For polar coordinates
    double *mag = new double[_mFFTLength];
    double *phase = new double[_mFFTLength];
    
    // Get complex vector absolute values
    vDSP_zvabsD(&_mSplitComplex, 1, mag, 1, _mFFTLength);
    // Get complex vector phase
    vDSP_zvphasD(&_mSplitComplex, 1, phase, 1, _mFFTLength);
    
    _mSplitComplex.realp = mag;
    _mSplitComplex.imagp = phase;
    
    // Convert from polar coords back to rectangular coords
    vDSP_ztocD(&_mSplitComplex, 1, complex, 2, _mFFTLength);
    vDSP_rectD((double *)complex, 2, (double *)complex, 2, _mFFTLength); // Polar to rectangular conversion
    vDSP_ctozD(complex, 2, &_mSplitComplex, 1, _mFFTLength);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the frequency domain to the time domain (inverse).
    vDSP_fft_zripD(_fftSetup, &_mSplitComplex, 1, _mLog2N, kFFTDirection_Inverse);
    
    // Unpack the result into a real vector
    vDSP_ztocD(&_mSplitComplex, 1, (DSPDoubleComplex *)_mAudioBuffer, 2, _mLog2N);
    
    // Compensate for scaling for both FFTs. See Apple's vDSP documentation
    // to see how much scaling is required.
    // NOTE: Apple scales values to increase the efficiency of their functions.
    // As such, we have to scale them back.
    vDSP_vsmulD(_mAudioBuffer, 1, &scale, _mAudioBuffer, 1, _mLog2N);
    
//    [self plotRange:ep - sp forPackets:_mAudioBuffer];
    
    free(_mSplitComplex.realp);
    free(_mSplitComplex.imagp);
}

- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev
{
    int range = ev - sv;
    double *window = new double[range];
    vDSP_hann_windowD(window, self.mNumberOfPackets, vDSP_HANN_NORM);
    
    for (int i = 0, index = sv; index < ev; i++, index++) {
        values[index] *= window[i];
    }
}

- (void)dealloc
{
    vDSP_destroy_fftsetupD(_fftSetup);
    free(_mAudioBuffer);
}

- (void)performFFT:(SInt16 *)packets PacketCount:(int)packetCount
{
    
    // initialize butter worth filter with 1k cuttoff
//    AAAButterLow *bLow = [[AAAButterLow alloc] initWithSampleFreq:44100.0f CutoffFreq:1000.0f];
//    int dataLength = byteCount/ sizeof(packets);
//    
//    float* floatData = (float*) malloc(dataLength * sizeof(float));
//    vDSP_vflt16(packets, 1, floatData, 1, dataLength);
//    float* filteredResult = (float*) malloc(dataLength * sizeof(float));
//
//    [bLow filterArray:floatData DataLength:dataLength ResultArray:filteredResult ResultLength:dataLength];
//    [self performCrossCorrelation:filteredResult NumSamples: dataLength];
    
    // ------------------------------------------------------------------------
}

- (void)plotRange:(int)range forPackets:(double *)packets
{
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    int resolution = (int) ceil(log2(range));
    resolution = CLAMP(resolution, 10, 20);
    for (int i = 0, index = 0; index < range; i++, index += 128) {
        xValues[i] = [NSNumber numberWithInt:index];
        yValues[i] = [NSNumber numberWithDouble:packets[index]];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
}

- (void)plotAllValuesStartingPacket:(int)sp EndWithPacket:(int)ep forPackets:(double *)packets
{
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    for (int i = 0, index = sp; index < ep; i++, index++) {
        xValues[i] = [NSNumber numberWithInt:index];
        yValues[i] = [NSNumber numberWithDouble:packets[index]];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
}

- (void)performCrossCorrelation:(float *)packets NumSamples:(int) numSamples
{    
    const int sampleCount = numSamples;
    // Populate *window with the values for a hamming window function
    float *hannWindow = (float *)malloc(sizeof(float) * numSamplesPerFrame);
    float *tempSamples = (float *)malloc(sizeof(float) * numSamplesPerFrame);
    float *subFrameSamples = (float *)malloc(sizeof(float) * numSamplesPerSubFrame);
    float *resultWindow = (float *)malloc(sizeof(float) * diffSubandFrame);
    
    vDSP_hamm_window(hannWindow, numSamplesPerFrame, 0);
    
    NSMutableArray *peakValues = [[NSMutableArray alloc] init];
    
    for (int i = 0; i < sampleCount - numSamplesPerFrame; i += numSamplesPerFrame) {
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

- (void)performVoiceAnalysisOn:(double *)frames
{
    const int p = 2048;
    const int lenResult = 2 * p - 1;
    double *resultsArray = new double[lenResult];

    vDSP_convD(frames, 1, frames, 1, resultsArray, 1, lenResult, p); // Auto Correlation
    NSLog(@"Autocorrelation:\n");
    for (int i = 0; i < lenResult; i++) {
        NSLog(@"%d = %.3f", i, resultsArray[i]);
    }
    
    // For graphing values
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    for (int i = 0; i < p; i++) {
        xValues[i] = [NSNumber numberWithInt:i];
        yValues[i] = [NSNumber numberWithDouble:resultsArray[i]];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
    
    
    long int N = p, NRHS = p, LDA = p, LDB = p, INFO;
    long int IPIV[p];
    // N -> Number of linear equations, i.e. number of rows for a square matrix A
    // NRHS -> Number of right-hand sides, i.e. number of columns for a matrix B
    // A -> Coefficient matrix A
    // LDA -> Dimension of array A, LDA >= max(1, N)
    // IPIV -> Integer array with dimension N. Pivot indices that define the permutation matrix P.
    // B -> The array B with dimension LDB, NRHS
    // LDB -> Leading dimension of array B, LDB >= max(1, N)
    // INFO -> Status
    //          = 0:  successful exit
    //          < 0:  if INFO = -i, the i-th argument had an illegal value
    //          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    //                has been completed, but the factor U is exactly
    //                singular, so the solution could not be computed.
//    dgesv_(&N, &NRHS, coeffArray, &LDA, IPIV, corrResults, &LDB, &INFO);
    
//    NSLog(@"Status: %ld", INFO);
    
//    NSLog(@"After:\n");
//    for (int a = 0; a < p; a++) {
//        NSLog(@"%.3f ", corrResults[a]);
//    }
}

@end