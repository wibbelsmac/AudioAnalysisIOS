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
#include <vector>
#include <stdlib.h>
#include "FFTManager.h"
#include "AAAButterLow.h"
#include <vector>

#import <Accelerate/Accelerate.h>
#import "GSGraph.h"

#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

@interface FFTManager(Private)

- (void)plotValues:(double *)values inRange:(int)range;
- (void)plotValues:(double *)values startingAt:(int)start endingAt:(int)end;

- (void)calculateVowelOnsetPointForSamples:(double *)samples withLength:(int)length;
- (void)applyHanningWindowOnSamples:(double *)samples startingAt:(int)start endingAt:(int)end;
- (void)performForwardFFTOnSamples:(double **)samples withLength:(int *)length;
- (void)performInverseFFTOnSamples:(double **)samples withLength:(int *)length;
- (void)performHilbertPhaseShiftsOnSamples:(double *)samples withLength:(int)length;
- (void)calculateHilbertEnvelopeForSamples:(double *)samples andTransformArray:(double *)transformArray withLength:(int)length;
- (void)convolveSamples:(double *)samples withLength:(int)samplesLength withFilter:(double *)filter withFilterLength:(int)filterLength
           putResultsIn:(double *)results forResultLength:(int)resultLength;

void ForwardLinearPrediction(std::vector<double> &coeffs, const std::vector<double> &x);

@end

@implementation FFTManager

const int sampleFreq = 44100;
const float sampleFrameTime = .02f;
int numSamplesPerFrame = sampleFrameTime * sampleFreq;
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
            _mAudioBuffer[k] = (double) (packets[k]);
        }
        
        [self calculateVowelOnsetPointForSamples:_mAudioBuffer withLength:_mNumberOfPackets];
    }
    free(packets);
    return self;
}

- (void)calculateVowelOnsetPointForSamples:(double *)samples withLength:(int)length {
    double hanning_window_time = .02;
    double overlap_scale = .5;
    int hanning_window_length =  hanning_window_time * sampleFreq;  // hamming window length
    
    // calculating linear prediction coefficients
    int total_window_size = pow(2, ceil(log2(hanning_window_length)));
    int lpc_order = 10;

    int number_overlap_windows = ((length - hanning_window_length) / hanning_window_length) /overlap_scale;
    int total_samples_processed = number_overlap_windows * hanning_window_length / 2;
    
    double* pWindow = (double*) calloc(total_window_size, sizeof(double));
    double* pLP_residual = (double *) calloc(total_window_size, sizeof(double));
    double* pHilbert_transform_LP_residual = (double *) calloc(total_window_size, sizeof(double));
    double* lpc_coefficient = (double*) malloc(sizeof(double) * lpc_order);
    double* result_data = (double*) calloc(length, sizeof(double));
    double* convolved_result_data = (double*) calloc(length, sizeof(double));

    for(int i = 0; i < number_overlap_windows; i ++) {
        // copy data over into temp window
        int window_start = i * hanning_window_length / 2;
        memcpy(pWindow, (samples + window_start), hanning_window_length * sizeof(double));
        std::vector<double> coeffs (lpc_order, 0.0);
        std::vector<double> dataVector (pWindow, pWindow + hanning_window_length);

        //  1.  calculate the LP Coefficients
        ForwardLinearPrediction(coeffs, dataVector);
        memcpy(lpc_coefficient, &(coeffs[0]),lpc_order * sizeof(double));
        
        //  2.  calculate the LP Residual
        [self inverseFilter:pWindow Result:pLP_residual FilterCoefficients:lpc_coefficient // construct LPC residual
                ArrayLength:hanning_window_length FilterLength:lpc_order];
        //  3. copy data over and perform hilbert transform
        memcpy(pHilbert_transform_LP_residual, pLP_residual, total_window_size * sizeof(double)); // copy lpc residual for Hilbert Trans of LPC Resid
        [self applyHanningWindowOnSamples:pWindow startingAt:0 endingAt:hanning_window_length];
        [self performForwardFFTOnSamples:&pHilbert_transform_LP_residual withLength:&total_window_size]; // forward fft
        [self performHilbertPhaseShiftsOnSamples:pHilbert_transform_LP_residual withLength:total_window_size]; // hilbert transform phase shifts
        [self performInverseFFTOnSamples:&pHilbert_transform_LP_residual withLength:&total_window_size]; // complete hilbert trans
        [self calculateHilbertEnvelopeForSamples:pLP_residual andTransformArray:pHilbert_transform_LP_residual withLength:hanning_window_length]; // envelope both values toghether
        
        for(int z = 0; z < hanning_window_length; z++) {  // copy over data to result
            if(z > hanning_window_length / 2) {
               result_data[window_start + z] +=  pHilbert_transform_LP_residual[z];  // scale to .5 for overlap
            } else {
                result_data[window_start + z] = pHilbert_transform_LP_residual[z];
            }
        }
    }
    
    // 50ms hanning window to smooth data after
    int convWindowLength = .05 * sampleFreq;
    double* hWindow = (double*) calloc(convWindowLength, sizeof(double)); // construct hamming window
    vDSP_hamm_windowD(hWindow, convWindowLength, 0);
//    [self convolveData:result_data DataLength:total_samples_processed Filter:hWindow FilterLength:convWindowLength
//          Result:convolved_result_data  ResultLength:(total_samples_processed - hanning_window_length)];
    
    [self scaleData:samples Length:length Max:1.0 Min:-1.0];
    [self scaleData:result_data Length:length Max:1.0 Min:-1.0];
    [self scaleData:convolved_result_data Length:length Max:1.0 Min:-1.0];
    for(int j = 0; j < (total_samples_processed); j++) {
        printf("%d\t%f\t%f\n", j, result_data[j], samples[j]);
    }
    
    [self plotValues:result_data startingAt:0 endingAt:6000];

}

- (void)applyHanningWindowOnSamples:(double *)samples startingAt:(int)start endingAt:(int)end {
    int range = end - start;
    double *window = new double[range];
    vDSP_hann_windowD(window, range, vDSP_HANN_NORM);
    
    for (int i = 0, index = start; index < end; i++, index++) {
        samples[index] *= window[i];
    }
}

- (void)dealloc {
    vDSP_destroy_fftsetupD(_fftSetup);
    free(_mAudioBuffer);
}

- (void)plotValues:(double *)values inRange:(int)range {
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    int resolution = (int) ceil(log2(range));
    resolution = CLAMP(resolution, 10, 20);
    for (int i = 0, index = 0; index < range; i++, index += resolution) {
        xValues[i] = [NSNumber numberWithInt:index];
        yValues[i] = [NSNumber numberWithDouble:values[index]];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
}

- (void)plotValues:(double *)values startingAt:(int)start endingAt:(int)end {
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    for (int i = 0, index = start; index < end; i++, index++) {
        xValues[i] = [NSNumber numberWithInt:index];
        yValues[i] = [NSNumber numberWithDouble:values[index]];
    }
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
}

-(int) findPitchPeakFreq:(double*)result NumFrames:(int) numFrames {
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

// performs inverse filter result length must equal arraylength
// modeled after this page http://iitg.vlab.co.in/?sub=59&brch=164&sim=616&cnt=1109
-(void) inverseFilter:(double*) data Result:(double*)result FilterCoefficients:(double*)filter
        ArrayLength:(int) al FilterLength:(int)fl {
    int inverseFilterLength = fl + 1;
    static double * temp = (double*) calloc(inverseFilterLength, sizeof(double));
    temp[0] = 1;
    memccpy(temp + 1, filter, fl, sizeof(double));
    // copy the last filter length of values 
    //memccpy(result + (al - inverseFilterLength), data + (al - inverseFilterLength), inverseFilterLength, sizeof(double));
    // the data vector must be atlest resultLength + filterLength
    [self convolveSamples:data withLength:al withFilter:filter withFilterLength:inverseFilterLength putResultsIn:result forResultLength:al];
}

- (void)performForwardFFTOnSamples:(double **)samples withLength:(int *)length {
    unsigned long log2N = ceil(log2(*length));
    int pow2length = (int) pow(2, log2N);
    if(pow2length != *length) {
        double* tempStorage = (double*) calloc(pow2length, sizeof(double));
        memcpy(tempStorage, samples, sizeof(double) * (*length));
        free(*samples);
        *samples = tempStorage;
        *length = pow2length;
    }

    int pow2lengthOver2 = pow2length / 2;;
    
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    DSPDoubleSplitComplex tempSplitComplex;
    tempSplitComplex.realp = (double *) calloc(pow2lengthOver2, sizeof(double));
    tempSplitComplex.imagp = (double *) calloc(pow2lengthOver2, sizeof(double));
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)*samples, 2, &tempSplitComplex, 1, pow2lengthOver2);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(_fftSetup, &tempSplitComplex, 1, log2N, kFFTDirection_Forward);
    
    // Unpack the result into a real vector
    vDSP_ztocD(&tempSplitComplex, 1, (DSPDoubleComplex *)*samples, 2, pow2lengthOver2);
    
    free(tempSplitComplex.realp);
    free(tempSplitComplex.imagp);
}

- (void)performInverseFFTOnSamples:(double **)samples withLength:(int *)length {
    unsigned long log2N = ceil(log2(*length));
    int pow2length = (int) pow(2, log2N);
    if(pow2length != *length) {
        double* tempStorage = (double*) calloc(pow2length, sizeof(double));
        memcpy(tempStorage, samples, sizeof(double) * (*length));
        free(*samples);
        *samples = tempStorage;
        *length = pow2length;
    }
    int pow2lengthOver2 = pow2length / 2;
    
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    DSPDoubleSplitComplex tempSplitComplex;
    tempSplitComplex.realp = (double *) calloc(pow2lengthOver2, sizeof(double));
    tempSplitComplex.imagp = (double *) calloc(pow2lengthOver2, sizeof(double));
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)*samples, 2, &tempSplitComplex, 1, pow2lengthOver2);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(_fftSetup, &tempSplitComplex, 1, log2N, kFFTDirection_Inverse);
    
    // Unpack the result into a real vector
    vDSP_ztocD(&tempSplitComplex, 1, (DSPDoubleComplex *)*samples, 2, pow2lengthOver2);
    double scale = (1.0 / (((double) *length) * 2));
    vDSP_vsmulD(*samples, 1, &scale, *samples, 1, *length); // account for 1/N scaling factor
    
    free(tempSplitComplex.realp);
    free(tempSplitComplex.imagp);
}

- (void)performHilbertPhaseShiftsOnSamples:(double *)samples withLength:(int)length {
    for(int i = 0; i < length/2; i++) {
        samples[i] *= -1.0000;
    }
    double temp;
    for(int i = 2; i<length - 1; i+=2) {
        temp = samples[i];
        samples[i] = -1.0 * samples[i+1];
        samples[i+1] = samples[i];
    }
}

// convolves data vector with window the length of data vector must be at least lResult + lFilter - 1.
- (void)convolveSamples:(double *)samples withLength:(int)samplesLength withFilter:(double *)filter withFilterLength:(int)filterLength
           putResultsIn:(double *)results forResultLength:(int)resultLength {
    // start at end of filter to perform convolution
    vDSP_convD(samples, 1, filter + filterLength -1, -1, results,  1, resultLength, filterLength);
}

// stores result of hiblert envelop in array1
- (void)calculateHilbertEnvelopeForSamples:(double *)samples andTransformArray:(double *)transformArray withLength:(int)length {
    for(int i = 0; i < length; i ++) {
        double arr0sqrd = pow(samples[i], 2.0);
        double arr1sqrd = pow(transformArray[i], 2.0);
        transformArray[i] = sqrt(arr0sqrd + arr1sqrd);
    }
}

- (void) scaleData: (double*) data Length:(int) length Max:(double) max Min:(double) min {
    double dMax, dMin;
    dMax = dMin = data[0];
    for(int i = 1; i < length; i ++) {
        if(data[i] > dMax) {
            dMax = data[i];
        } else if(data[i] < dMin) {
            dMin = data[i];
        }
    }
    for(int i = 0; i < length; i ++) {
        data[i] = (((max - min) * (data[i] - dMin)) / (dMax - dMin)) + min;
    }
}

void ForwardLinearPrediction(std::vector<double> &coeffs, const std::vector<double> &x) {
    // Get size from input vectors
    size_t N = x.size() - 1;
    size_t m = coeffs.size();
    
    // Initialize r with autocorrelation coefficients
    std::vector<double> R( m + 1, 0.0 );
    for ( size_t i = 0; i <= m; i++ ) {
        for ( size_t j = 0; j <= N - i; j++ ) {
            R[ i ] += x[ j ] * x[ j + i ];
        }
    }
    
    // initialize Ak
    std::vector<double> Ak( m + 1, 0.0 );
    Ak[ 0 ] = 1.0;

    // initialize Ek
    double Ek = R[ 0 ];

    // Levinson-Durbin Recursion
    for ( size_t k = 0; k < m; k++)
    {
        // Compute lambda
        double lambda = 0.0;
        for ( size_t j = 0; j <= k; j++) { 
            lambda -= Ak[ j ] * R[ k + 1 - j ];
        } 
        lambda /= Ek;
        
        // Update Ak
        for ( size_t n = 0; n <= ( k + 1 ) / 2; n++) {
            double temp = Ak[ k + 1 - n ] + lambda * Ak[ n ];
            Ak[ n ] = Ak[ n ] + lambda * Ak[ k + 1 - n ];
            Ak[ k + 1 - n ] = temp;
        }
        
        // Update Ek
        Ek *= 1.0 - lambda * lambda;
    }
    
    // Assign coefficients
    coeffs.assign( ++Ak.begin(), Ak.end() ); 
}

@end