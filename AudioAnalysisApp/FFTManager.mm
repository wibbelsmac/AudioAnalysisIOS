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


using namespace std;
@interface FFTManager(Private)

<<<<<<< HEAD
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define localMin(a, dAta) (dAta[a -1] > dAta[a] && dAta[a + 1] > dAta[a])
#define localMax(a, dAta) (dAta[a -1] < dAta[a] && dAta[a + 1] < dAta[a])
#define mag(a) (sqrt(pow(a, 2.0)))

typedef struct {
    int minStartIndex;
    int maxIndex;
    int minEndIndex;
    double totalSlope;
} vowelSeg;

typedef struct {
    int start;
    int end;
} vowelBoundaries;

- (void)plotRange:(int)range forPackets:(double *)packets;
- (void)plotAllValuesStartingPacket:(int)sp EndWithPacket:(int)ep forPackets:(double *)packets;
- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev;
-(void) FindLinearPredictionCoeffwithFrameArray:(double*) frame ResultArray:(double*) result Order:(int) p FrameLength:(int) frameLength;

-(double)performCrossCorrelation:(double *)packets NumSamples:(int) numSamples;
-(void) performFFTForward:(double**)data Length:(int*) length;
-(void) performFFTInverse:(double**)data Length:(int*) length;
-(void) performHilbertPhaseShifts:(double*)packedData Length:(int)length;
-(void)CalculateVOP: (double*) data PacketCount:(int) pCount;
-(void) hilberEnvelope: (double*) array0 Transform: (double*) array1 Length:(int) length;
-(void) convolveData: (double*) data DataLength:(int) dLength Filter:(double*) filter FilterLength:(int)
fLength Result:(double*) result ResultLength:(int) resultLength;
- (vector<vowelBoundaries>) pickPeaksandValleys:(double*)data Length:(int) length PeakThreshold:(double) pThresh ValleyThreshold:(double) vThresh;
-(double) SumOfSlopeAroundPoint:(int) index FodData:(double*) fodData FodLength:(int) length SecondsFrom:(double) secs;
-(void) CalculateFormantsFromLPC:(double*)data Length:(int) length;
=======
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

>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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
        
<<<<<<< HEAD
=======
        // An opaque type that contains setup information for a given double-precision FFT transform.
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
        _fftSetup = vDSP_create_fftsetupD(_mLog2N, kFFTRadix2);
        
        // Initialize the input buffer with packets from the audio file
        _mAudioBuffer = (double *) calloc(_mNumberOfPackets, sizeof(double));
        for (int k = 0; k < _mNumberOfPackets; k++) {
            _mAudioBuffer[k] = (double) (packets[k]);
        }
<<<<<<< HEAD
        
        double * filteredData = (double*) calloc(_mNumberOfPackets, sizeof(double));
        AAAButterLow *blow = [[AAAButterLow alloc] initWithSampleFreq:sampleFreq CutoffFreq:5500];
        [blow filterArray:_mAudioBuffer DataLength:_mNumberOfPackets ResultArray:filteredData ResultLength:_mNumberOfPackets];
        _mAudioBuffer = filteredData;
        
        [self CalculateVOP:_mAudioBuffer PacketCount:_mNumberOfPackets];
        
  
=======
        
        [self calculateVowelOnsetPointForSamples:_mAudioBuffer withLength:_mNumberOfPackets];
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
    }
    free(packets);
    return self;
}

- (void)calculateVowelOnsetPointForSamples:(double *)samples withLength:(int)length {
    double hanning_window_time = .02;
    double overlap_scale = .5;
<<<<<<< HEAD
    int windowLength =  hWindow_time * sampleFreq;

    // calculating linear prediction coefficients
    int totalwindowSize = pow(2, ceil(log2(windowLength)));
    int LPCorder = 44;

    int numoverlapWindows = ((pCount - windowLength) / windowLength) /overlap_scale;
    int totalSampleProcessed = numoverlapWindows * windowLength / 2;
    
    double * window = (double*) calloc(totalwindowSize, sizeof(double));
    double * LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * HilTrans_LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * lpc_coef = (double*) malloc(sizeof(double) * LPCorder);
    double * resultData = (double*) calloc(pCount, sizeof(double));
    double * conv_resultData = (double*) calloc(pCount, sizeof(double));
    double * pitchEst = (double*) calloc(pCount, sizeof(double));
    double * VOP_Evidience = (double*) calloc(pCount, sizeof(double));
    double* frameWindow = (double*) calloc(windowLength, sizeof(double)); // construct hamming window
    
    vDSP_hamm_windowD(frameWindow, windowLength, 0);
    bool calculateFormants = false;
    [self scaleData:data Length:pCount Max:1.0 Min:-1.0];
    
    for(int i = 0; i < numoverlapWindows; i ++) {
        
        // copy data over into temp window
        int window_start = i * windowLength / 2;
//        vDSP_vmulD(frameWindow, 1, (data + window_start), 1, window, 1, windowLength);
        memcpy(window, (data + window_start), windowLength * sizeof(double));

        //  1.  calculate the LP Coefficients
//        [self ForwardLinearPredictionwithData:window DataLength:windowLength Coeffs:lpc_coef + 1 Order:LPCorder];
        [self FindLinearPredictionCoeffwithFrameArray:window ResultArray:lpc_coef Order:LPCorder FrameLength: windowLength];
//        memcpy(lpc_coef, &(coeffs[0]),LPCorder * sizeof(double));
        if(window_start > 36000 && !calculateFormants) {
            [self CalculateFormantsFromLPC:lpc_coef Length:LPCorder];
            calculateFormants = true;
        }
        //  2.  calculate the LP Residual
        [self inverseFilter:window Result:LP_Residual FilterCoefficients:lpc_coef // construct LPC residual
                ArrayLength:windowLength FilterLength:LPCorder];
        //  3. copy data over and perform hilbert transform
        memcpy(HilTrans_LP_Residual, LP_Residual, totalwindowSize * sizeof(double)); // copy lpc residual for Hilbert Trans of LPC Resid
        [self applyHanningWindowTo:window StartingAt:0 EndingAt:windowLength];
        
        [self performFFTForward:&HilTrans_LP_Residual Length:&totalwindowSize]; // forward fft
        [self performHilbertPhaseShifts:HilTrans_LP_Residual Length:totalwindowSize]; // hilbert transform phase shifts
        [self performFFTInverse:&HilTrans_LP_Residual Length:&totalwindowSize]; // complete hilbert trans
        
        [self hilberEnvelope:LP_Residual Transform:HilTrans_LP_Residual Length:windowLength]; // envelope both values toghether
        for(int z = 0; z < windowLength; z++) {  // copy over data to result
            if(z > windowLength / 2) {
               resultData[window_start + z] +=  HilTrans_LP_Residual[z];  // scale to .5 for overlap
=======
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
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
            } else {
                result_data[window_start + z] = pHilbert_transform_LP_residual[z];
            }
        }
    }
<<<<<<< HEAD
=======
    
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
    // 50ms hanning window to smooth data after
    int convWindowLength = .05 * sampleFreq;
    double* hWindow = (double*) calloc(convWindowLength, sizeof(double)); // construct hamming window
    vDSP_hamm_windowD(hWindow, convWindowLength, 0);
<<<<<<< HEAD
    [self convolveData:resultData DataLength:totalSampleProcessed Filter:hWindow FilterLength:convWindowLength
          Result:conv_resultData  ResultLength:totalSampleProcessed];

    [self scaleData:resultData Length:totalSampleProcessed Max:1.0 Min:-1.0];
    [self scaleData:conv_resultData Length:totalSampleProcessed Max:1.0 Min:0];
    
    [self findVowelSegments:conv_resultData Length:totalSampleProcessed Result:VOP_Evidience];

//    for(int j = 0; j < (totalSampleProcessed)/10; j++) {
//        printf("%d\t%f\t%f\t%f\n", j, conv_resultData[j], data[j], (VOP_Evidience[j]));
//    }

    [self plotAllValuesStartingPacket:0 EndWithPacket:6000 forPackets:resultData];

}


- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev
{
    int range = ev - sv;
=======
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
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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

<<<<<<< HEAD
- (double)performCrossCorrelation:(double *)packets NumSamples:(int) numSamples
{    
    const int sampleCount = numSamples;
    // Populate *window with the values for a hamming window function
    double *hannWindow = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *tempSamples = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *subFrameSamples = (double *)malloc(sizeof(double) * numSamplesPerSubFrame);
    double *resultWindow = (double *)malloc(sizeof(double) * diffSubandFrame);
    
    vDSP_hamm_windowD(hannWindow, numSamplesPerFrame, 0);
    

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
        vDSP_convD(tempSamples,1, subFrameSamples,1,resultWindow, 1, diffSubandFrame, numSamplesPerSubFrame);
        return [self findPitchPeakFreq:resultWindow NumFrames:diffSubandFrame];
    }
    return 0.0;
}

-(double) findPitchPeakFreq:(double*)result NumFrames:(int) numFrames {
=======
-(int) findPitchPeakFreq:(double*)result NumFrames:(int) numFrames {
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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

<<<<<<< HEAD
-(void) FindLinearPredictionCoeffwithFrameArray:(double*) frame ResultArray:(double*) result Order:(int) p FrameLength:(int) frameLength{
    int col = 0;
    double* corMatrix = (double*) calloc(p * p, sizeof(double));
    memset(result, 0,p * sizeof(double));
    vDSP_convD(frame + 1, 1, frame, 1,  result, 1, p, frameLength);
    for(int i = 0; i < p; i++){
        result[i] = -1 * result[i];
    }
    for(int row =0; row < p; row++) {
        vDSP_convD(frame, 1, frame, 1, (corMatrix + row * p + col) , 1, p - col, frameLength);
            for(int temp = 1; temp < p - col; temp++) {
                corMatrix[row + col * p + temp * p] = corMatrix[row *p + col + temp];
            }
        col++;
    }
    // solve the system of equations for the linear prediction coefficients
    char trans = 'N';
    long int dim = p;
    long int nrhs = 1;
    long int LDA = dim;
    long int LDB = dim;
    long int info = 0;
    long int ipiv[p * p -1];
    dgetrf_(&dim, &dim, corMatrix, &LDA, ipiv, &info);
    dgetrs_(&trans, &dim, &nrhs, corMatrix, &LDA, ipiv, result, &LDB, &info);
}


=======
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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
<<<<<<< HEAD
-(void) convolveData: (double*) data DataLength:(int) dLength Filter:(double*) filter FilterLength:(int)
        fLength Result:(double*) result ResultLength:(int) resultLength {
    int filterOver2 = fLength / 2;
    // start at end of filter to perform convolution
    memcpy(result, data, filterOver2 * sizeof(double));
    memcpy(result + (resultLength - filterOver2), data + (resultLength - filterOver2), filterOver2 * sizeof(double));
    vDSP_convD(data, 1, filter + fLength -1, -1, result + filterOver2,  1, resultLength - filterOver2,fLength);
=======
- (void)convolveSamples:(double *)samples withLength:(int)samplesLength withFilter:(double *)filter withFilterLength:(int)filterLength
           putResultsIn:(double *)results forResultLength:(int)resultLength {
    // start at end of filter to perform convolution
    vDSP_convD(samples, 1, filter + filterLength -1, -1, results,  1, resultLength, filterLength);
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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
// computes first order central difference on data array
- (void) firstOrderDifference:(double*)data Length:(int)length Result:(double*)result Mean:(double*) mean SampleWindowTime:(double) time {
    *mean = 0.0;
    int sampleWindow = sampleFreq * time;
    int samWindowOver2 = sampleWindow >> 1;
    for(int i = samWindowOver2; i < length - samWindowOver2; i ++) {
        for(int j = i - samWindowOver2; j < i + samWindowOver2; j++) {
            if(i != j)
                result[i] += (data[i] - data[j]) / (i - j);
        }
        *mean += abs(result[i]);
    }
    (*mean) /= length;
}

- (void) findVowelSegments:(double*)data Length:(int)length Result:(double*)result {

    double pThresh = .125;
    double halfSlopeWindow = .005;
    double mean;

    [self firstOrderDifference:data Length:length Result:result Mean:&mean SampleWindowTime:.0001];
    int minBeetwPeaks = .05 * sampleFreq;
    vector<vowelSeg> tempSegs = [self locatePeakswithFODData:result Data:data Length:length Threshold:pThresh MeanSlope:mean];
   for(int i = 0; i < tempSegs.size(); i++) {
        result[tempSegs[i].maxIndex] = 1.0;
    }
//    [self scaleData:result Length:length Max:1.0 Min:0];
    vowelSeg tempSeg, prevSeg;
    bool minStart, max, minEnd, storedPrevSeg;
    
//    for(int i = 1; i < length - 1; i ++) {
        // if starting find minStart, max, and then minEnd
//        if((localMin(i, data)) && (!max && !minEnd)) { // set minStart and set minstart to true
//            tempSeg.minStartIndex = i;
//            minStart = true;
//        } else if((localMax(i, data) && data[i] &&(!max || data[i] > data[tempSeg.maxIndex]) && minStart && !minEnd) && (data[i] > pThresh)) {
//            tempSeg.totalSlope = [self SumOfSlopeAroundMax:i FodData:result FodLength:length SecondsFrom:halfSlopeWindow];
//            if((tempSeg.totalSlope / (sampleFreq * halfSlopeWindow * 2.0)) > (mean * .75)) {
//                tempSeg.maxIndex = i;
//                max = true;
//            }
//        } else if(!minEnd && localMin(i, data) && minStart && max  < vThresh) {
//            minEnd = i;
//            minEnd = true;
//        }
//        
//        // located potential vowel segment with peak and two local minima
//        if (minStart && max && minEnd) {
//            if(!storedPrevSeg || ((tempSeg.maxIndex - prevSeg.minEndIndex) >= minBeetwPeaks)) { // not within 50 ms
//                tempSegs.push_back(prevSeg);
//                prevSeg = tempSeg;
//                storedPrevSeg = true;
//            } else if(prevSeg.totalSlope > tempSeg.totalSlope) { // conflict but temp seg has higher slope
//                tempSegs.push_back(tempSeg);
//                prevSeg = tempSeg;
//            } else { // prev seg has higher slope
//                 tempSegs.push_back(prevSeg);
//            }
//            tempSeg.minStartIndex = tempSeg.minEndIndex;  // start min of temp seg is endmin of last tempSef
//            max = minEnd = false;
//        }
//    }
//    for(int i =0 ; i < tempSegs.size(); i++) {
//        vowelSeg temp = tempSegs[i];
//        [self scaleData:(data + temp.minStartIndex) Length:(temp.minEndIndex - temp.minStartIndex) Max:1.0 Min:0.0];
//    }
//    double* VOPEvidience = (double*) calloc(length, sizeof(double));
//    int hundredmsWindowLength = .1 * sampleFreq;
//    double* window = (double*) calloc(hundredmsWindowLength, sizeof(double));
//    [self firstOderGaussianWindow:window Length:hundredmsWindowLength];
//    [self convolveData:data DataLength:length Filter:window FilterLength:hundredmsWindowLength Result:VOPEvidience ResultLength:length];
//    [self scaleData:VOPEvidience Length:length Max:1.0 Min:-1.0];
//    vector<vowelBoundaries>  vowelSegs = [self pickPeaksandValleys:VOPEvidience Length:length PeakThreshold:.5 ValleyThreshold:-.25];

}

-(vector<vowelSeg>) locatePeakswithFODData:(double*) fodData  Data:(double*) data Length:(int) length Threshold:(double) thresh MeanSlope:(double) mSlope{
    vector<vowelSeg> tempVect;
    vowelSeg tempSeg;
    bool foundPrevSeg = false;
    int minBeetwPeaks = .05 * sampleFreq;
    int sumofSlopeWindow = 2 * .005 * sampleFreq;
    for(int i = 0; i < length; i++) {
        if(fodData[i] <= 0 && fodData[i -1] > 0) {
            double totalSlope = [self SumOfSlopeAroundPoint:i FodData:fodData FodLength:length SecondsFrom:.005];
            if(abs(totalSlope) > (thresh * sumofSlopeWindow * mSlope)) {
                // potential local max
                if(!foundPrevSeg) {
                    tempSeg.maxIndex = i;
                    tempSeg.totalSlope = abs(totalSlope);
                    foundPrevSeg = true;
                } else if((i - tempSeg.maxIndex) > minBeetwPeaks) {
                    tempVect.push_back(tempSeg);
                    tempSeg.maxIndex = i;
                    tempSeg.totalSlope = abs(totalSlope);
                } else if(mag(totalSlope) > tempSeg.totalSlope) { // if conflict between two segs
                    tempSeg.maxIndex = i;
                    tempSeg.totalSlope = abs(totalSlope);
                }
            }
        }
    }
    return tempVect;
}

-(double) SumOfSlopeAroundPoint:(int) index FodData:(double*) fodData FodLength:(int) length SecondsFrom:(double) secs {
    int halfSlopeSumLength =  secs * sampleFreq;
    int start = (index - halfSlopeSumLength) > 1 ? (index - halfSlopeSumLength) : 1;
    int end = (index + halfSlopeSumLength) < length ? (index + halfSlopeSumLength) : length - 1;
    
    double sum = 0.0;
    for(int i = start; i <= end; i++) {
//        if(index != i)
//            sum += (fodData[index] - fodData[i]) /  (index - i);
        sum += abs(fodData[index]);
    }
    return sum;
}
- (void) firstOderGaussianWindow:(double*)window Length:(int) length {
    double signma = 234;
    double omega  = .00206;
    for(int i =0; i < length ; i++) {
        window[i] = -1.0 / (sqrt(2 * PI) * signma) * exp(-1.0 * pow(i,2) / (2 * pow(signma, 2))) * cos(omega * i);
    }
}
- (vector<vowelBoundaries>) pickPeaksandValleys:(double*)data Length:(int) length PeakThreshold:(double) pThresh ValleyThreshold:(double) vThresh {
    vector<vowelBoundaries> tempVect;
    vowelBoundaries vb;
    bool fndlMax, fndLMin;
    for(int i = 0; i < length ; i++) {
        if(localMax(i, data) && (data[i] > pThresh)) {
            vb.start = i;
            fndlMax = true;
        } else if (fndlMax && localMin(i, data) && data[i] < vThresh) {
            vb.end = i;
            fndLMin = true;
        }
        
        if(fndlMax && fndLMin) {
            tempVect.push_back(vb);
            fndlMax = fndLMin = false;
        }
    }
    return tempVect;
}

<<<<<<< HEAD
//- (void)performFFTOnRangeStartingAt:(int)sp EndingAt:(int)ep
//{
//    _mLog2N = ceil(log2(ep - sp));
//    _mNumberOfPackets = (int) pow(2, _mLog2N);
//    _mFFTLength = _mNumberOfPackets / 2;
//
//    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
//    // by getting rid of redundant values in real and imaginary parts of complex numbers.
//    _mSplitComplex.realp = (double *) calloc(_mFFTLength, sizeof(double));
//    _mSplitComplex.imagp = (double *) calloc(_mFFTLength, sizeof(double));
//
//    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
//    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
//    // splitComplex stores the imaginary and real parts in separate arrays
//    vDSP_ctozD((DSPDoubleComplex *)_mAudioBuffer, 2, &_mSplitComplex, 1, _mFFTLength);
//
//    // Computes an in-place single-precision real discrete Fourier transform,
//    // from the time domain to the frequency domain (forward).
//    vDSP_fft_zripD(self.fftSetup, &_mSplitComplex, 1, _mLog2N, kFFTDirection_Forward);
//
//    [self plotRange:ep - sp forPackets:_mAudioBuffer];
//
//    double scale = 1.0 / (2 * _mNumberOfPackets);
//    vDSP_vsmulD(_mSplitComplex.realp, 1, &scale, _mSplitComplex.realp, 1, _mFFTLength);
//    vDSP_vsmulD(_mSplitComplex.imagp, 1, &scale, _mSplitComplex.imagp, 1, _mFFTLength);
//
//    _mSplitComplex.realp[0] = 0; // zeroing out the DC value
//    _mSplitComplex.imagp[0] = 0; // zeroing out the Nyquist value
//
//    double tmpData[_mFFTLength];
//    double mAdjust0BD = 1.5849e-13;
//    double one = 1.0;
//
//    // Convert from complex/rectangular (real & imaginary) coordinates
//    // to polar (magnitude & phase) coordinates
//    vDSP_zvmagsD(&_mSplitComplex, 1, tmpData, 1, _mFFTLength);
//    vDSP_vsaddD(tmpData, 1, &mAdjust0BD, tmpData, 1, _mFFTLength);
//    vDSP_vdbconD(tmpData, 1, &one, tmpData, 1, _mFFTLength, 0);
//
//    for (int i = 0; i < _mFFTLength; i++) {
//        if (tmpData[i] < 0.0) {
//            tmpData[i] = 0;
//        }
//    }
//
////    [self plotRange:_mFFTLength forPackets:tmpData];
//
//    DSPDoubleComplex *complex = new DSPDoubleComplex[_mFFTLength];
//
//    // For polar coordinates
//    double *mag = new double[_mFFTLength];
//    double *phase = new double[_mFFTLength];
//
//    // Get complex vector absolute values
//    vDSP_zvabsD(&_mSplitComplex, 1, mag, 1, _mFFTLength);
//    // Get complex vector phase
//    vDSP_zvphasD(&_mSplitComplex, 1, phase, 1, _mFFTLength);
//
//    _mSplitComplex.realp = mag;
//    _mSplitComplex.imagp = phase;
//
//    // Convert from polar coords back to rectangular coords
//    vDSP_ztocD(&_mSplitComplex, 1, complex, 2, _mFFTLength);
//    vDSP_rectD((double *)complex, 2, (double *)complex, 2, _mFFTLength); // Polar to rectangular conversion
//    vDSP_ctozD(complex, 2, &_mSplitComplex, 1, _mFFTLength);
//
//    // Computes an in-place single-precision real discrete Fourier transform,
//    // from the frequency domain to the time domain (inverse).
//    vDSP_fft_zripD(_fftSetup, &_mSplitComplex, 1, _mLog2N, kFFTDirection_Inverse);
//
//    // Unpack the result into a real vector
//    vDSP_ztocD(&_mSplitComplex, 1, (DSPDoubleComplex *)_mAudioBuffer, 2, _mLog2N);
//
//    // Compensate for scaling for both FFTs. See Apple's vDSP documentation
//    // to see how much scaling is required.
//    // NOTE: Apple scales values to increase the efficiency of their functions.
//    // As such, we have to scale them back.
//    vDSP_vsmulD(_mAudioBuffer, 1, &scale, _mAudioBuffer, 1, _mLog2N);
//
////    [self plotRange:ep - sp forPackets:_mAudioBuffer];
//    free(_mSplitComplex.realp);
//    free(_mSplitComplex.imagp);
//}

-(void) ForwardLinearPredictionwithData:(double*)x DataLength:(int) length Coeffs:(double*) coeffs Order:(int) order{
    // GET SIZE FROM INPUT VECTORS
    size_t N = length - 1;
    size_t m = order;
    // INITIALIZE R WITH AUTOCORRELATION COEFFICIENTS
=======
void ForwardLinearPrediction(std::vector<double> &coeffs, const std::vector<double> &x) {
    // Get size from input vectors
    size_t N = x.size() - 1;
    size_t m = coeffs.size();
    
    // Initialize r with autocorrelation coefficients
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
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
<<<<<<< HEAD
    } 
    // ASSIGN COEFFICIENTS 
    memcpy(coeffs, &(Ak[0]),  Ak.size() * sizeof(double));
}
-(void) CalculateFormantsFromLPC:(double*)data Length:(int) length{
    unsigned long log2N = ceil(log2(sampleFreq));
    int pow2length = (int) pow(2, log2N);
    //    [self applyHanningWindowTo:*data StartingAt:0 EndingAt:length -1];
    int pow2lengthOver2 = pow2length / 2;;
    double* data1 = (double*) malloc((length + 1) * sizeof(double));
    memcpy(data1 + 1, data, length * sizeof(double));
    data1[0] = 1.0;
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    DSPDoubleSplitComplex tempSplitComplex;
    tempSplitComplex.realp = (double *) calloc(pow2lengthOver2, sizeof(double));
    tempSplitComplex.imagp = (double *) calloc(pow2lengthOver2, sizeof(double));
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)data1, 2, &tempSplitComplex, 1, length + 1);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(_fftSetup, &tempSplitComplex, 1, log2N, kFFTDirection_Forward);
    
    for(int i = 0; i < pow2lengthOver2 /2; i ++) {
        double arr0sqrd = pow(tempSplitComplex.realp[i], 2.0);
        double arr1sqrd = pow(tempSplitComplex.imagp[i], 2.0);
        double tempDub = 1 / sqrt(arr0sqrd + arr1sqrd);
        printf("%d\t%f\t%f\t%f\n", i, tempDub, tempDub, tempDub);
    }
    
    free(tempSplitComplex.realp);
    free(tempSplitComplex.imagp);
=======
    }
    
    // Assign coefficients
    coeffs.assign( ++Ak.begin(), Ak.end() ); 
>>>>>>> 15b0188d3d36e2854d9889933a31201fba964f23
}

@end