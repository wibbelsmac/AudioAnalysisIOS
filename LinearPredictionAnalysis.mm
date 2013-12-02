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
#include "LinearPredictionAnalysis.h"
#include "AAAButterLow.h"

#import <Accelerate/Accelerate.h>
#import "GSGraph.h"


using namespace std;
@interface LinearPredictionAnalysis(Private)

#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define localMin(a, dAta) (dAta[a -1] > dAta[a] && dAta[a + 1] > dAta[a])
#define localMax(a, dAta) (dAta[a -1] < dAta[a] && dAta[a + 1] < dAta[a])
#define mag(a) (sqrt(pow(a, 2.0)))
#define numFormants 4
typedef struct {
    int wordStartIndex;
    int vowelStartIndex;
    int maxIndex;
    int wordEndIndex;
    int vowelEndIndex;
    double totalSlope;
} audioSegment;

typedef struct {
    double formant[numFormants];
} formantData;

typedef struct {
    double pitch;
    int startIndex;
    int endIndex;
    double pitchCorrelation;
} PitchStruct;


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
-(double) SumOfSlopeAroundPoint:(int) index FodData:(double*) fodData FodLength:(int) length SecondsFrom:(double) secs;
-(void) CalculateFormantsFromLPC:(double*)data Length:(int) length;
@end

@implementation LinearPredictionAnalysis

static const int sampleFreq = 44100;
static const float sampleFrameTime = .02f;
static int numSamplesPerFrame = sampleFrameTime * sampleFreq;
static const float subFrameFract = .60f;
static const int numSamplesPerSubFrame = numSamplesPerFrame * subFrameFract;
static const int diffSubandFrame = numSamplesPerFrame - numSamplesPerSubFrame;
static const int LPCorder = 44;
- (id)initWithPackets:(SInt16 *)packets PacketCount:(int)packetCount
{
    self = [super init];
    if (self) {
        _mLog2N = ceil(log2(packetCount));
        _mNumberOfPackets = (int) pow(2, _mLog2N);
        _mFFTLength = _mNumberOfPackets / 2;
        NSLog(@"n = %d", _mNumberOfPackets);
        
        _fftSetup = vDSP_create_fftsetupD(_mLog2N, kFFTRadix2);
        
        // Initialize the input buffer with packets from the audio file
        _mAudioBuffer = (double *) calloc(_mNumberOfPackets, sizeof(double));
        for (int k = 0; k < _mNumberOfPackets; k++) {
            _mAudioBuffer[k] = (double) (packets[k]);
        }
        
        double * filteredData = (double*) calloc(_mNumberOfPackets, sizeof(double));
        // initialize a butterworth low pass filter with a cutoff frequency of 8KHz
        AAAButterLow *blow = [[AAAButterLow alloc] initWithSampleFreq:sampleFreq CutoffFreq:8000];
        [blow filterArray:_mAudioBuffer DataLength:_mNumberOfPackets ResultArray:filteredData ResultLength:_mNumberOfPackets];
        _mAudioBuffer = filteredData;
        
        [self CalculateVOP:_mAudioBuffer PacketCount:_mNumberOfPackets];
        
        
    }
    free(packets);
    return self;
}

- (void)CalculateVOP: (double*) data PacketCount:(int) pCount
{
    double hWindow_time = .02;
    double overlap_scale = .5;
    int windowLength =  hWindow_time * sampleFreq;
    
    int totalwindowSize = pow(2, ceil(log2(windowLength)));

    
    int numoverlapWindows = ((pCount - windowLength) / windowLength) /overlap_scale;
    int totalSampleProcessed = numoverlapWindows * windowLength / 2;
    // teporary storage for meta data
    double * window = (double*) calloc(totalwindowSize, sizeof(double));
    double * LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * HilTrans_LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * lpc_coef = (double*) malloc(sizeof(double) * LPCorder);
    // storage for result data
    double * resultData = (double*) calloc(pCount, sizeof(double));
    double * conv_resultData = (double*) calloc(pCount, sizeof(double));
    double * VOP_Evidience = (double*) calloc(pCount, sizeof(double));
    double* frameWindow = (double*) calloc(windowLength, sizeof(double)); // construct hamming window
    
    vDSP_hamm_windowD(frameWindow, windowLength, 0);
    [self scaleData:data Length:pCount Max:1.0 Min:-1.0];
    for(int i = 0; i < numoverlapWindows; i ++) {
        
        // copy data over into temp window
        int window_start = i * windowLength / 2;
        //vDSP_vmulD(frameWindow, 1, (data + window_start), 1, window, 1, windowLength);
        memcpy(
            window,
            (data + window_start),
            windowLength * sizeof(double)
        );
        
        //  1.  calculate the LP Coefficients
        //        [self ForwardLinearPredictionwithData:window DataLength:windowLength Coeffs:lpc_coef + 1 Order:LPCorder];
        [
            self
            FindLinearPredictionCoeffwithFrameArray:window
            ResultArray:lpc_coef
            Order:LPCorder
            FrameLength: windowLength
        ];
        
        //  2.  calculate the LP Residual
        [
            self
            inverseFilter:window
            Result:LP_Residual
            FilterCoefficients:lpc_coef // construct LPC residual
            ArrayLength:windowLength
            FilterLength:LPCorder
        ];
        
        //  3. copy data over and perform hilbert transform
        memcpy(HilTrans_LP_Residual, LP_Residual, totalwindowSize * sizeof(double)); // copy lpc residual for Hilbert Trans of LPC Resid
        [
            self
            applyHanningWindowTo:window
            StartingAt:0 EndingAt:windowLength
        ];
        
        [
            self
            performFFTForward:&HilTrans_LP_Residual
            Length:&totalwindowSize
        ]; // forward fft
        [
            self
            performHilbertPhaseShifts:HilTrans_LP_Residual
            Length:totalwindowSize
        ]; // hilbert transform phase shifts
        
        [
            self
            performFFTInverse:&HilTrans_LP_Residual
            Length:&totalwindowSize
        ]; // complete hilbert trans
        [self hilberEnvelope:LP_Residual Transform:HilTrans_LP_Residual Length:windowLength]; // envelope both values toghether
        
        for(int z = 0; z < windowLength; z++) {  // copy over data to result
            if(z > windowLength / 2) {
                resultData[window_start + z] +=  HilTrans_LP_Residual[z];  // scale to .5 for overlap
            } else {
                resultData[window_start + z] = HilTrans_LP_Residual[z];
            }
        }
    }
    // 50ms hamming window to smooth data
    int convWindowLength = .05 * sampleFreq;
    double* hWindow = (double*) calloc(convWindowLength, sizeof(double)); // construct hamming window
    vDSP_hamm_windowD(hWindow, convWindowLength, 0);
    [self convolveData:resultData DataLength:totalSampleProcessed Filter:hWindow FilterLength:convWindowLength
                Result:conv_resultData  ResultLength:totalSampleProcessed];
    
    [self scaleData:resultData Length:totalSampleProcessed Max:1.0 Min:-1.0];
    [self scaleData:conv_resultData Length:totalSampleProcessed Max:1.0 Min:0];

//    // segments audio into word and vowel segments
    double* fodData = (double*) malloc(sizeof(double) * pCount);
    double* segData = (double*) calloc(pCount, sizeof(double));
    double* segData1 = (double*) calloc(pCount, sizeof(double));
    memcpy(segData1, conv_resultData, totalSampleProcessed * sizeof(double));
    vector<audioSegment> segments = [self locateWordPeaksandBoundsWithData:conv_resultData Length:totalSampleProcessed];
    [self enchanceHilbertEnvelopeWithData:conv_resultData Length:totalSampleProcessed WordSegments:segments];
    [self generateVOPEvidiencewithData:conv_resultData Length:totalSampleProcessed WindowLength:(sampleFreq * .1)];
    for(int i = 0; i < segments.size(); i++) {
        [self print_temp_seg:&(segments[i]) Index:i];
       
    }
    [self locateVowelBoundariesWithVOPEvidience:conv_resultData Length:totalSampleProcessed WordSegments:segments];
        printf("\n");
    for(int i = 0; i < segments.size(); i++) {
        [self print_temp_seg:&(segments[i]) Index:i];
        for(int j = segments[i].vowelStartIndex; j < segments[i].vowelEndIndex; j += (int)(sampleFreq * .02)) {
            printf("%f\n",[self performCrossCorrelation:(data + j) NumSamples:(int)(sampleFreq * .02)]);
        }
        segData[segments[i].wordStartIndex] = .5;
        segData[segments[i].vowelStartIndex] = .75;
        segData[segments[i].maxIndex] = 1.0;
        segData[segments[i].vowelEndIndex] = .75;
        segData[segments[i].wordEndIndex] = .5;
    }
    
//    for(int j = 0; j < (totalSampleProcessed); j++) {
//        printf("%f\t%f\t%f\t%f\t%f\n", segData1[j], conv_resultData[j], data[j], fodData[j], segData[j]);
//    }
    
  //  [self plotAllValuesStartingPacket:0 EndWithPacket:6000 forPackets:resultData];
    
}
-(void) hilberTransformwithData:(double*) inputData ResultData:(double*) resultData Length:(int*) length {
    memcpy(resultData, inputData, (*length) * sizeof(double)); // copy lpc residual for Hilbert Trans of LPC Resid
    [self applyHanningWindowTo:resultData StartingAt:0 EndingAt:*length];
    
    [self performFFTForward:&resultData Length:length]; // forward fft
    [self performHilbertPhaseShifts:resultData Length:*length]; // hilbert transform phase shifts
    [self performFFTInverse:&resultData Length:length]; // complete hilbert trans
}


- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev
{
    int range = ev - sv;
    double *window = new double[range];
    vDSP_hann_windowD(window, range, vDSP_HANN_NORM);
    
    for (int i = 0, index = sv; index < ev; i++, index++) {
        values[index] *= window[i];
    }
}

- (void)dealloc
{
    vDSP_destroy_fftsetupD(_fftSetup);
    free(_mAudioBuffer);
}

- (void)plotRange:(int)range forPackets:(double *)packets
{
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    int resolution = (int) ceil(log2(range));
    resolution = CLAMP(resolution, 10, 20);
    for (int i = 0, index = 0; index < range; i++, index += resolution) {
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

-(vector<PitchStruct>) calculatePitchWithData:(double*) data Length:(int) len MinPitch:(double) minP MaxPitch:(double) maxP
                       WindowTime:(double) winTime WindowShift:(double) winShift AudioSegment:(audioSegment) audSeg {
    vector<PitchStruct> pitch;
    int winLength = winTime * sampleFreq;
    int winShiftLen  = winShift * sampleFreq;
    int start = max(0, (audSeg.vowelStartIndex - (winLength /2)));
    int end = min(len, (audSeg.vowelEndIndex + (winLength /2)));
    for(int i = start; i < end - (winLength / 2); i += winShift) {
        
    }
    return pitch;
}

- (double)performCrossCorrelation:(double *)packets NumSamples:(int) numSamples
{
    const int sampleCount = numSamples;
    // Populate *window with the values for a hamming window function
    double *hannWindow = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *tempSamples = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *subFrameSamples = (double *)malloc(sizeof(double) * numSamplesPerSubFrame);
    double *resultWindow = (double *)malloc(sizeof(double) * diffSubandFrame);
    
    vDSP_hamm_windowD(hannWindow, numSamplesPerFrame, 0);
    
    
    for (int i = 0; i <= sampleCount - numSamplesPerFrame; i += numSamplesPerFrame) {
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
//        for(int i = 0; i < )
        //return [self findPitchPeakFreq:resultWindow NumFrames:diffSubandFrame];
    }
    return 0.0;
}

- (int) findMinorMaxinRange:(double*) data Length:(int) len Start:(int) sIndex End:(int) eIndex MaxorMin:(bool) maxOrMin {
    int maxIndex = sIndex;
    double max = data[maxIndex];
    for(int i = sIndex; i < eIndex; i ++) {
        if(data[i] > max) {
            maxIndex = i;
            max = data[maxIndex];
        }
    }
    return maxIndex;
}

-(double) findPitchPeakFreq:(double*)result NumFrames:(int) numFrames {
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

-(void) FindLinearPredictionCoeffwithFrameArray:(double*) frame
                                    ResultArray:(double*) result
                                          Order:(int) p
                                    FrameLength:(int) frameLength{
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
    dgetrf_(&dim, &dim, corMatrix, &LDA, ipiv, &info); // invert matrix
    dgetrs_(&trans, &dim, &nrhs, corMatrix, &LDA, ipiv, result, &LDB, &info); // mutiply inverted matrix to solve lpc coefficients
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
    // the data vector must be atlest resultLength + filterLength
    [self convolveData:data DataLength:al Filter:filter FilterLength:inverseFilterLength Result:result ResultLength:(al)];
    
}

-(void) performFFTForward:(double**)data Length:(int*) length {
    unsigned long log2N = ceil(log2(*length));
    int pow2length = (int) pow(2, log2N);
    if(pow2length != *length) {
        double* tempStorage = (double*) calloc(pow2length, sizeof(double));
        memcpy(tempStorage, data, sizeof(double) * (*length));
        free(*data);
        *data = tempStorage;
        *length = pow2length;
    }
    //    [self applyHanningWindowTo:*data StartingAt:0 EndingAt:length -1];
    int pow2lengthOver2 = pow2length / 2;;
    
    // We need complex buffers for real and imaginary parts. These buffers efficiently store values
    // by getting rid of redundant values in real and imaginary parts of complex numbers.
    DSPDoubleSplitComplex tempSplitComplex;
    tempSplitComplex.realp = (double *) calloc(pow2lengthOver2, sizeof(double));
    tempSplitComplex.imagp = (double *) calloc(pow2lengthOver2, sizeof(double));
    
    // Transforms the real array input = {A[0],..., A[n]} into an even-odd
    // array splitComplex = {A[0], A[2],..., A[n-1], A[1], A[3],..., A[n]}.
    // splitComplex stores the imaginary and real parts in separate arrays
    vDSP_ctozD((DSPDoubleComplex *)*data, 2, &tempSplitComplex, 1, pow2lengthOver2);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(_fftSetup, &tempSplitComplex, 1, log2N, kFFTDirection_Forward);
    
    // Unpack the result into a real vector
    vDSP_ztocD(&tempSplitComplex, 1, (DSPDoubleComplex *)*data, 2, pow2lengthOver2);
    
    free(tempSplitComplex.realp);
    free(tempSplitComplex.imagp);
}

-(void) performFFTInverse:(double**)data Length:(int*) length {
    unsigned long log2N = ceil(log2(*length));
    int pow2length = (int) pow(2, log2N);
    if(pow2length != *length) {
        double* tempStorage = (double*) calloc(pow2length, sizeof(double));
        memcpy(tempStorage, data, sizeof(double) * (*length));
        free(*data);
        *data = tempStorage;
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
    vDSP_ctozD((DSPDoubleComplex *)*data, 2, &tempSplitComplex, 1, pow2lengthOver2);
    
    // Computes an in-place single-precision real discrete Fourier transform,
    // from the time domain to the frequency domain (forward).
    vDSP_fft_zripD(_fftSetup, &tempSplitComplex, 1, log2N, kFFTDirection_Inverse);
    
    // Unpack the result into a real vector
    vDSP_ztocD(&tempSplitComplex, 1, (DSPDoubleComplex *)*data, 2, pow2lengthOver2);
    double scale = (1.0 / (((double) *length) * 2));
    vDSP_vsmulD(*data, 1, &scale, *data, 1, *length); // account for 1/N scaling factor
    
    free(tempSplitComplex.realp);
    free(tempSplitComplex.imagp);
}
// estimates hilber transform by multipling entries by -j for ( index < (N/2) - 1)
// and by j for ( index > (N/2) - 1)
-(void) performHilbertPhaseShifts:(double*)packedData Length:(int)length {
    for(int i = 0; i < length/2; i++) {
        packedData[i]*= -1.0000;
    }
    double temp;
    for(int i = 2; i<length - 1; i+=2) {
        temp = packedData[i];
        packedData[i] = -1.0 * packedData[i+1];
        packedData[i+1] = packedData[i];
    }
}

// convolves data vector with window the length of data vector must be at least lResult + lFilter - 1.
-(void) convolveData: (double*) data DataLength:(int) dLength Filter:(double*) filter FilterLength:(int)
fLength Result:(double*) result ResultLength:(int) resultLength {
    // Convolution of edge of vector
    static int pad_length = fLength * 1.5 + 1;
    static double* conv_pad = (double*) malloc(pad_length * sizeof(double));
    if(pad_length < ((fLength * 1.5) -1)) {
        free(conv_pad);
        pad_length = fLength * 1.5;
        conv_pad = (double*) malloc(pad_length * sizeof(double));
    }
    int filterOver2 = fLength / 2;
    // lots of logic for the edge cases
    memset(conv_pad, 0, pad_length * sizeof(double)); // set padding to 0
    memcpy((conv_pad + filterOver2), data, fLength * sizeof(double)); // copy over first filter over2 values at end of pad
                                                                      //to avoid seg fault and convolve all values
    vDSP_convD(conv_pad, 1, filter + fLength -1, -1, result,  1, filterOver2,fLength); // perform conv on padded values
    
    memset(conv_pad, 0, pad_length * sizeof(double)); // set padding to 0
    memcpy(conv_pad, data + (resultLength - fLength), fLength * sizeof(double)); // copy over first filter over 2 values at begin of pad
                                                                                 // to avoid seg fault and convolve all values
    vDSP_convD(conv_pad, 1, filter + fLength -1, -1, result + (resultLength - filterOver2),  1, filterOver2,fLength);
    
    
    // Bulk of Convolution
    vDSP_convD(data, 1, filter + fLength -1, -1, result + filterOver2,  1, resultLength - fLength,fLength);
}

// stores result of hiblert envelop in array1
- (void) hilberEnvelope: (double*) array0 Transform: (double*) array1 Length:(int) length {
    for(int i = 0; i < length; i ++) {
        double arr0sqrd = pow(array0[i], 2.0);
        double arr1sqrd = pow(array1[i], 2.0);
        array1[i] = sqrt(arr0sqrd + arr1sqrd);
    }
}

// scales the data in an array between min and max
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
    if(dMax != dMin) {
        for(int i = 0; i < length; i ++) {
            data[i] = (((max - min) * (data[i] - dMin)) / (dMax - dMin)) + min;
        }
    } else {
        for(int i = 0; i < length; i ++) {
            data[i] = (max + min) / 2;
        }
    }
}
// computes first order central difference on data array
- (void) firstOrderDifference:(double*)data Length:(int)length Result:(double*)result Mean:(double*) mean SampleWindowTime:(double) time {
    *mean = 0.0;
    int sampleWindow = sampleFreq * time;
    int samWindowOver2 = sampleWindow >> 1;
    int start, end;
    for(int i = 0; i < length; i ++) {
        start = (i - samWindowOver2) > 0 ? (i - samWindowOver2) : 0;
        end =   (i + samWindowOver2) < length ? (i + samWindowOver2) : length;
        for(int j = start; j < end; j++) {
            if(i != j)
                result[i] += (data[i] - data[j]) / (i - j);
        }
        result[i] /= (end - start);
        *mean += abs(result[i]);
        
    }
    (*mean) /= length;
}

- (vector<audioSegment>) locateWordPeaksandBoundsWithData:(double*)data Length:(int)length {
    
    double pThresh = .5;
    double mean;
    // locate peaks and work boundaries
    double* fodData = (double*) malloc(sizeof(double) * length);
    [self firstOrderDifference:data Length:length Result:fodData Mean:&mean SampleWindowTime:.001];
    vector<audioSegment> tempSegs = [self locatePeakswithFODData:fodData Data:data Length:length Threshold:pThresh MeanSlope:mean];
    [self wordSegmentBoundaries:fodData DataLength:length AroundPeaks:tempSegs];
    free(fodData);
    return tempSegs;
}

// enchances hilber envelope by scaling word segments from 1.0
- (void) enchanceHilbertEnvelopeWithData:(double*) data Length:(int) length WordSegments:(vector<audioSegment>&) wordSegs{
    for(int i = 0; i < wordSegs.size(); i++) {
        int start = wordSegs[i].wordStartIndex;
        int end = wordSegs[i].wordEndIndex;
        [self scaleData:(data + start) Length:(end - start) Max:1.0 Min:0];
    }
}

// generate VOP Evidience by convolving hilbert envelope of data with first order Gaussian Differentiator
- (void) generateVOPEvidiencewithData:(double*) data Length:(int) length WindowLength:(int) winLength{
    double* tempData = (double*) malloc(length * sizeof(double));
    double* window = (double*) malloc(winLength * sizeof(double));
    [self firstOderGaussianWindow:window Length:winLength];
    memcpy(tempData, data, length * sizeof(double));
    [self convolveData:tempData DataLength:length Filter:window FilterLength:winLength
                Result:data  ResultLength:length];
    [self scaleData:data Length:length Max:1.0 Min:-1.0];
    free(tempData);
    free(window);
}

// locate peaks in VOP evidience waveform data 
- (void) locateVowelBoundariesWithVOPEvidience:(double*) vopEvd Length:(int) length WordSegments:(vector<audioSegment>&) wordSegs{
    int segment = 0;
    int index = 0;
    bool foundVoweStart, foundVowelEnd;
    foundVoweStart = foundVowelEnd = false;
    double secondsRange = .0005;
    double threshold = .5;
    while ((segment < wordSegs.size()) && (index < length - 1)) {
        
        if((vopEvd[index] > threshold) && localMax(index, vopEvd)
           && (index >= wordSegs[segment].wordStartIndex)) {
            if([self localMinOrMaxwithData:vopEvd Length:length AtIndex:index WithinSeconds:secondsRange MinOrMax:false]) {
                wordSegs[segment].vowelStartIndex = index;
                foundVoweStart = true;
            }
        } else if (foundVoweStart && vopEvd[index] < 0
                   && (index <= wordSegs[segment].wordEndIndex)) {
            if([self localMinOrMaxwithData:vopEvd Length:length AtIndex:index WithinSeconds:secondsRange  MinOrMax:true]) {
                wordSegs[segment].vowelEndIndex = index;
                foundVowelEnd = true;
            }
        }
        
        if(foundVoweStart && foundVowelEnd) {
            segment++;
            foundVoweStart = foundVowelEnd = false;
        } else if((segment < (wordSegs.size() - 1)) && index >  wordSegs[segment + 1].wordStartIndex) {
            printf("Eliminating word at:%d after reaching index:%d\n", segment, index);
            if(!foundVoweStart) {
                printf("NO Vowel Start above thresh\n");
                wordSegs[segment].vowelStartIndex = 0;
                wordSegs[segment].vowelEndIndex = 0;
            } else if(!foundVowelEnd) {
                printf("NO Vowel end below thresh\n");
                  wordSegs[segment].vowelEndIndex = 0;
            }
            [self print_temp_seg:&(wordSegs[segment]) Index:segment];
            //wordSegs.erase(wordSegs.begin() + segment);
            segment++;
            foundVoweStart = foundVowelEnd = false;
        }
        index++;
    }
}

// false = MAX true = MIN 
-(bool) localMinOrMaxwithData:(double*) data Length:(int) length AtIndex:(int) index WithinSeconds:(double) secs MinOrMax:(bool) min{
    int range = secs * sampleFreq;
    int start = (index - range) > 0 ? (index - range) : 0;
    int end =   (index + range) < length ? (index + range) : length;
    if(end == length || start == 0) {
        printf("trimmed summation window\n");
    }
    if(!min) {
        for(int i = start; i < end; i++) {
            if(data[i] > data[index])
                return false;
        }
    } else {
        for(int i = start; i < end; i++) {
            if(data[i] < data[index])
                return false;
        }
    }
    return true;
}

// false = MAX true = MIN
-(bool) localMinOrMaxwithFODData:(double*) fodData Length:(int) length AtIndex:(int) index WithinSeconds:(double) secs MinOrMax:(bool) min{
    int range = secs * sampleFreq;
    int start = (index - range) > 0 ? (index - range) : 0;
    int end =   (index + range) < length ? (index + range) : length;
    if(end == length || start == 0) {
        printf("trimmed summation window\n");
    }
    if(!min) {
        for(int i = start; i < index; i++) {
            if(fodData[i] > 0)
                return false;
        }
        for(int i = index + 1; i < end; i++) {
            if(fodData[i] < 0)
                return false;
        }
    } else {
        for(int i = start; i < index; i++) {
            if(fodData[i] < 0)
                return false;
        }
        for(int i = index + 1; i < end; i++) {
            if(fodData[i] > 0)
                return false;
        }
    }
    return true;
}

-(vector<audioSegment>) locatePeakswithFODData:(double*) fodData  Data:(double*) data Length:(int) length Threshold:(double) thresh MeanSlope:(double) mSlope{
    vector<audioSegment> tempVect;
    audioSegment tempSeg;
    bool foundPrevSeg = false;
    int minBeetwPeaks = .05 * sampleFreq;
    double sumofSlopeTimeOver2 = .005;
    int sumofSlopeWindow = 2 * sumofSlopeTimeOver2 * sampleFreq;
    for(int i = 0; i < length; i++) {
        if(fodData[i] <= 0 && fodData[i -1] > 0) {
            double totalSlope = [self SumOfSlopeAroundPoint:i FodData:fodData FodLength:length SecondsFrom:sumofSlopeTimeOver2];
            if(abs(totalSlope) > (thresh * sumofSlopeWindow * mSlope) && totalSlope > 0) {
                // potential local max
                if(!foundPrevSeg) { // if this is the first segment located set tempseg index and total slope
                    tempSeg.maxIndex = i;
                    tempSeg.totalSlope = abs(totalSlope);
                    foundPrevSeg = true;
                } else if((i - tempSeg.maxIndex) > minBeetwPeaks) { // segments are farther than min distance away
                    tempVect.push_back(tempSeg); // push back first temp set
                    tempSeg.maxIndex = i; // set temp seg with index and slope
                    tempSeg.totalSlope = abs(totalSlope);
                } else if(abs(totalSlope) > tempSeg.totalSlope) { // if conflict between two segs and slope of new index is higher
                    tempSeg.maxIndex = i; // update temp seg
                    tempSeg.totalSlope = abs(totalSlope);
                }
            }
        }
    }
    if(foundPrevSeg) // push on pending segment to vector
        tempVect.push_back(tempSeg);
    return tempVect;
}
-(void) wordSegmentBoundaries:(double*) fodData DataLength:(int) length AroundPeaks:(vector<audioSegment>&) tempSegs {
    for(int i =0; i < tempSegs.size(); i++) {
        [self locateNearestNegtoZeroCrossingAround:tempSegs[i] FODData:fodData FODLength:length];
    }
}

-(void) locateNearestNegtoZeroCrossingAround:(audioSegment&) tempSeg FODData:(double*) data FODLength:(int) length {
    double windowThresh = .002;
    int startIndex, endIndex;
    startIndex = endIndex = tempSeg.maxIndex;
    startIndex = tempSeg.maxIndex;
    endIndex = tempSeg.maxIndex;
    while (startIndex > 0) {
        startIndex--;
        if((data[startIndex] < 0) && (data[startIndex + 1] >0)
           && [self localMinOrMaxwithFODData:data Length:length AtIndex:startIndex WithinSeconds:windowThresh MinOrMax:false]) {
            printf("Found min start: slope:%f to:%f\n", data[startIndex], data[startIndex + 1]);
            break;
        }
    }
    while (endIndex < length) {
        endIndex++;
        if((data[endIndex] < 0) && (data[endIndex + 1] >0)
           && [self localMinOrMaxwithFODData:data Length:length AtIndex:endIndex WithinSeconds:windowThresh MinOrMax:false]) {
            printf("Found min end: slope:%f to:%f\n", data[endIndex], data[endIndex + 1]);
            break;
        }
    }
    tempSeg.wordStartIndex = startIndex;
    tempSeg.wordEndIndex = endIndex;

}

-(double) SumOfSlopeAroundPoint:(int) index FodData:(double*) fodData FodLength:(int) length SecondsFrom:(double) secs {
    int halfSlopeSumLength =  secs * sampleFreq;
    int start = (index - halfSlopeSumLength) > 1 ? (index - halfSlopeSumLength) : 1;
    int end = (index + halfSlopeSumLength) < length ? (index + halfSlopeSumLength) : length - 1;
    if(end == length || start == 0) {
        printf("trimmed summation window");
    }
    double sum = 0.0;
    for(int i = start; i <= end; i++) {
        sum += (i < index) ? fodData[i] : (-1.0 * fodData[i]);
    }
    return sum;
}
- (void) firstOderGaussianWindow:(double*)window Length:(int) length {
    double signma = 551;
//    double omega  = .00206;
    double last, current;
    last = window[0] = (1.0 / (sqrt(2 * PI) * signma)) * exp( -1.0 * pow((0 - (length/2)),2.0) / (2 * pow(signma, 2.0)));
    for(int i =1; i < length ; i++) {
        current = (1.0 / (sqrt(2 * PI) * signma)) * exp( -1.0 * pow((i - (length/2)),2.0) / (2 * pow(signma, 2.0)));
        window[i] = current - last;
        last = current;
    }
}

-(void) ForwardLinearPredictionwithData:(double*)x DataLength:(int) length Coeffs:(double*) coeffs Order:(int) order{
    // GET SIZE FROM INPUT VECTORS
    size_t N = length - 1;
    size_t m = order;
    // INITIALIZE R WITH AUTOCORRELATION COEFFICIENTS
    std::vector<double> R( m + 1, 0.0 );
    for ( size_t i = 0; i <= m; i++ )
    {
        for ( size_t j = 0; j <= N - i; j++ )
        {
            R[ i ] += x[ j ] * x[ j + i ];
        }
    }
    // INITIALIZE Ak
    std::vector<double> Ak( m + 1, 0.0 );
    Ak[ 0 ] = 1.0;
    // INITIALIZE Ek
    double Ek = R[ 0 ];
    // LEVINSON-DURBIN RECURSION
    for ( size_t k = 0; k < m; k++ )
    {
        // COMPUTE LAMBDA
        double lambda = 0.0;
        for ( size_t j = 0; j <= k; j++ )
        {
            lambda -= Ak[ j ] * R[ k + 1 - j ];
        }
        lambda /= Ek;
        // UPDATE Ak
        for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ )
        {
            double temp = Ak[ k + 1 - n ] + lambda * Ak[ n ];
            Ak[ n ] = Ak[ n ] + lambda * Ak[ k + 1 - n ];
            Ak[ k + 1 - n ] = temp;
        }
        // UPDATE Ek
        Ek *= 1.0 - lambda * lambda;
    }
    // ASSIGN COEFFICIENTS
    memcpy(coeffs, &(Ak[0]),  Ak.size() * sizeof(double));
}
-(void) CalculateFormantsWithLPC:(double*) lpcData LPCLength:(int) lpcLength LPCWindowTime:(double) lpcWinTime
                                 WindowTime:(double) winTime WindowOverlap:(double) overlap Result:(formantData*) result ResultLength:(int) rLength {
    int lpcWinLength = sampleFreq * lpcWinTime;
    int formantCount = lpcLength / LPCorder;
    int start = 0;
    for (int i = 0; i < formantCount; i++) {
        start = i * LPCorder;
        [self CalculateFormantsFromLPC:(lpcData + start) Length:LPCorder];
    }
}

-(formantData) CalculateFormantsFromLPC:(double*)data Length:(int) length {
    static int pow2length = 0;
    static DSPDoubleSplitComplex tempSplitComplex;
    static long log2N = 0;
    if(pow2length < length) {
        if(pow2length != 0) {
            free(tempSplitComplex.realp);
            free(tempSplitComplex.imagp);
        }
        log2N = ceil(log2(sampleFreq));
        pow2length = (int) pow(2, log2N);
        tempSplitComplex.realp = (double *) calloc(pow2length>>1, sizeof(double));
        tempSplitComplex.imagp = (double *) calloc(pow2length>>1, sizeof(double));
    }
    [self applyHanningWindowTo:data StartingAt:0 EndingAt:length -1];
    int pow2lengthOver2 = pow2length / 2;
    double* data1 = (double*) malloc((length + 1) * sizeof(double));
    memcpy(data1 + 1, data, length * sizeof(double));
    data1[0] = 1.0;

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
        tempSplitComplex.realp[i] = 1 / sqrt(arr0sqrd + arr1sqrd);
    }
    
int index, formant;
index = formant = 0;
    formantData fdTemp;
while(index < pow2length && formant < 4) {
    if([self localMinOrMaxwithData:tempSplitComplex.realp Length:pow2lengthOver2 AtIndex:index WithinSeconds:.00025 MinOrMax:false]) {
        fdTemp.formant[formant] = ((index * sampleFreq) / pow2length);
    }
}
return fdTemp;
}

-(void) print_temp_seg:(audioSegment*) audSeg Index:(int) index{
    double wStart, wEnd, maxIndex, vStart, vEnd;
    wStart = audSeg->wordStartIndex;
    wEnd = audSeg->wordEndIndex;
    maxIndex = audSeg->maxIndex;
    vStart = audSeg->vowelStartIndex;
    vEnd   = audSeg->vowelEndIndex;
    printf("Seg:%d start:%E vStart:%E peak:%E vEnd:%E end:%E wLeng:%E vLeng%E\n", index, wStart, vStart, maxIndex, vEnd, wEnd, ((wEnd - wStart)/sampleFreq), (vEnd - vStart)/sampleFreq);
}
@end