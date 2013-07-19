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

- (void)plotRange:(int)range forPackets:(double *)packets;
- (void)plotAllValuesStartingPacket:(int)sp EndWithPacket:(int)ep forPackets:(double *)packets;
- (void)applyHanningWindowTo:(double *)values StartingAt:(int)sv EndingAt:(int)ev;
-(void) FindLinearPredictionCoeffwithFrameArray:(double*) frame ResultArray:(double*) result Order:(int) p FrameLength:(int) frameLength;
-(double*) crreatZeroBufferedwindowNonZeroLength:(int) nonZerloLengh ZeroLeng:(int)zeroPadding Data:(double*) data;
- (void)performCrossCorrelation:(double *)packets NumSamples:(int) numSamples;
-(void) performFFTForward:(double**)data Length:(int*) length;
-(void) performFFTInverse:(double**)data Length:(int*) length;
-(void) performHilbertPhaseShifts:(double*)packedData Length:(int)length;
-(void)CalculateVOP: (double*) data PacketCount:(int) pCount;
-(void) hilberEnvelope: (double*) array0 Transform: (double*) array1 Length:(int) length;
-(void) convolveData: (double*) data DataLength:(int) dLength Filter:(double*) filter FilterLength:(int)
fLength Result:(double*) result ResultLength:(int) resultLength;
void ForwardLinearPrediction( std::vector<double> &coeffs, const std::vector<double> &x );
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
        [self CalculateVOP:_mAudioBuffer PacketCount:_mNumberOfPackets];
        
        //    double * filteredData = (double*) calloc(pCount, sizeof(double));
        //    AAAButterLow *blow = [[AAAButterLow alloc] initWithSampleFreq:sampleFreq CutoffFreq:2500];
        //    [blow filterArray:data DataLength:pCount ResultArray:filteredData ResultLength:pCount];
        //    data = filteredData;
    }
    free(packets);
    return self;
}
- (void)CalculateVOP: (double*) data PacketCount:(int) pCount
{
    double hWindow_time = .02;
    double overlap_scale = .5;
    int hWindowLength =  hWindow_time * sampleFreq;  // hamming window length
    // calculating linear prediction coefficients
    int totalwindowSize = pow(2, ceil(log2(hWindowLength)));
    int LPCorder = 10;

    int numoverlapWindows = ((pCount - hWindowLength) / hWindowLength) /overlap_scale;
    int totalSampleProcessed = numoverlapWindows * hWindowLength / 2;
    
    double * window = (double*) calloc(totalwindowSize, sizeof(double));
    double * LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * HilTrans_LP_Residual = (double *) calloc(totalwindowSize, sizeof(double));
    double * lpc_coef = (double*) malloc(sizeof(double) * LPCorder);
    double * resultData = (double*) calloc(pCount, sizeof(double));
    double * conv_resultData = (double*) calloc(pCount, sizeof(double));

    for(int i = 0; i < numoverlapWindows; i ++) {
        // copy data over into temp window
        int window_start = i * hWindowLength / 2;
        memcpy(window, (data + window_start),hWindowLength * sizeof(double));
        std::vector<double> coeffs (LPCorder, 0.0);
        std::vector<double> dataVector (window, window + hWindowLength);

        //  1.  calculate the LP Coefficients
        ForwardLinearPrediction(coeffs, dataVector);
//        [self FindLinearPredictionCoeffwithFrameArray:window ResultArray:lpc_coef Order:LPCorder FrameLength: hWindowLength];
        memcpy(lpc_coef, &(coeffs[0]),LPCorder * sizeof(double));
        
        //  2.  calculate the LP Residual
        [self inverseFilter:window Result:LP_Residual FilterCoefficients:lpc_coef // construct LPC residual
                ArrayLength:hWindowLength FilterLength:LPCorder];
        //  3. copy data over and perform hilbert transform
        memcpy(HilTrans_LP_Residual, LP_Residual, totalwindowSize * sizeof(double)); // copy lpc residual for Hilbert Trans of LPC Resid
        [self applyHanningWindowTo:window StartingAt:0 EndingAt:hWindowLength];
        [self performFFTForward:&HilTrans_LP_Residual Length:&totalwindowSize]; // forward fft
        [self performHilbertPhaseShifts:HilTrans_LP_Residual Length:totalwindowSize]; // hilbert transform phase shifts
        [self performFFTInverse:&HilTrans_LP_Residual Length:&totalwindowSize]; // complete hilbert trans
        [self hilberEnvelope:LP_Residual Transform:HilTrans_LP_Residual Length:hWindowLength]; // envelope both values toghether
        
        for(int z = 0; z < hWindowLength; z++) {  // copy over data to result
            if(z > hWindowLength / 2) {
               resultData[window_start + z] +=  HilTrans_LP_Residual[z];  // scale to .5 for overlap
            } else {
                resultData[window_start + z] = HilTrans_LP_Residual[z];
            }
        }
    }
      // 50ms hanning window to smooth data after
    int convWindowLength = .05 * sampleFreq;
    double* hWindow = (double*) calloc(convWindowLength, sizeof(double)); // construct hamming window
    vDSP_hamm_windowD(hWindow, convWindowLength, 0);
//    [self convolveData:resultData DataLength:totalSampleProcessed Filter:hWindow FilterLength:convWindowLength
//          Result:conv_resultData  ResultLength:(totalSampleProcessed - hWindowLength)];
    
    [self scaleData:data Length:pCount Max:1.0 Min:-1.0];
    [self scaleData:resultData Length:pCount Max:1.0 Min:-1.0];
    [self scaleData:conv_resultData Length:pCount Max:1.0 Min:-1.0];
    for(int j = 0; j < (totalSampleProcessed); j++) {
        printf("%d\t%f\t%f\n", j, resultData[j], data[j]);
    }
    
    [self plotAllValuesStartingPacket:0 EndWithPacket:6000 forPackets:resultData];

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

- (void)performCrossCorrelation:(double *)packets NumSamples:(int) numSamples
{    
    const int sampleCount = numSamples;
    // Populate *window with the values for a hamming window function
    double *hannWindow = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *tempSamples = (double *)malloc(sizeof(double) * numSamplesPerFrame);
    double *subFrameSamples = (double *)malloc(sizeof(double) * numSamplesPerSubFrame);
    double *resultWindow = (double *)malloc(sizeof(double) * diffSubandFrame);
    
    vDSP_hamm_windowD(hannWindow, numSamplesPerFrame, 0);
    
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
        vDSP_convD(tempSamples,1, subFrameSamples,1,resultWindow, 1, diffSubandFrame, numSamplesPerSubFrame);
        [peakValues addObject:[NSNumber numberWithInt:[self findPitchPeakFreq:resultWindow NumFrames:diffSubandFrame]]];
    }
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

-(double*) crreatZeroBufferedwindowNonZeroLength:(int) nonZerloLengh ZeroLeng:(int)zeroPadding Data:(double *) data{
    double * buff = (double*) calloc(nonZerloLengh + zeroPadding, sizeof(double));
    memcpy(buff, data, sizeof(double) * nonZerloLengh);
    return buff;
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
    // start at end of filter to perform convolution
    vDSP_convD(data, 1, filter + fLength -1, -1, result,  1, resultLength,fLength);
}

// stores result of hiblert envelop in array1
- (void) hilberEnvelope: (double*) array0 Transform: (double*) array1 Length:(int) length {
    for(int i = 0; i < length; i ++) {
        double arr0sqrd = pow(array0[i], 2.0);
        double arr1sqrd = pow(array1[i], 2.0);
        array1[i] = sqrt(arr0sqrd + arr1sqrd);
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

void ForwardLinearPrediction( std::vector<double> &coeffs, const std::vector<double> &x )
{
    // GET SIZE FROM INPUT VECTORS
    size_t N = x.size() - 1;
    size_t m = coeffs.size();
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
    coeffs.assign( ++Ak.begin(), Ak.end() ); 
}
@end