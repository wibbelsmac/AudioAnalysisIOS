//
//  FFTManager.h
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#include <stdio.h>
#include <Accelerate/Accelerate.h>

@interface FFTManager : NSObject

@property (nonatomic) FFTSetupD fftSetup;
@property (nonatomic) int mLog2N;
@property (nonatomic) int mNumberOfPackets;
@property (nonatomic) int mFFTLength;
@property (nonatomic) double *mAudioBuffer;
@property (nonatomic) DSPDoubleSplitComplex mSplitComplex;

- (id)initWithPackets:(SInt16 *)packets PacketCount:(int)packetCount;
- (void)performFFTOnRangeStartingAt:(int)sp EndingAt:(int)ep;
- (void)performVoiceAnalysisOn:(double *)frames;

@end