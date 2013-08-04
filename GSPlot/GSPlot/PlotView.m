//
//  PlotView.m
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "PlotView.h"
#import "AAAAudioFile.h"
#import "FFTManager.h"
#import "LinearPredictionAnalysis.h"

@implementation PlotView

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code
    }
    return self;
}

- (void)drawRect:(CGRect)rect
{
    
    AAAAudioFile *audioFile = [[AAAAudioFile alloc] init];
//    SInt16* packets = [audioFile open:[[NSBundle mainBundle] pathForResource:@"250Hz_44100Hz_16bit_05sec" ofType:@"wav"] ofType:@"wav"];
    SInt16* packets = [audioFile open:[[NSBundle mainBundle] pathForResource:@"test1" ofType:@"wav"] ofType:@"wav"];
    
    LinearPredictionAnalysis *lpcManager = [[LinearPredictionAnalysis alloc] initWithPackets:packets PacketCount:audioFile.packetCount];
//    [fftManager performFFTOnRangeStartingAt:(int) 44100.0 * 0.8 EndingAt:(int) 44100.0 * 1.1];
//    [fftManager performFFTOnRangeStartingAt:0 EndingAt:fftManager.mNumberOfPackets];
}

@end
