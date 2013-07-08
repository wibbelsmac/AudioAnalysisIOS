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

@implementation PlotView

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code
    }
    return self;
}

// Only override drawRect: if you perform custom drawing.
// An empty implementation adversely affects performance during animation.
- (void)drawRect:(CGRect)rect
{
//    NSMutableArray *xValues = [[NSMutableArray alloc] init];
//    NSMutableArray *yValues = [[NSMutableArray alloc] init];
//    for (double i = 0; i < 10.0; i += 0.01) {
//        [xValues addObject:[NSNumber numberWithDouble:i]];
//        [yValues addObject:[NSNumber numberWithDouble:sin(i)]];
//    }
    
    AAAAudioFile *audioFile = [[AAAAudioFile alloc] init];
    SInt16* packets = [audioFile open:[[NSBundle mainBundle] pathForResource:@"250Hz_44100Hz_16bit_05sec" ofType:@"wav"] ofType:@"wav"];
//    SInt16* packets = [audioFile open:[[NSBundle mainBundle] pathForResource:@"test2" ofType:@"wav"] ofType:@"wav"];
    
    FFTManager *fftManager = [[FFTManager alloc] init];
    [fftManager performFFT:packets ByteCount:audioFile.byteCount];
    
//    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
//    [graph drawAxisLines];
//    [graph plotXAndYValues];
}

@end
