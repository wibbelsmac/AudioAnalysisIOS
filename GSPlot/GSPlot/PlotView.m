//
//  PlotView.m
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "PlotView.h"

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
    NSMutableArray *xValues = [[NSMutableArray alloc] init];
    NSMutableArray *yValues = [[NSMutableArray alloc] init];
    for (double i = 0; i < 10.0; i += 0.01) {
        [xValues addObject:[NSNumber numberWithDouble:i]];
        [yValues addObject:[NSNumber numberWithDouble:sin(i)]];
    }
    
    GSGraph *graph = [[GSGraph alloc] initWithXValues:xValues YValues:yValues withContext:UIGraphicsGetCurrentContext()];
    [graph drawAxisLines];
    [graph plotXAndYValues];
}

@end
