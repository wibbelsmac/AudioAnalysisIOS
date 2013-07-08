//
//  GSGraphLine.m
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "GSAxisLine.h"

@implementation GSAxisLine

@synthesize minorIndicatorSize, majorIndicatorSize, minValue, maxValue;

- (id)initWithStartPoint:(CGPoint)start endPoint:(CGPoint)end inContext:(CGContextRef)context
{
    self = [super init];
    if (self) {
        [self setStartPoint:start];
        [self setEndPoint:end];
        [self setMinValue:0.0f];
        [self setMaxValue:0.0f];
        [self setMinorIndicatorSize:0.0f];
        [self setMajorIndicatorSize:0.0f];
        [self setCurrentContext:context];
    }
    return self;
}

- (void)drawAxisLine
{
    CGContextSetLineWidth(self.currentContext, self.width);
    CGContextSetStrokeColorWithColor(self.currentContext, [UIColor blackColor].CGColor);
    CGContextMoveToPoint(self.currentContext, self.startPoint.x, self.startPoint.y);
    CGContextAddLineToPoint(self.currentContext, self.endPoint.x, self.endPoint.y);
    CGContextStrokePath(self.currentContext);
}

@end
