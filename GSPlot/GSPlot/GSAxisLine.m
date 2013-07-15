//
//  GSGraphLine.m
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "GSAxisLine.h"

@interface GSAxisLine(Private)

- (void)drawAxisLine;
- (void)drawMaxAndMinIndicator;

@end

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

- (void)drawAxis
{
    CGContextSetLineWidth(self.currentContext, self.width);
    CGContextSetStrokeColorWithColor(self.currentContext, [UIColor blackColor].CGColor);
    
    [self drawAxisLine];
    [self drawMaxAndMinIndicator];
}

- (void)drawAxisLine
{
    CGContextMoveToPoint(self.currentContext, self.startPoint.x, self.startPoint.y);
    CGContextAddLineToPoint(self.currentContext, self.endPoint.x, self.endPoint.y);
    CGContextStrokePath(self.currentContext);
}

- (void)drawMaxAndMinIndicator
{
    if (self.startPoint.x == self.endPoint.x) {
        // Y-axis
        CGContextMoveToPoint(self.currentContext, self.startPoint.x - 5.0, self.startPoint.y);
        CGContextAddLineToPoint(self.currentContext, self.startPoint.x + 5.0, self.startPoint.y);
        CGContextStrokePath(self.currentContext);
        
        CGContextMoveToPoint(self.currentContext, self.endPoint.x - 5.0, self.endPoint.y);
        CGContextAddLineToPoint(self.currentContext, self.endPoint.x + 5.0, self.endPoint.y);
        CGContextStrokePath(self.currentContext);
        
        NSString *label = [NSString stringWithFormat:@"%.0f", self.minValue];
        [label drawAtPoint:CGPointMake(self.startPoint.x - 15.0, self.startPoint.y - 8.0) withFont:[UIFont fontWithName:@"Helvetica-Bold" size:12.0f]];
        
        label = [NSString stringWithFormat:@"%.0f", self.maxValue];
        [label drawAtPoint:CGPointMake(self.endPoint.x - 15.0, self.endPoint.y - 8.0) withFont:[UIFont fontWithName:@"Helvetica-Bold" size:12.0f]];
    } else {
        CGContextMoveToPoint(self.currentContext, self.startPoint.x, self.startPoint.y - 5.0);
        CGContextAddLineToPoint(self.currentContext, self.startPoint.x, self.startPoint.y + 5.0);
        CGContextStrokePath(self.currentContext);
        
        CGContextMoveToPoint(self.currentContext, self.endPoint.x, self.endPoint.y - 5.0);
        CGContextAddLineToPoint(self.currentContext, self.endPoint.x, self.endPoint.y + 5.0);
        CGContextStrokePath(self.currentContext);
        
        NSString *label = [NSString stringWithFormat:@"%.0f", self.minValue];
        [label drawAtPoint:CGPointMake(self.startPoint.x - 3.0, self.startPoint.y + 5.0) withFont:[UIFont fontWithName:@"Helvetica-Bold" size:12.0f]];
        
        label = [NSString stringWithFormat:@"%.0f", self.maxValue];
        [label drawAtPoint:CGPointMake(self.endPoint.x - 10.0, self.endPoint.y + 5.0) withFont:[UIFont fontWithName:@"Helvetica-Bold" size:12.0f]];
    }
}

- (void)validate
{
    if (self.startPoint.x == self.endPoint.x) {
        double range = self.startPoint.y - self.endPoint.y;
        [self setMajorIndicatorSize:range / 10.0];
        [self setMinorIndicatorSize:range / 20.0];
    } else {
        double range = self.endPoint.x - self.startPoint.x;
        [self setMajorIndicatorSize:range / 15.0];
        [self setMinorIndicatorSize:range / 30.0];
    }
}

@end
