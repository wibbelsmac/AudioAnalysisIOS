//
//  GSGraph.m
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "GSGraph.h"

const double kGraphMargin = 10.0;

@interface GSGraph(Private)

- (double)getScreenRangeForNegativeValuesWith:(double)firstNum and:(double)secondNum forRange:(double)range;
- (double)scaleValueToViewCoordinates:(double)value forAxis:(AxisType)axis;

@end

@implementation GSGraph

- (id)initWithXValues:(NSArray *)xVals YValues:(NSArray *)yVals withContext:(CGContextRef)context
{
    self = [super init];
    if (self) {
        if (self.xValues.count != self.yValues.count) {
            NSLog(@"Both arrays must have the same number of elements!");
        }
        
        [self setXValues:xVals];
        [self setYValues:yVals];
        [self setCurrentContext:context];
        
        // Find the smallest and largest X and Y values
        double smallestX = [[self.xValues objectAtIndex:0] doubleValue];
        double largestX = [[self.xValues objectAtIndex:0] doubleValue];
        double smallestY = [[self.yValues objectAtIndex:0] doubleValue];
        double largestY = [[self.yValues objectAtIndex:0] doubleValue];
        for (int i = 1; i < self.xValues.count; i++) {
            double currentXValue = [[self.xValues objectAtIndex:i] doubleValue];
            double currentYValue = [[self.yValues objectAtIndex:i] doubleValue];
            if (currentXValue < smallestX) {
                smallestX = currentXValue;
            } else if (currentXValue > largestX) {
                largestX = currentXValue;
            }
            
            if (currentYValue < smallestY) {
                smallestY = currentYValue;
            } else if (currentYValue > largestY) {
                largestY = currentYValue;
            }
        }
        
        CGSize viewSize = [UIApplication sharedApplication].keyWindow.rootViewController.view.bounds.size;
        
        // Create the x-axis based on the input values.
        if (smallestY >= 0) {
            self.xAxisLine = [[GSAxisLine alloc] initWithStartPoint:CGPointMake(kGraphMargin, viewSize.height - kGraphMargin) endPoint:CGPointMake(viewSize.width - kGraphMargin, viewSize.height - kGraphMargin) inContext:self.currentContext];
        } else {
            double screenRangeForNegativeValues = [self getScreenRangeForNegativeValuesWith:smallestY and:largestY forRange:viewSize.height];
            self.xAxisLine = [[GSAxisLine alloc] initWithStartPoint:CGPointMake(kGraphMargin, viewSize.height - kGraphMargin - screenRangeForNegativeValues) endPoint:CGPointMake(viewSize.width - kGraphMargin, viewSize.height - kGraphMargin - screenRangeForNegativeValues) inContext:self.currentContext];
        }
        [self.xAxisLine setMinValue:smallestX];
        [self.xAxisLine setMaxValue:largestX];
        [self.xAxisLine setWidth:2.0f];
        
        // Create the y-axis based on the input values.
        if (smallestX >= 0) {
            self.yAxisLine = [[GSAxisLine alloc] initWithStartPoint:CGPointMake(kGraphMargin, viewSize.height - kGraphMargin) endPoint:CGPointMake(kGraphMargin, kGraphMargin) inContext:self.currentContext];
        } else {
            double screenRangeForNegativeValues = [self getScreenRangeForNegativeValuesWith:smallestX and:largestX forRange:viewSize.width];
            self.yAxisLine = [[GSAxisLine alloc] initWithStartPoint:CGPointMake(screenRangeForNegativeValues + kGraphMargin, viewSize.height - kGraphMargin) endPoint:CGPointMake(screenRangeForNegativeValues + kGraphMargin, kGraphMargin) inContext:self.currentContext];
        }
        [self.yAxisLine setMinValue:smallestY];
        [self.yAxisLine setMaxValue:largestY];
        [self.yAxisLine setWidth:2.0f];
        
        self.zeroPoint = CGPointMake([self.yAxisLine startPoint].x, [self.xAxisLine startPoint].y);
    }
    return self;
}

- (double)getScreenRangeForNegativeValuesWith:(double)firstNum and:(double)secondNum forRange:(double)range
{
    double actualRange = fabsf(firstNum) + fabsf(secondNum);
    double actualToScreenRatio = actualRange / range;
    return fabsf(firstNum) / actualToScreenRatio;
}

- (void)drawAxisLines
{
    [self.xAxisLine drawAxisLine];
    [self.yAxisLine drawAxisLine];
}

- (void)plotXAndYValues
{
    CGContextSetStrokeColorWithColor(self.currentContext, [UIColor redColor].CGColor);
    
    CGContextMoveToPoint(self.currentContext, [self scaleValueToViewCoordinates:[self.xValues[0] floatValue] forAxis:kXAxis], [self scaleValueToViewCoordinates:[self.yValues[0] floatValue] forAxis:kYAxis]);
    for (int i = 1; i < self.xValues.count; i++) {
        CGContextAddLineToPoint(self.currentContext, [self scaleValueToViewCoordinates:[self.xValues[i] floatValue] forAxis:kXAxis], [self scaleValueToViewCoordinates:[self.yValues[i] floatValue] forAxis:kYAxis]);
    }
    CGContextStrokePath(self.currentContext);
}

- (double)scaleValueToViewCoordinates:(double)value forAxis:(AxisType)axis
{
    double numerator, denominator;
    switch (axis) {
        case kXAxis:
            numerator = ([self.xAxisLine endPoint].x - [self.xAxisLine startPoint].x) * (value - [self.xAxisLine minValue]);
            denominator = [self.xAxisLine maxValue] - [self.xAxisLine minValue];
            return (numerator / denominator) + [self.xAxisLine startPoint].x;
            
        default:
            numerator = ([self.yAxisLine endPoint].y - [self.yAxisLine startPoint].y) * (value - [self.yAxisLine minValue]);
            denominator = [self.yAxisLine maxValue] - [self.yAxisLine minValue];
            return (numerator / denominator) + [self.yAxisLine startPoint].y;
    }
}

@end
