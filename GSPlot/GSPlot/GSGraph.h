//
//  GSGraph.h
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "GSAxisLine.h"

typedef enum AxisType {
    kXAxis,
    kYAxis
} AxisType;

@interface GSGraph : NSObject

@property (nonatomic) CGContextRef currentContext;
@property (nonatomic, strong) NSArray *xValues;
@property (nonatomic, strong) NSArray *yValues;
@property (nonatomic, strong) GSAxisLine *xAxisLine;
@property (nonatomic, strong) GSAxisLine *yAxisLine;
@property (nonatomic) CGPoint zeroPoint;

- (id)initWithXValues:(NSArray *)xVals YValues:(NSArray *)yVals withContext:(CGContextRef)context;
- (void)drawAxisLines;
- (void)plotXAndYValues;

@end
