//
//  GSGraphLine.h
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "GSLine.h"

@interface GSAxisLine : GSLine

@property (nonatomic) double minorIndicatorSize;
@property (nonatomic) double majorIndicatorSize;
@property (nonatomic) double maxValue;
@property (nonatomic) double minValue;

- (id)initWithStartPoint:(CGPoint)start endPoint:(CGPoint)end inContext:(CGContextRef)context;
- (void)drawAxis;
- (void)validate;

@end
