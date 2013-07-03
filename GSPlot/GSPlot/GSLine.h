//
//  GSLine.h
//  GSPlot
//
//  Created by Arda Tugay on 7/2/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface GSLine : NSObject

@property (nonatomic) CGPoint startPoint;
@property (nonatomic) CGPoint endPoint;
@property (nonatomic) CGFloat width;
@property (nonatomic) CGContextRef currentContext;

@end
