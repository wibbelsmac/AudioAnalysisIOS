//
//  AAAButterLow.h
//  AudioAnalysisApp
//
//  Created by Mackenzie Wibbels on 7/8/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import <Foundation/Foundation.h>
#define PI 3.14159265
@interface AAAButterLow : NSObject {

}
@property (nonatomic) double fr;
@property (nonatomic) double a0;
@property (nonatomic) double omegaPrime;
@property (nonatomic) double a1;
@property (nonatomic) double a2;
@property (nonatomic) double b1;
@property (nonatomic) double b2;
@property (nonatomic) double c;

-(id)initWithSampleFreq:(double)sampleFreq_ CutoffFreq:(double)cuttofff;
-(void) filterArray: (double*)dataArray DataLength:(int)dataLength ResultArray:(double*)resultArray ResultLength:(int)resultLength;
@end
