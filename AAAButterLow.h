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
@property (nonatomic) float fr;
@property (nonatomic) float a0;
@property (nonatomic) float omegaPrime;
@property (nonatomic) float a1;
@property (nonatomic) float a2;
@property (nonatomic) float b1;
@property (nonatomic) float b2;
@property (nonatomic) float c;

-(id)initWithSampleFreq:(float)sampleFreq_ CutoffFreq:(float)cuttofff;
-(void) filterArray: (float*)dataArray DataLength:(int)dataLength ResultArray:(float*)resultArray ResultLength:(int)resultLength;
@end
