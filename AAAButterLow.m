//
//  AAAButterLow.m
//  AudioAnalysisApp
//
//  Created by Mackenzie Wibbels on 7/8/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import "AAAButterLow.h"
#import <Accelerate/Accelerate.h>


@implementation AAAButterLow
@synthesize fr;
@synthesize omegaPrime;
@synthesize a0;
@synthesize a1;
@synthesize a2;
@synthesize b1;
@synthesize b2;
@synthesize c;

-(id)initWithSampleFreq:(double)sampleFreq_ CutoffFreq:(double)cuttofff {
    // setting up the sample coefficients from following example:
    // http://www.kwon3d.com/theory/filtering/fil.html
    
    self.fr = sampleFreq_ / cuttofff;
    self.omegaPrime = tan (PI / fr);
    self.c = 1.0f + 2.0f * cos(PI/ 4.0f) * omegaPrime + powf(omegaPrime, 2.0f);
    self.a0 = self.a2 = powf(omegaPrime, 2.0f) / self.c;
    self.a1 = 2.0f * self.a0;
    self.b1 = 2.0f * (powf(omegaPrime, 2.0f) - 1)/ self.c;
    self.b2 =  (1.0f - 2.0f * cos(PI/ 4.0f) * omegaPrime + powf(omegaPrime, 2.0f)) / self.c;
    
    printf("a0: %f\n", a0);
    printf("a1: %f\n", a1);
    printf("a2: %f\n", a2);
    printf("b1: %f\n", b1);
    printf("b2: %f\n", b2);
    printf("total: %f\n", a0 + a1 + a2 - b1 - b2);
    return self;
}

-(void) filterArray: (double*)dataArray DataLength:(int)dataLength ResultArray:(double*)resultArray ResultLength:(int)resultLength {
    double *coefficients = (double *) malloc(sizeof(double) * 5);
    
    coefficients[0] = a0;
    coefficients[1] = a1;
    coefficients[2] = a2;
    coefficients[3] = b1;
    coefficients[4] = b2;

    // initializing result vector for recusive filtering
    resultArray[0] = dataArray [0];
    resultArray[1] = dataArray [1];
    
    vDSP_deq22D(dataArray, 1, coefficients, resultArray, 1, (resultLength - 2));
}
@end
