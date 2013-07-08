//
//  FFTManager.h
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#include <stdio.h>
#include <Accelerate/Accelerate.h>

@interface FFTManager : NSObject

- (void)performFFT:(SInt16 *)packets ByteCount:(int) byteCount;

@end