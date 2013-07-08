//
//  AAAAudioFile.h
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface AAAAudioFile : NSObject {

}
@property (nonatomic) int byteCount;
- (SInt16 *)open:(NSString *)fileName ofType:(NSString *)fileType;

@end
