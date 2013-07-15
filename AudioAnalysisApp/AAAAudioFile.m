//
//  AAAAudioFile.m
//  AudioAnalysisApp
//
//  Created by Arda Tugay on 6/28/13.
//  Copyright (c) 2013 Arda Tugay. All rights reserved.
//

#import <AudioToolbox/AudioToolbox.h>

#import "AAAAudioFile.h"

@implementation AAAAudioFile

@synthesize packetCount;

- (SInt16 *)open:(NSString *)fileName ofType:(NSString *)fileType
{
    OSStatus result = -1;
    
    CFStringRef filePath = (__bridge CFStringRef)(fileName);
    
    CFURLRef audioFileURL = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, (CFStringRef)filePath, kCFURLPOSIXPathStyle, false);
 
    AudioFileID mAudioFile;
    // Open Audio File
    result = AudioFileOpenURL(audioFileURL, kAudioFileReadPermission, 0, &mAudioFile);
    if (result != noErr) {
        NSLog(@"Could not open file: %@", filePath);
    }
    AudioStreamBasicDescription mASBD;
    UInt32 size = sizeof(mASBD);
    result = AudioFileGetProperty(mAudioFile, kAudioFilePropertyDataFormat, &size, &mASBD);
    
    UInt64 fileSize;
    size = sizeof(fileSize);
    result = AudioFileGetProperty(mAudioFile, kAudioFilePropertyAudioDataByteCount, &size, &fileSize);
    NSLog(@"File size is %d bytes or %.2f megabytes", (unsigned int)fileSize, (float)fileSize * 9.53674e-7);
    
    UInt32 packets = (UInt32) fileSize / mASBD.mBytesPerPacket;
    const int log2n = ceil(log2(packets));
    const int sizeOfPacketArray = (int) pow(2, log2n);
    
    UInt32 packetsRead = packets;
    UInt32 numBytesRead = -1;
    SInt16 *audioData = nil;
    if (packets > 0) {
        // Allocate Buffer
        audioData = (SInt16 *) malloc(sizeOfPacketArray * sizeof(SInt16));

        // Read the Packets
        result = AudioFileReadPackets(mAudioFile, false, &numBytesRead, NULL, 0, &packetsRead, audioData);
        self.packetCount = packetsRead;
        NSLog(@"Read %d bytes or %d packets", (unsigned int)numBytesRead, (unsigned int)packetsRead);
        
        for (int j = packetsRead; j < sizeOfPacketArray; j++) {
            audioData[j] = 0;
        }
    }
    
    CFRelease(audioFileURL);
    return audioData;
}

@end
