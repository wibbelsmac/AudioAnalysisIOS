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
    
    UInt32 packetCount;
    UInt32 dataSize = sizeof(packetCount);
    result = AudioFileGetProperty(mAudioFile, kAudioFilePropertyAudioDataPacketCount, &dataSize, &packetCount);
    NSLog(@"File Opened, packet count is %d", (unsigned int)packetCount);
    
    UInt32 packetsRead = packetCount;
    UInt32 numBytesRead = -1;
    SInt16 *audioData = nil;
    if (packetCount > 0) {
        // Allocate Buffer
        audioData = (SInt16 *) malloc(2 * packetCount);
        // Read the Packets
        result = AudioFileReadPackets(mAudioFile, false, &numBytesRead, NULL, 0, &packetsRead, audioData);
        NSLog(@"Read %d bytes, %d packets", (unsigned int)numBytesRead, (unsigned int)packetsRead);
    }
        
    CFRelease(audioFileURL);
    return audioData;
}

@end
