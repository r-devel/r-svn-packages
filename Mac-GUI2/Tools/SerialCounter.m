//
//  SerialCounter.m
//  R
//
//  Created by Simon Urbanek on 12/28/07.
//  Copyright 2007 Simon Urbanek. All rights reserved.
//

#import "SerialCounter.h"

double defaultSharedInterval = 0.050; // default global slice 50ms (20Hz)

counter_t global_serial;
static double global_interval;
static SerialCounter* sharedCounter = nil;

#define WQ_SLOTS 16
static id wq_target[WQ_SLOTS];
static id wq_data[WQ_SLOTS];
static SEL wq_sel[WQ_SLOTS];
static counter_t wq_wake[WQ_SLOTS];
static int wq_flags[WQ_SLOTS];
static unsigned int wq_mask;

BOOL setGlobalWakeQueueSlot(unsigned int slot, SEL selector, id target, id userData, double interval, int flags)
{
    if (slot >= WQ_SLOTS) return NO;
    wq_wake[slot] = 0; // first set wake to 0 for thread-safety
    wq_data[slot] = userData;
    wq_target[slot] = target;
    wq_sel[slot] = selector;
    wq_flags[slot] = flags;
    wq_wake[slot] = global_serial + (int)((interval/global_interval)+0.5); // now set the wake
    wq_mask |= 1 << slot;
    if (!sharedCounter) [SerialCounter sharedCounter]; // start the global counter if not present
    return YES;
}

void resetGlobalWakeQueueSlot(unsigned int slot)
{
    if (slot >= WQ_SLOTS) return;
    wq_mask &= 0xffff ^ (1 << slot);
    wq_wake[slot] = 0;
    wq_target[slot] = nil;
    return;
}

@implementation SerialCounter

- (counter_t) value
{
    return counter;
}

- (void) run: (id) userData
{
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    thread = [NSThread currentThread];
    while (active) {
        [NSThread sleepForTimeInterval:interval];
        counter++;
    }
    [pool release];
}

- (void) runShared: (id) userData
{
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    thread = [NSThread currentThread];
    [thread setName:@"global serial counter thread"];
    while (active) {
        [NSThread sleepForTimeInterval:interval];
        counter++;
        global_serial++;
        if (wq_mask) { // check for any pending wake-ups
            int i = 0;
            while (i < WQ_SLOTS) {
                counter_t w = wq_wake[i];
                if (w > 0 && w <= global_serial) {
                    id target = wq_target[i];
                    id data = wq_data[i];
                    SEL sel = wq_sel[i];
                    int flags = wq_flags[i];
                    if (wq_wake[i] == w) { // still same wake? for thread-safety
                        wq_wake[i] = 0; // remove from the queue
                        if (flags & WQF_ON_MAIN_THREAD)
                            [target performSelectorOnMainThread:sel withObject:data waitUntilDone:NO];
                        else // FIXME: we hope that the method doesn't take too long to screw our timer
                            [target performSelector:sel withObject:data];
                    }
                }
                i++;
            }
        }
    }
    thread = nil;
    [pool release];
}

- (void) dealloc
{
    if (active) {
        active = NO;
        thread = nil;
    }
    [super dealloc];
}


- (id) initWithInterval: (double) anInterval
{
    self = [super init];
    if (self) {
        counter = 0;
        interval = anInterval;
        thread = nil;
        active = YES;
        [NSThread detachNewThreadSelector:@selector(run:) toTarget:self withObject:nil];
    }
    return self;
}

// private initializer for the shared counter
- (id) initSharedWithInterval: (double) anInterval
{
    if (sharedCounter) {
        [super release];
        return sharedCounter;
    }
    self = [super init];
    if (self) {
        counter = 0;
        interval = anInterval;
        thread = nil;
        active = YES;
        sharedCounter = self;
        [NSThread detachNewThreadSelector:@selector(runShared:) toTarget:self withObject:nil];
    }
    return self;
}

+ (SerialCounter*) counterWithInterval: (double) anInterval
{
    SerialCounter* counter = [[SerialCounter alloc] initWithInterval:anInterval];
    return counter?[counter autorelease]:nil;
}

+ (SerialCounter*) sharedCounter
{
    if (sharedCounter) return sharedCounter;
    sharedCounter = [[SerialCounter alloc] initSharedWithInterval:defaultSharedInterval];
    return sharedCounter;
}

+ (void) destroySharedCounter
{
    if (sharedCounter) {
        [sharedCounter release];
        sharedCounter = nil;
    }
}

@end
