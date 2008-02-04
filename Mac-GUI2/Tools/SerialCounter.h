//
//  SerialCounter.h
//  R
//
//  Created by Simon Urbanek on 12/28/07.
//  Copyright 2007 Simon Urbanek. All rights reserved.
//
//  This is an implementation of a serial counter - i.e. a counter that is
//  increased periodically after a speicified time interval. Any number
//  of instances can be created, but note that each counter uses its own
//  thread.
//
//  In addition, this class also implements a single, application-wide,
//  global counter that can be used for low-overhead firing of events.
//  The counters are not necessarily precise, their main purpose is
//  usually buffer flushing and timeout signalling with very little
//  setup overhead.
//
//  The implementation is highly efficient on the enqueuing side, i.e.
//  is it ok to have a rapid rate of activations, usually many per
//  time slice.
//

#ifdef __OBJC__
#import <Foundation/Foundation.h>
#endif

typedef unsigned int counter_t;

// the counter can be accessed by any thread but must be treated read-only!
extern counter_t global_serial;

// use this flag in setGlobalWakeQueueSlot to specify that the selector
// should be fired on the main thread. If not specified, the selector will be
// run on the counter thread. IMPORTANT: if run on the counter thread,
// the action performed must be *very short* to not delay the counter.
// All longer actions must be performed on the main or dedicated thread.
#define WQF_ON_MAIN_THREAD 0x0001

// (the following is C code for efficiency)
// the slot can be any number between 0 and 15. all other numbers are ignored.
// the selector will be fired on the target with the specified user data
// roughly after interval seconds. The interval granularity depends on the
// interval of the global counter. Flags are OR-ed WQF_xx values.
// The intended use is to call setGlobalWakeQueueSlot for each activity
// and takes constant time and needs no locking.
BOOL setGlobalWakeQueueSlot(unsigned int slot, SEL selector, id target,
                            id userData, double interval, int flags);
// faster version of the above that assumes a pre-populated slot
// BOOL refreshGlobalWakeQueueSlot(unsigned int slot, double interval);
void resetGlobalWakeQueueSlot(unsigned int slot); 


#ifdef __OBJC__
// the ObjC counter class
@interface SerialCounter : NSObject {
    counter_t counter;
    double interval;
    NSThread *thread;
    BOOL active;
}

- (id) initWithInterval: (double) anInterval;
- (counter_t) value;

+ (SerialCounter*) counterWithInterval: (double) anInterval;

+ (SerialCounter*) sharedCounter;

+ (void) destroySharedCounter; // should be used ONLY on application exit during cleanup to release the global counter

@end
#endif
