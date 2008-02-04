//
//  RDocumentController.m
//  R
//
//  Created by Simon Urbanek on 2/4/08.
//  Copyright 2008 R-foundation. All rights reserved.
//

#import "RDocumentController.h"


@implementation RDocumentController

- (id) init
{
	self = [super init];
	if (self != nil) {
		NSLog(@"RDocumentController.init (self=%@)", self);
	}
	return self;
}


- (IBAction) newDocument: (id) sender
{
	NSArray *docs = [self documents];
	NSLog(@"RDocumentController.newDocument: docs.count=%d", [docs count]);
	// if no documents are open, open the console; resort to default behavior otherwise
	NSError *err = nil;
	NSDocument * doc;
	if ([docs count] == 0)
		doc = [self makeUntitledDocumentOfType:@"R Console" error:&err];
	else
		doc = [self makeUntitledDocumentOfType:@"R Source File" error:&err];
	if (err)
		NSLog(@" - FAILED: %@", err);
	if (doc) {
		NSLog(@" - doc = %@", doc);
		[self addDocument:doc];
		[doc makeWindowControllers];
		[doc showWindows];
	}	
}

#pragma mark --- application delegate methods ---

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
	[self setAutosavingDelay:4*59.0]; // enable autosave, ca. 4 min
}

- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)theApplication
{
	return YES;
}

@end
