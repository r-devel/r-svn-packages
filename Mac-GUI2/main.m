//
//  main.m
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright __MyCompanyName__ 2007 . All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import "REngine/REngine.h"
#import <Rversion.h>
#import "RGUI.h"

NSString *Rapp_R_version_short;
NSString *Rapp_R_version;

int main(int argc, char *argv[])
{
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

	/* setup R version strings */
	Rapp_R_version_short = [[NSString alloc] initWithFormat:@"%d.%d", (R_VERSION >> 16), (R_VERSION >> 8)&255];
	Rapp_R_version = [[NSString alloc] initWithFormat:@"%s.%s", R_MAJOR, R_MINOR];

	SLog(@" - set APP VERSION (%s) and REVISION (%@)", R_GUI_VERSION_STR,
		 [[[NSBundle mainBundle] infoDictionary] objectForKey:@"CFBundleVersion"]);
	setenv("R_GUI_APP_VERSION", R_GUI_VERSION_STR, 1);
	setenv("R_GUI_APP_REVISION", [(NSString*)[[[NSBundle mainBundle] infoDictionary] objectForKey:@"CFBundleVersion"] UTF8String], 1);
	
	[NSApplication sharedApplication];
	[NSBundle loadNibNamed:@"MainMenu" owner:NSApp];

	/*
	SLog(@" - initalizing R");
	if (![[REngine mainEngine] activate]) {
		NSRunAlertPanel(NLS(@"Cannot start R"),[NSString stringWithFormat:NLS(@"Unable to start R: %@"), [[REngine mainEngine] lastError]],NLS(@"OK"),nil,nil);
		exit(-1);
	}*/
	
	SLog(@"main: finish launching");
	[NSApp finishLaunching];

	{ /* process all pending events before proceeding  to R's REPL */
		NSEvent *event;
		while (event = [NSApp nextEventMatchingMask:NSAnyEventMask
										  untilDate:nil
											 inMode:NSDefaultRunLoopMode 
											dequeue:YES])
			[NSApp sendEvent:event];
	}
	
	// ready to rock
	SLog(@"main: entering REPL");
	[[REngine mainEngine] runREPL];
	
	SLog(@"main: returned from REPL");
	[pool release];
	
	SLog(@"main: exiting with status 0");
	return 0;
}

/* return NSApplicationMain(argc, (const char **) argv); */
