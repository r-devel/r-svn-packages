//
//  RConsoleWinCtrl.h
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import "REngine.h"
#import "History.h"

@interface RConsoleWinCtrl : NSWindowController <REPLHandler, CocoaHandler> {
	/* UI components */
	IBOutlet NSTextView *textView;
	IBOutlet NSTextField *statusText;
	
	/* underlying R engine */
	REngine *engine;
	
	/* text editing sentinels */
	unsigned committedLength; // any text before this position cannot be edited by the user
    unsigned promptPosition;  // the last prompt is positioned at this position
	unsigned outputPosition;  // any output (stdxx or consWrite) is to be place here, if -1 then the text can be appended	

	/* console input history */
	History* hist;
	
	/* console font and its size */
	float currentFontSize;
	NSFont *textFont;

	/* quaue for incoming input */
	NSMutableArray *consoleInputQueue;
	NSString *currentConsoleInput;

	/* buffer for console output */
	char *writeBuffer;
	char *writeBufferPos;
	int  writeBufferLen;
	
	/* buffor used to talk to R via ReadConsole */
	char *readConsTransBuffer; // transfer buffer returned by handleReadConsole
	int readConsTransBufferSize; // size of the above buffer

}

- (void) writeConsoleDirectly: (NSString*) txt withColor: (NSColor*) color;

@end
