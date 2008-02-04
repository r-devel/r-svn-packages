//
//  RDocumentWinCtrl.h
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface RDocumentWinCtrl : NSWindowController {
	IBOutlet NSTextView *textView;
	IBOutlet NSTextField *statusTextField;
	
	BOOL edited;
}

- (IBAction) executeSelection: (id) sender;

- (NSString*) string;
- (BOOL) isEdited;
- (void) replaceContentsWithString: (NSString*) aString;

@end
