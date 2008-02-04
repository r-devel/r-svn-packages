//
//  RDocumentWinCtrl.m
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "RDocumentWinCtrl.h"
#import "RDocument.h"

@implementation RDocumentWinCtrl

- (void)windowDidLoad
{
	RDocument *doc = [self document];
	NSString *str = [doc string];
	if (str)
		[self replaceContentsWithString:str];
}

- (IBAction) executeSelection: (id) sender
{
	
}

- (NSString*) string
{
	return [textView string];
}

- (BOOL) isEdited
{
	return edited;
}

- (void) replaceContentsWithString: (NSString*) aString
{
	NSLog(@"RDocumentWinCtrl.replaceContentsWithString: (tv=%@)", textView);
	if (textView)
		[textView setString:aString];
	edited = NO;
}

- (void)textDidChange:(NSNotification *)aNotification
{
	edited = YES;
}

@end
