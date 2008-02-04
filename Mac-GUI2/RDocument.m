//
//  MyDocument.m
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright __MyCompanyName__ 2007 . All rights reserved.
//

#import "RDocument.h"
#import "RDocumentWinCtrl.h"

#include <string.h>

@implementation RDocument

- (id)init
{
    self = [super init];
    if (self) {
		contents = nil;
    }
    return self;
}

- (void)makeWindowControllers
{
	NSLog(@"RDocument.makeWindowControllers");
	RDocumentWinCtrl *ctrl = [[RDocumentWinCtrl alloc] initWithWindowNibName:@"RDocument"];
	[self addWindowController: ctrl];
	[ctrl release];
}

- (NSString*) string
{
	return contents;
}

- (BOOL) isDocumentEdited
{
	NSArray *a = [self windowControllers];
	if ([a count] < 1) return NO;
	BOOL dirty = NO;
	for (RDocumentWinCtrl *ctrl in a)
		dirty |= [ctrl isEdited];	
	return dirty;
}

- (NSData *)dataOfType:(NSString *)typeName error:(NSError **)outError
{
	NSLog(@"RDocument.dataOfType:%@", typeName);
	NSString *str = contents;
	NSArray *a = [self windowControllers];
	if ([a count] > 0)
		str = [(RDocumentWinCtrl*)[a objectAtIndex:0] string];
	
	if (str) {
		const char *c = [str UTF8String];
		return [NSData dataWithBytes:c length:strlen(c)];
	}
    if ( outError != NULL ) {
		*outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
	}
	return nil;
}

- (BOOL)readFromData:(NSData *)data ofType:(NSString *)typeName error:(NSError **)outError
{
	NSLog(@"RDocument.readFromData:ofType:%@", typeName);
	
	NSString *str = [[NSString alloc] initWithData:data encoding: NSUTF8StringEncoding];
	if (!str) { // if UTF-8 fails, try Latin1
		str = [[NSString alloc] initWithData:data encoding: NSISOLatin1StringEncoding];
		if (!str) {
			if (outError)
				*outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:kTextUnsupportedEncodingErr userInfo:NULL];
			return NO;
		}
	}
	
	if (contents) [contents release];
	contents = str;
	NSArray * a = [self windowControllers];
	for (RDocumentWinCtrl *ctrl in a)
		[ctrl replaceContentsWithString: str];
	NSLog(@" - controllers = %d", [a count]);

    return YES;
}

- (void) changeTitle: (NSString *)title
{
}

- (void) setEditable: (BOOL) editable
{
	// FIXME: implement setEditable
}

@end
