//
//  MyDocument.h
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright __MyCompanyName__ 2007 . All rights reserved.
//


#import <Cocoa/Cocoa.h>

#define docTypeRSource @"R Source File";

@interface RDocument : NSDocument
{
	NSString *contents;
}

- (void) changeTitle: (NSString*) title;
- (void) setEditable: (BOOL) editable;

- (NSString*) string;

@end
