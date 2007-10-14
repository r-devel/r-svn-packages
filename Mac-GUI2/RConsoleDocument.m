//
//  RConsoleDocument.m
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "RConsoleDocument.h"
#import "RConsoleWinCtrl.h"
#import "RGUI.h"

@implementation RConsoleDocument

- (id)init
{
    self = [super init];
    if (self) {
		
        // Add your subclass-specific initialization here.
        // If an error occurs here, send a [self release] message and return nil.
		
    }
    return self;
}

- (void)makeWindowControllers
{
	SLog(@"RConsoleDocument: makeWindowControllers");
	RConsoleWinCtrl *wctrl = [[RConsoleWinCtrl alloc] initWithWindowNibName:@"RConsole"];
	[self addWindowController:wctrl];
}

- (void)windowControllerDidLoadNib:(NSWindowController *) aController
{
    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController has loaded the document's window.
	SLog(@"RConsoleDocument: windowControllerDidLoadNib");
}

- (NSData *)dataOfType:(NSString *)typeName error:(NSError **)outError
{
	SLog(@"RConsoleDocument: dataOfType: %@", typeName);
    // Insert code here to write your document to data of the specified type. If the given outError != NULL, ensure that you set *outError when returning nil.

    // You can also choose to override -fileWrapperOfType:error:, -writeToURL:ofType:error:, or -writeToURL:ofType:forSaveOperation:originalContentsURL:error: instead.

    // For applications targeted for Panther or earlier systems, you should use the deprecated API -dataRepresentationOfType:. In this case you can also choose to override -fileWrapperRepresentationOfType: or -writeToFile:ofType: instead.

    return nil;
}

- (BOOL)readFromData:(NSData *)data ofType:(NSString *)typeName error:(NSError **)outError
{
    // Insert code here to read your document from the given data of the specified type.  If the given outError != NULL, ensure that you set *outError when returning NO.

    // You can also choose to override -readFromFileWrapper:ofType:error: or -readFromURL:ofType:error: instead. 
    
    // For applications targeted for Panther or earlier systems, you should use the deprecated API -loadDataRepresentation:ofType. In this case you can also choose to override -readFromFile:ofType: or -loadFileWrapperRepresentation:ofType: instead.
    
    return YES;
}

@end
