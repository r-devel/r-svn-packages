//
//  RConsoleWinCtrl.m
//  R
//
//  Created by Simon Urbanek on 10/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "RGUI.h"
#import "RConsoleWinCtrl.h"

#import "Preferences.h"
#import "RDocument.h"

// size of the console output cache buffer
#define DEFAULT_WRITE_BUFFER_SIZE 32768
// high water-mark of the buffer - it's [length - x] where x is the smallest possible size to be flushed before a new string will be split.
#define writeBufferHighWaterMark  (DEFAULT_WRITE_BUFFER_SIZE-4096)
// low water-mark of the buffer - if less than the water mark is available then the buffer will be flushed
#define writeBufferLowWaterMark   2048

// class-level colors
static NSColor *outputColor, *promptColor, *inputColor;

@implementation RConsoleWinCtrl

#pragma mark --- Initialization ---

+ (void) initialize
{
	NSLog(@"RConsoleWinCtrl.initialize (allocate default colors)");
	outputColor = [[NSColor colorWithCalibratedRed: 0.0 green: 0.0 blue: 0.5 alpha: 1.0] retain];
	promptColor = [[NSColor colorWithCalibratedRed: 0.5 green: 0.0 blue: 0.5 alpha: 1.0] retain];
	inputColor  = [[NSColor colorWithCalibratedRed: 0.0 green: 0.0 blue: 0.0 alpha: 1.0] retain];
}

- (id)initWithWindowNibName:(NSString *)windowNibName
{
	self = [super initWithWindowNibName:windowNibName];
	if (self) {
		outputPosition = promptPosition = committedLength = 0;
		
		consoleInputQueue = [[NSMutableArray alloc] initWithCapacity:8];
		currentConsoleInput = nil;
		
		writeBufferLen = DEFAULT_WRITE_BUFFER_SIZE;
		writeBufferPos = writeBuffer = (char*) malloc(writeBufferLen);
		
		readConsTransBufferSize = 1024; // initial size - will grow as needed
		readConsTransBuffer = (char*) malloc(readConsTransBufferSize);
		
		hist = [[History alloc] init];

		currentFontSize = [Preferences floatForKey: kConsoleFontSize withDefault: 11.0];
		textFont = [[NSFont userFixedPitchFontOfSize:currentFontSize] retain];
	}
	NSLog(@"RConsoleWinCtrl%@: initWithWindowNibName: %@", self, windowNibName);
	return self;
}

/* initializes R engine used by this console instance, sets engine cv accordingly */
- (void) initEngine
{
	SLog(@" - init R_LIBS");	
	{
		NSString *prefStr = [Preferences stringForKey:kAddUserLibPath withDefault:nil];
		BOOL flag = !isAdmin(); // the default is YES for users and NO for admins
		if (prefStr)
			flag=[prefStr isEqualToString: @"YES"];
		if (flag) {
			char *cRLIBS = getenv("R_LIBS");
			NSString *addPath = [[NSString stringWithFormat:@"~/Library/R/%@/library", Rapp_R_version_short] stringByExpandingTildeInPath];
			if (![[NSFileManager defaultManager] fileExistsAtPath:addPath]) { // make sure the directory exists
				[[NSFileManager defaultManager] createDirectoryAtPath:[@"~/Library/R" stringByExpandingTildeInPath] attributes:nil];
				[[NSFileManager defaultManager] createDirectoryAtPath:[[NSString stringWithFormat:@"~/Library/R/%@", Rapp_R_version_short] stringByExpandingTildeInPath] attributes:nil];
				[[NSFileManager defaultManager] createDirectoryAtPath:addPath attributes:nil];
			}
			if (cRLIBS && *cRLIBS)
				addPath = [NSString stringWithFormat: @"%s:%@", cRLIBS, addPath];
			setenv("R_LIBS", [addPath UTF8String], 1);
			SLog(@" - setting R_LIBS=%s", [addPath UTF8String]);
		}
	}
	
	SLog(@" - starting REngine");
	{
		char *args[5]={ "R", "--no-save", "--no-restore-data", "--gui=cocoa", 0 };
		SLog(@" - initialize REngine");
		engine = [[REngine alloc] initWithHandler:self arguments:args];
		[engine setCocoaHandler:self];
		
		SLog(@" - activate R");
		if (![[REngine mainEngine] activate]) {
			NSRunAlertPanel(NLS(@"Cannot start R"),[NSString stringWithFormat:NLS(@"Unable to start R: %@"), [[REngine mainEngine] lastError]],NLS(@"OK"),nil,nil);
			exit(-1);
		}
            // start R on a separate thread
            SLog(@" - start REPL thread");
            [NSThread detachNewThreadSelector:@selector(run:) toTarget:engine withObject:self];
	}
}


#pragma mark --- NSWindowController delegate methods ---

- (void)windowDidLoad
{
	SLog(@"RConsoleWinCtrl: windowDidLoad, name=%@, window=%@", [self windowNibName], [self window]);
	SLog(@" - font = %@", textFont);
	[textView setFont:textFont];
	//[textView setTypingAttributes:[NSDictionary dictionaryWithObject:textFont forKey:@"NSFont"]];
	[textView setString:@"\n"];
	[self initEngine];
}

- (NSString *)windowTitleForDocumentDisplayName:(NSString *)displayName
{
	// the title is always "R Console" regardless of the file name
	return @"R Console";
}


#pragma mark --- console API ---

/* console input - the string passed here is handled as if it was typed on the console */
- (void) consoleInput: (NSString*) cmd interactive: (BOOL) inter
{
	//@synchronized(textViewSync) {
	if (!inter) {
		int textLength = [[textView textStorage] length];
		if (textLength>committedLength)
			[textView replaceCharactersInRange:NSMakeRange(committedLength,textLength-committedLength) withString:@""];
		[textView setSelectedRange:NSMakeRange(committedLength,0)];
		[textView insertText: cmd];
		textLength = [[textView textStorage] length];
		[textView setTextColor:inputColor range:NSMakeRange(committedLength,textLength-committedLength)];
	}
		
	if (inter) {
		if ([cmd characterAtIndex:[cmd length]-1]!='\n') cmd=[cmd stringByAppendingString: @"\n"];
		[consoleInputQueue addObject:[[NSString alloc] initWithString:cmd]];
		[self handleShowInfo:@""]; //  FIXME: possibly don't use handleShowInfo
	}
	//}
}

- (BOOL) processSingleEventBlocking: (BOOL) blocking
{
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	NSEvent *event = [NSApp nextEventMatchingMask:NSAnyEventMask
										untilDate:blocking?[NSDate distantFuture]:nil
										   inMode:NSDefaultRunLoopMode 
										  dequeue:YES];
	if (event)
		[NSApp sendEvent:event];
	[pool release];
	return event?YES:NO;
}

/* this writes R output to the Console window directly, i.e. without using a buffer. Use handleWriteConsole: for the regular way. */
- (void) writeConsoleDirectly: (NSString*) txt withColor: (NSColor*) color {
	//@synchronized(textViewSync) {
	NSTextStorage *textStorage = [textView textStorage];
	NSRange origSel = [textView selectedRange];
	unsigned tl = [txt length];
	if (tl>0) {
		unsigned oldCL=committedLength;
		/* NSLog(@"original: %d:%d, insertion: %d, length: %d, prompt: %d, commit: %d", origSel.location,
		 origSel.length, outputPosition, tl, promptPosition, committedLength); */
		SLog(@"RConsoleWinCtrl writeConsoleDirectly, beginEditing");
		[textStorage beginEditing];
		committedLength=0;
		[textStorage replaceCharactersInRange: NSMakeRange(outputPosition,0) withString: txt];
		[textStorage addAttribute:@"NSColor" value:color range: NSMakeRange(outputPosition, tl)];
		//[textStorage addAttribute:@"NSFont" value:[[RController sharedController] currentFont] range: NSMakeRange(index, [text length])];
		if (outputPosition<=promptPosition) promptPosition+=tl;
		committedLength=oldCL;
		if (outputPosition<=committedLength) committedLength+=tl;
		if (outputPosition<=origSel.location) origSel.location+=tl;
		outputPosition+=tl;
		[textStorage endEditing];
		SLog(@"RConsoleWinCtrl writeConsoleDirectly, endEditing");
		[textView setSelectedRange:origSel];
		[textView scrollRangeToVisible:origSel];
	}
	//}
}


#pragma mark --- REngine callbacks ---
/* --- REngine callbacks --- */

- (int) handleChooseFile:(char *)buf len:(int)len isNew:(int)isNew
{
	const char *fn;
	int answer;
	NSSavePanel *sp;
	NSOpenPanel *op;
	
	*buf = 0;
	if(isNew==1){
		sp = [NSSavePanel savePanel];
		[sp setTitle:NLS(@"Choose New File Name")];
		answer = [sp runModalForDirectory:nil file:nil];
		
		if(answer == NSOKButton) {
			if([sp filename] != nil){
				fn = [[sp filename] UTF8String];
				if (strlen(fn)>=len) {
					SLog(@"** handleChooseFile: bufer too small, truncating");
					memcpy(buf, fn, len-1);
					buf[len-1]=0;
				} else
					strcpy(buf, fn);
			}
		}
	} else {
		op = [NSOpenPanel openPanel];
		[op setTitle:NLS(@"Choose File")];
		answer = [op runModalForDirectory:nil file:nil];
		
		if(answer == NSOKButton) {
			if([op filename] != nil){
				fn = [[op filename] UTF8String];
				if (strlen(fn)>=len) {
					SLog(@"** handleChooseFile: bufer too small, truncating");
					memcpy(buf, fn, len-1);
					buf[len-1]=0;
				} else
					strcpy(buf, fn);
			}
		}
	}
	[[self window] makeKeyWindow];
	return strlen(buf); // is is used? it's potentially incorrect...
}

- (void) flushConsole {
	if (writeBuffer!=writeBufferPos) {
		[self writeConsoleDirectly:[NSString stringWithUTF8String:writeBuffer] withColor:outputColor];
		writeBufferPos=writeBuffer;
	}
}

- (int) handleEdit: (char*) file
{
	//NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	NSString *fn = [NSString stringWithUTF8String:file];
	if (fn) fn = [fn stringByExpandingTildeInPath];
	if (!fn) Rf_error("Invalid file name.");
	
	SLog(@"RController.handleEdit: %s", file);
	
	if (![[NSFileManager defaultManager] isReadableFileAtPath:fn]) {
		//[pool release];
		return 0;
	}
	NSURL *url = [[NSURL alloc] initFileURLWithPath:fn];
	NSError *theError;
	RDocument *document = [[NSDocumentController sharedDocumentController] openDocumentWithContentsOfURL:url display:YES error:&theError];
	[document setREditFlag: YES];
	
	NSArray *wcs = [document windowControllers];
	if ([wcs count]<1) {
		SLog(@"handleEdit: WARNING, no window controllers for newly created document!");
	} else {
		NSWindowController *wc = (NSWindowController*)[wcs objectAtIndex:0];
		NSWindow *win = [wc window];
		if (win) [NSApp runModalForWindow:win];
		else { SLog(@"handleEdit: WARNING, window is null!"); };
		if ([wcs count]>1) {
			SLog(@"handleEdit: WARNING, there is more than one window controller, ignoring all but the first one.");
		}
	}
	
	//[pool release];
	return(0);
}

/* FIXME: the filename is not set for newly created files */
- (int) handleEditFiles: (int) nfile withNames: (char**) file titles: (char**) wtitle pager: (char*) pager
{
	int    	i;
    
	SLog(@"RController.handleEditFiles (%d of them, pager %s)", nfile, pager);
    if (nfile <=0) return 1;
	
    for (i = 0; i < nfile; i++) {
		NSString *fn = [NSString stringWithUTF8String:file[i]];
		if (fn) fn = [fn stringByExpandingTildeInPath];
		if (!fn) {
			if([[NSFileManager defaultManager] fileExistsAtPath:fn]) {
				NSURL *url = [[NSURL alloc] initFileURLWithPath:fn];
				NSError *theError;
				[[NSDocumentController sharedDocumentController] openDocumentWithContentsOfURL:url display:YES error:&theError];
			} else
				[[NSDocumentController sharedDocumentController] newDocument: self];
			
			NSDocument *document = [[NSDocumentController sharedDocumentController] currentDocument];
			if(wtitle[i]!=nil)
				[(RDocument*)document changeTitle: [NSString stringWithUTF8String:wtitle[i]]];
		}
    }
	return 1;
}

- (void) handleFlushConsole {
	[self flushConsole];
	//[self flushStdConsole];
}

/* this writes R output to the Console window, but indirectly by using a buffer */
- (void) handleWriteConsole: (NSString*) txt {
	if (!txt) return;
	const char *s = [txt UTF8String];
	int sl = strlen(s);
	int fits = writeBufferLen-(writeBufferPos-writeBuffer)-1;
	
	// let's flush the buffer if the new string is large and it would, but the buffer should be occupied
	if (fits<sl && fits>writeBufferHighWaterMark) {
		// for efficiency we're not using handleFlushConsole, because that would trigger stdxx flush, too
		[self writeConsoleDirectly:[NSString stringWithUTF8String:writeBuffer] withColor:outputColor];
		writeBufferPos=writeBuffer;
		fits = writeBufferLen-1;
	}
	
	while (fits<sl) {	// ok, we're in a situation where we must split the string
		memcpy(writeBufferPos, s, fits);
		writeBufferPos[writeBufferLen-1]=0;
		[self writeConsoleDirectly:[NSString stringWithUTF8String:writeBuffer] withColor:outputColor];
		sl-=fits; s+=fits;
		writeBufferPos=writeBuffer;
		fits=writeBufferLen-1;
	}
	
	strcpy(writeBufferPos, s);
	writeBufferPos+=sl;
	
	// flush the buffer if the low watermark is reached
	if (fits-sl<writeBufferLowWaterMark) {
		[self writeConsoleDirectly:[NSString stringWithUTF8String:writeBuffer] withColor:outputColor];
		writeBufferPos=writeBuffer;
	}
}

/* Just writes the prompt in a different color */
- (void)handleWritePrompt: (NSString*) prompt {
    [self handleFlushConsole];
	//@synchronized(textViewSync) {
	NSTextStorage *textStorage = [textView textStorage];
	unsigned textLength = [textStorage length];
	int promptLength=[prompt length];
	//		NSLog(@"Prompt: %@", prompt);
	NSRange lr = [[textStorage string] lineRangeForRange:NSMakeRange(textLength,0)];
	SLog(@"RController handleWritePrompt: '%@', beginEditing", prompt);
	[textStorage beginEditing];
	promptPosition=textLength;
	if (lr.location!=textLength) { // the prompt must be on the beginning of the line
		//[textStorage insertText: @"\n" atIndex: textLength withColor:[consoleColors objectAtIndex:iPromptColor]];
		[textStorage replaceCharactersInRange: NSMakeRange(promptPosition,0) withString: @"\n"];
		[textStorage addAttribute:@"NSColor" value:promptColor range: NSMakeRange(promptPosition, 1)];
		textLength = [textStorage length];
		promptLength++;
	}
		
	if (promptLength>0) {
		//[textStorage insertText:prompt atIndex: textLength withColor:[consoleColors objectAtIndex:iPromptColor]];
		[textStorage replaceCharactersInRange: NSMakeRange(promptPosition,0) withString: prompt];
		[textStorage addAttribute:@"NSColor" value:promptColor range: NSMakeRange(promptPosition, promptLength)];
		if (promptLength>1) // this is a trick to make sure that the insertion color doesn't change at the prompt
			[textStorage addAttribute:@"NSColor" value:inputColor range:NSMakeRange(promptPosition+promptLength-1, 1)];
		committedLength=promptPosition+promptLength;
	}
	committedLength=promptPosition+promptLength;
	[textStorage endEditing];
	SLog(@"RController handleWritePrompt: '%@', endEditing", prompt);
	{
		NSRange targetRange = NSMakeRange(committedLength,0);
		[textView setSelectedRange:targetRange];
		[textView scrollRangeToVisible:targetRange];
	}
	//}
}

- (void) handleProcessEvents
{
	while ([self processSingleEventBlocking:NO]) {};
}

- (void)  handleProcessingInput: (char*) cmd
{
	NSString *s = [[NSString alloc] initWithUTF8String:cmd];
	
	//@synchronized(textViewSync) {
	unsigned textLength = [[textView textStorage] length];
		
	[textView setSelectedRange:NSMakeRange(committedLength, textLength-committedLength)];
	[textView insertText:s];
	textLength = [[textView textStorage] length];
	[textView setTextColor:inputColor range:NSMakeRange(committedLength, textLength-committedLength)];
	outputPosition=committedLength=textLength;
		
	// remove undo actions to prevent undo across prompts
	[[textView undoManager] removeAllActions];
	//}
	
	[s release];	
}

- (char*) handleReadConsole: (int) addtohist
{
	if (currentConsoleInput) {
		[currentConsoleInput release];
		currentConsoleInput=nil;
	}
	
	while ([consoleInputQueue count]==0)
		[self processSingleEventBlocking: YES];
	
	currentConsoleInput = [consoleInputQueue objectAtIndex:0];
	[consoleInputQueue removeObjectAtIndex:0];
	
	if (addtohist) {
		//		Figure out how to get hold of ParseStatus here!
		[hist commit:currentConsoleInput];
		// FIXME: add historyView support
		//[historyView reloadData];
	}
	
	{
		const char *c = [currentConsoleInput UTF8String];
		if (!c) return 0;
		if (strlen(c)>readConsTransBufferSize-1) { // grow as necessary
			free(readConsTransBuffer);
			readConsTransBufferSize = (strlen(c)+2048)&0xfffffc00;
			readConsTransBuffer = (char*) malloc(readConsTransBufferSize);
		} // we don't shrink the buffer if gets too large - we may want to think about that ...
		
		strcpy(readConsTransBuffer, c);
	}
	return readConsTransBuffer;
}

- (int) handleSystemCommand: (char*) cmd
{	
	int cstat=-1;
	pid_t pid;
	
	// FIXME: enable root system commands again
	/*
	if ([self getRootFlag]) {
		FILE *f;
		char *argv[3] = { "-c", cmd, 0 };
		int res;
 		NSBundle *b = [NSBundle mainBundle];
		char *sushPath=0;
		if (b) {
			NSString *sush=[[b resourcePath] stringByAppendingString:@"/sush"];
			sushPath = (char*) malloc([sush cStringLength]+1);
			[sush getCString:sushPath maxLength:[sush cStringLength]];
		}
		
		res = runRootScript(sushPath?sushPath:"/bin/sh",argv,&f,1);
		if (!res && f) {		
			int fd = fileno(f);
			if (fd != -1) {
				struct timespec peSleep = { 0, 50000000 }; // 50ms sleep
				[self setRootFD:fileno(f)];
				
				while (rootFD!=-1) { // readThread will reset rootFD to -1 when reaching EOF
					nanosleep(&peSleep, 0); // sleep at least 50ms between PE calls (they're expensive)
					Re_ProcessEvents();
				}
			}
		}
		if (sushPath) free(sushPath);
		return res;
	}*/
	
	pid=fork();
	if (pid==0) {
		// int sr;
		// reset signal handlers
		signal(SIGINT, SIG_DFL);
		signal(SIGTERM, SIG_DFL);
		signal(SIGQUIT, SIG_DFL);
		signal(SIGALRM, SIG_DFL);
		signal(SIGCHLD, SIG_DFL);
		execl("/bin/sh", "/bin/sh", "-c", cmd, NULL);
		exit(-1);
		//sr=system(cmd);
		//exit(WEXITSTATUS(sr));
	}
	if (pid==-1) return -1;
	
	{
		struct timespec peSleep = { 0, 100000000 }; // 100ms sleep
		while (1) {
			pid_t w = waitpid(pid, &cstat, WNOHANG);
			if (w!=0) break;
			nanosleep(&peSleep, 0); // sleep at least 50ms between PE calls (they're expensive)
			Re_ProcessEvents();
		}
	}
	// FIXME: [[RController sharedController] rmChildProcess: pid];
	return cstat;
}	



- (int) handleShowFiles: (int) nfile withNames: (char**) file headers: (char**) headers windowTitle: (char*) wtitle pager: (char*) pages andDelete: (BOOL) del
{
	int    	i;
    
    if (nfile <=0) return 1;
	SLog(@"RController.handleShowFiles (%d of them, title %s, pager %s)", nfile, wtitle, pages);
	
    for (i = 0; i < nfile; i++){
		NSString *fn = [NSString stringWithUTF8String:file[i]];
		if (fn) fn =[fn stringByExpandingTildeInPath];
		if (fn) {
			NSURL *url = [[NSURL alloc] initFileURLWithPath:fn];
			NSError *theError;
			RDocument *document = [[NSDocumentController sharedDocumentController] openDocumentWithContentsOfURL:url display:YES error:&theError];
			// don't display - we need to prevent the window controller from using highlighting
			if (document) {
				NSArray *wcs = [document windowControllers];
				if (wcs && [wcs count]>0) {
					SLog(@" - Disabling syntax highlighting for this document");
					// FIXME: add setPlain
					//[(RDocumentWinCtrl*) [wcs objectAtIndex:0] setPlain:YES];
				}
				if (wtitle)
					[document changeTitle: [NSString stringWithUTF8String:wtitle]];
				[document setEditable: NO];
				SLog(@" - finally show the document window");
				[document showWindows];
			}
		}
    }
	return 1;
}

- (void) handleBusy: (BOOL) isBusy {
	/*
    if (isBusy)
        [progressWheel startAnimation:self];
    else
        [progressWheel stopAnimation:self];
	
	busyRFlag = isBusy;
	if (toolbarStopItem) {
		if (isBusy || childPID>0)
			[toolbarStopItem setEnabled:YES];
		else
			[toolbarStopItem setEnabled:NO];
	}
	 */
}

- (void)  handleShowMessage: (char*) msg
{
	NSRunAlertPanel(NLS(@"R Message"),[NSString stringWithUTF8String:msg],NLS(@"OK"),nil,nil);
}

- (void)  handleShowInfo: (NSString*) msg
{
	[statusText setStringValue:msg?msg:@""];
}

#pragma mark --- Console text view delegate methods ---

/* Allow changes only for uncommitted text */
- (BOOL)textView:(NSTextView *)aTextView shouldChangeTextInRange:(NSRange)affectedCharRange replacementString:(NSString *)replacementString {
	if (aTextView != textView) return YES; /* this should never happen */
	if (replacementString && /* on font change we get nil replacementString which is ok to pass through */
		affectedCharRange.location < committedLength) { /* if the insertion is outside editable scope, append at the end */
		[textView setSelectedRange:NSMakeRange([[textView textStorage] length],0)];
		[textView insertText:replacementString];
		return NO;
	}
	return YES;
}

- (BOOL)textView:(NSTextView *)aTextView doCommandBySelector:(SEL)commandSelector {
    BOOL retval = NO;

	if (aTextView != textView) return NO; /* this should never happen */

	SLog(@"RConsoleWinCtrl textView: doCommandBySelector: %@\n", NSStringFromSelector(commandSelector));
	
    if (@selector(insertNewline:) == commandSelector) {
        unsigned textLength = [[textView textStorage] length];
		[textView setSelectedRange:NSMakeRange(textLength,0)];
        if (textLength >= committedLength) {
			[textView insertText:@"\n"];
			textLength = [[textView textStorage] length];
			[self consoleInput: [[textView attributedSubstringFromRange:NSMakeRange(committedLength, textLength - committedLength)] string] interactive: YES];
			return(YES);
        }
        retval = YES;
    }
	
	// ---- history browsing ----
	if (@selector(moveUp:) == commandSelector) {
        unsigned textLength = [[textView textStorage] length];        
        NSRange sr=[textView selectedRange];
        if (sr.location==committedLength || sr.location==textLength) {
            NSRange rr=NSMakeRange(committedLength, textLength-committedLength);
            NSString *text = [[textView attributedSubstringFromRange:rr] string];
            if ([hist isDirty]) {
                [hist updateDirty: text];
            }
            NSString *news = [hist prev];
            if (news!=nil) {
                [news retain];
                sr.length=0; sr.location=committedLength;
                [textView setSelectedRange:sr];
                [textView replaceCharactersInRange:rr withString:news];
                [textView insertText:@""];
                [news release];
            }
            retval = YES;
        }
    }
    if (@selector(moveDown:) == commandSelector) {
        unsigned textLength = [[textView textStorage] length];        
        NSRange sr=[textView selectedRange];
        if ((sr.location==committedLength || sr.location==textLength) && ![hist isDirty]) {
            NSRange rr=NSMakeRange(committedLength, textLength-committedLength);
            NSString *news = [hist next];
            if (news==nil) news=@""; else [news retain];
            sr.length=0; sr.location=committedLength;
            [textView setSelectedRange:sr];
            [textView replaceCharactersInRange:rr withString:news];
            [textView insertText:@""];
            [news release];
            retval = YES;
        }
    }
    
	// ---- make sure the user won't accidentally get out of the input line ----
	
	if (@selector(moveToBeginningOfParagraph:) == commandSelector || @selector(moveToBeginningOfLine:) == commandSelector) {
        [textView setSelectedRange: NSMakeRange(committedLength,0)];
        retval = YES;
    }
	
	if (@selector(moveToBeginningOfParagraphAndModifySelection:) == commandSelector || @selector(moveToBeginningOfLineAndModifySelection:) == commandSelector) {
		// FIXME: this kills the selection - we should retain it ...
        [textView setSelectedRange: NSMakeRange(committedLength,0)];
        retval = YES;
    }
	
	if (@selector(moveWordLeft:) == commandSelector || @selector(moveLeft:) == commandSelector ||
		@selector(moveWordLeftAndModifySelection:) == commandSelector || @selector(moveLeftAndModifySelection:) == commandSelector) {
        NSRange sr=[textView selectedRange];
		if (sr.location==committedLength) return YES;
	}
	
	// ---- code/file completion ----
	
	if (@selector(insertTab:) == commandSelector) {
		[textView complete:self];
		retval = YES;
	}
	
	// ---- cancel ---
	
	if (@selector(cancel:) == commandSelector || @selector(cancelOperation:) == commandSelector) {
		[self breakR:self];
		retval = YES;
	}
    
	return retval;
}


@end
