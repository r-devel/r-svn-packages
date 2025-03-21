/*
 *  R.app : a Cocoa front end to: "R A Computer Language for Statistical Data Analysis"
 *  
 *  R.app Copyright notes:
 *                     Copyright (C) 2004-5  The R Foundation
 *                     written by Stefano M. Iacus and Simon Urbanek
 *
 *                  
 *  R Copyright notes:
 *                     Copyright (C) 1995-1996   Robert Gentleman and Ross Ihaka
 *                     Copyright (C) 1998-2001   The R Development Core Team
 *                     Copyright (C) 2002-2004   The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 59 Temple Place,
 *  Suite 330, Boston, MA  02111-1307  USA.
 *
 *  Created by Simon Urbanek on Wed Dec 10 2003.
 *
 */

#import <Cocoa/Cocoa.h>
#include "Rinit.h"
#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
#import "REngine.h"

#if R_VERSION < R_Version(2,5,0)
#define RS_ParseVector R_ParseVector
#else
#define RS_ParseVector(A,B,C) R_ParseVector(A,B,C,R_NilValue)
#endif

#define DO_RENG_EVAL_STATUS(S) if (cocoaHandler) [cocoaHandler handleShowInfo:[NSString stringWithFormat:@"executing: %@", S]];
#define DONE_RENG_EVAL_STATUS() if (cocoaHandler) [cocoaHandler handleShowInfo:nil];

/* this is also provided in RGUI.h, but we want to be independent */
#ifndef SLog
#if defined DEBUG_RGUI && defined PLAIN_STDERR
#define SLog(X,...) NSLog(X, ## __VA_ARGS__)
#else
#define SLog(X,...)
#endif
#endif

static REngine* mainRengine=nil;

// this flag causes some parts of the code to not use REngine if that would cause re-entrance
// it is meant for the user-level code, not for REngine itself - such that the UI can react and display appropriate warnings
BOOL preventReentrance = NO;

@implementation REngine

+ (REngine*) mainEngine
{
    if (mainRengine==nil)
        mainRengine=[[REngine alloc] init];
    return mainRengine;
}

+ (id <REPLHandler>) mainHandler
{
	return [mainRengine handler];
}

+ (id <CocoaHandler>) cocoaHandler
{
	return [mainRengine cocoaHandler];
}

- (id) init
{
    return [self initWithHandler:nil];
}

- (id) initWithHandler: (id <REPLHandler>) hand
{
    char *args[4]={ "R", "--no-save", "--gui=cocoa", 0 };
	return [self initWithHandler: hand arguments: args];
}

- (id) initWithHandler: (id <REPLHandler>) hand arguments: (char**) args
{
	int i=0;
	argc=0;
	while (args[argc]) argc++;
	
	argv = (char**) malloc(sizeof(char*) * (argc+1));
	while (i<argc) {
		argv[i]=(char*) malloc(strlen(args[i])+1);
		strcpy(argv[i], args[i]);
		i++;
	}
	argv[i]=0;
	
    replHandler=hand;
	cocoaHandler=nil; // cocoaHandlier is optional
    mainRengine=self;
    loopRunning=NO;
	active=NO;
	insideR=0;
	maskEvents=0;
	saveAction=@"ask";
	
	/* --- set R_HOME --- */
	SLog(@" - set R home");
	if (!getenv("R_HOME")) {
		NSBundle *rfb = [NSBundle bundleWithIdentifier:@"org.r-project.R-framework"];
		if (!rfb) {
			SLog(@" * problem: R_HOME is not set and I can't find the framework bundle");
			if ([[NSFileManager defaultManager] fileExistsAtPath:@"/Library/Frameworks/R.framework/Resources/bin/R"]) {
				SLog(@" * I'm being desperate and I found R at /Library/Frameworks/R.framework - so I'll use it, wish me luck");
				setenv("R_HOME", "/Library/Frameworks/R.framework/Resources", 1);
			} else
				SLog(@" * I didn't even found R framework in the default location, I'm giving up - you're on your own");
		} else {
			SLog(@"   %s", [[rfb resourcePath] UTF8String]);
			setenv("R_HOME", [[rfb resourcePath] UTF8String], 1);
		}
	}
	
	/* --- set R_xxx_DIR --- */
	{
		char tp[1024];
		/* since 2.2.0 those are set in the R shell script, so we need to set them as well */
		/* FIXME: possible buffer-overflow attack by over-long R_HOME */
		if (!getenv("R_INCLUDE_DIR")) {
			strcpy(tp, getenv("R_HOME")); strcat(tp, "/include"); setenv("R_INCLUDE_DIR", tp, 1);
		}
		if (!getenv("R_SHARE_DIR")) {
			strcpy(tp, getenv("R_HOME")); strcat(tp, "/share"); setenv("R_SHARE_DIR", tp, 1);
		}
		if (!getenv("R_DOC_DIR")) {
			strcpy(tp, getenv("R_HOME")); strcat(tp, "/doc"); setenv("R_DOC_DIR", tp, 1);
		}
	}
	
	/* --- set R_ARCH --- */
#ifdef __ppc__
#define arch_lib_nss @"/lib/ppc"
#define arch_str "/ppc"
#else
#ifdef __ppc64__
#define arch_lib_nss @"/lib/ppc64"
#define arch_str "/ppc64"
#else
#ifdef __i386__
#define arch_lib_nss @"/lib/i386"
#define arch_str "/i386"
#else
#ifdef __x86_64__
#define arch_lib_nss @"/lib/x86_64"
#define arch_str "/x86_64"
#endif
#endif
#endif
#endif
#ifdef arch_lib_nss
	if (!getenv("R_ARCH")) {
		if ([[NSFileManager defaultManager] fileExistsAtPath:[[NSString stringWithUTF8String:getenv("R_HOME")] stringByAppendingString: arch_lib_nss]])
			setenv("R_ARCH", arch_str, 1);
	}
#else
#warning "Unknown architecture, R_ARCH won't be set automatically."
#endif
    
	return self;
}

- (BOOL) activate
{
	SLog(@"REngine.activate: starting R ...");
	RENGINE_BEGIN;
	{
		int res = initR(argc, argv, [saveAction isEqual:@"yes"]?Rinit_save_yes:([saveAction isEqual:@"no"]?Rinit_save_no:Rinit_save_ask));
		active = (res==0)?YES:NO;
	}
	RENGINE_END;
	if (lastInitRError) {
		if (lastError) [lastError release];
		lastError = [[NSString alloc] initWithUTF8String:lastInitRError];
	} else lastError=nil;
	SLog(@"REngine.activate: %@", (lastError)?lastError:@"R started with no error");
    return active;
}

- (NSString*) lastError
{
	return lastError;
}

- (BOOL) isActive { return active; }
- (BOOL) isLoopRunning { return loopRunning; }

- (BOOL) allowEvents { return (maskEvents==0); }

- (BOOL) beginProtected {
	SLog(@"REngine.beginProtected, maskEvents=%d, protectedMode=%d", maskEvents, (int)protectedMode);
	if (protectedMode) return NO;
	maskEvents++;
	protectedMode=YES;
	return YES;
}

- (void) endProtected {
	SLog(@"REngine.endProtected, maskEvents=%d, protectedMode=%d", maskEvents, (int)protectedMode);
	maskEvents--;
	protectedMode=NO;
}

- (void) run: (id) sender
{
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    replThread = [NSThread currentThread];
    [self runREPL];
    replThread = nil;
    [pool release];
}

- (void) runREPL
{
	BOOL keepInLoop = YES;
	if (!active) return;
	loopRunning=YES;
	while (keepInLoop) {
		insideR++;
		@try {
			run_REngineRmainloop(0);
			insideR--;
			keepInLoop = NO; // voluntary exit, break the loop
		}
		@catch (NSException *foo) {
			insideR--;
			NSLog(@"*** REngine.runREPL: caught ObjC exception in the main loop. Update to the latest GUI version and consider reporting this properly (see FAQ) if it persists and is not known. \n*** reason: %@\n*** name: %@, info: %@\n*** Version: R %s.%s (%s) R.app %@%s\nConsider saving your work soon in case this develops into a problem.", [foo reason], [foo name], [foo userInfo], R_MAJOR, R_MINOR, R_SVN_REVISION, [[NSBundle mainBundle] objectForInfoDictionaryKey:@"CFBundleShortVersionString"], getenv("R_ARCH"));
		}
	}
	loopRunning=NO;	
}

- (void) runDelayedREPL
{
	if (!active) return;
	loopRunning=YES;
	insideR++;
    run_REngineRmainloop(1);
	insideR--;
	/* in fact loopRunning is not determinable, because later longjmp may have re-started the loop, so we just keep it at YES */
}

- (id <REPLHandler>) handler
{
    return replHandler;
}

- (id <CocoaHandler>) cocoaHandler
{
	return cocoaHandler;
}

- (void) setCocoaHandler: (id <CocoaHandler>) ch
{
	cocoaHandler=ch;
}

- (void) setSaveAction: (NSString*) action
{
	saveAction = action?action:@"ask";
}

- (NSString*) saveAction
{
	return saveAction;
}

- (void) disableRSignalHandlers: (BOOL) disable
{
	setRSignalHandlers(disable?0:1);
}

- (void) begin
{
	// FIXME: we should set a lock here
	[replHandler handleBusy:YES];
	if (insideR) SLog(@"***********> REngine.begin: expected insideR to be 0, but it's %d", insideR);
	insideR++;
}

- (void) end
{
	// FIXME: we should release a lock here
	insideR--;
	if (insideR) SLog(@"***********> REngine.end: expected insideR to be 0, but it's %d", insideR);
	[replHandler handleBusy:NO];
}

- (RSEXP*) parse: (NSString*) str
{
    return [self parse: str withParts: 1];
}

- (RSEXP*) parse: (NSString*) str withParts: (int) count
{
    ParseStatus ps;
    SEXP pstr, cv;

	if (!active) return nil;
	RENGINE_BEGIN;
    PROTECT(cv=allocVector(STRSXP, 1));
    SET_STRING_ELT(cv, 0, mkChar([str UTF8String]));    
    pstr=RS_ParseVector(cv, count, &ps);
    UNPROTECT(1);
	RENGINE_END;
    //NSLog(@"parse status: %d, SEXP: %x, type: %d\n", ps, pstr, TYPEOF(pstr));
	return pstr?[[RSEXP alloc] initWithSEXP: pstr]:nil;
}

- (RSEXP*) evaluateExpressions: (RSEXP*) expr
{
    SEXP es=0;
    int er=0;
    int i=0,l;

    //NSLog(@"evaluateExpressions: %@", expr);
	if (!active) return nil;
	RENGINE_BEGIN;
    // if we have an entire expression list, evaluate its contents one-by-one and return only the last one
    if ([expr type]==EXPRSXP) {
        l=[expr length];
        while (i<l) {
            //NSLog(@"expression %d: %@", i, [expr elementAt: i]);
            es=R_tryEval([[expr elementAt:i] directSEXP], R_GlobalEnv, &er);
			//NSLog(@"Eval result: %d [es=%x]\n",er,es);
            i++;
        }
    } else
        es=R_tryEval([expr directSEXP], R_GlobalEnv, &er);
	RENGINE_END;
        
    return es?[[RSEXP alloc] initWithSEXP: es]:nil;
}

- (RSEXP*) evaluateString: (NSString*) str
{
    RSEXP *ps, *xr;
	SLog(@"REngine.evaluateString:\"%@\"", str);
	if (!active) return nil;
    ps=[self parse: str];
    if (ps==nil) return nil;
	if([ps type]==NILSXP) { [ps release]; return nil; }
	DO_RENG_EVAL_STATUS(str);
    xr=[self evaluateExpressions: ps];
	DONE_RENG_EVAL_STATUS();
	[ps release];
	SLog(@" - result: %@", xr);
	return xr;
}

- (RSEXP*) evaluateString: (NSString*) str withParts: (int) count
{
    RSEXP *ps, *xr;
	SLog(@"REngine.evaluateString:\"%@\" withParts:%d", str, count);
	if (!active) return nil;
    ps=[self parse: str withParts: count];
    if (ps==nil) return nil;
	if([ps type]==NILSXP) { [ps release]; return nil; }
	DO_RENG_EVAL_STATUS(str);
    xr=[self evaluateExpressions: ps];
	DONE_RENG_EVAL_STATUS();
	[ps release];
	SLog(@" - result: %@", xr);
	return xr;
}

- (BOOL) executeString: (NSString*) str
{
    RSEXP *ps, *xr;
	BOOL success=NO;
	SLog(@"REngine.executeString:\"%@\"", str);
	if (!active) return NO;
    ps=[self parse: str];
    if (ps==nil) return NO;
	DO_RENG_EVAL_STATUS(str);
    xr=[self evaluateExpressions: ps];
	DONE_RENG_EVAL_STATUS();
	[ps release];
	if (xr!=nil) success=YES;
	if (xr) [xr release];
	SLog(@" - success: %@", success?@"YES":@"NO");
	return success;
}

@end
