/*
 *  R.app : a Cocoa front end to: "R A Computer Language for Statistical Data Analysis"
 *  
 *  R.app Copyright notes:
 *                     Copyright (C) 2004  The R Foundation
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
 */

#import <Cocoa/Cocoa.h>
#include "Rinit.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
#import "REngine.h"

static REngine* mainRengine=nil;

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
	
    //setenv("R_HOME","/Library/Frameworks/R.framework/Resources",1);
    //setenv("DYLD_LIBRARY_PATH","/Library/Frameworks/R.framework/Resources/lib",1);
    
	return self;
}

- (BOOL) activate
{
	RENGINE_BEGIN;
	{
		int res = initR(argc, argv);
		active = (res==0)?YES:NO;
	}
	RENGINE_END;
	if (lastInitRError) {
		if (lastError) [lastError release];
		lastError = [[NSString alloc] initWithCString:lastInitRError];
	} else lastError=nil;
    return active;
}

- (NSString*) lastError
{
	return lastError;
}

- (BOOL) isActive { return active; }
- (BOOL) isLoopRunning { return loopRunning; }

- (void) runREPL
{
	if (!active) return;
	loopRunning=YES;
    run_REngineRmainloop();
	loopRunning=NO;	
}

- (id) handler
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

- (void) begin
{
	// FIXME: we should set a lock here
	[replHandler handleBusy:YES];
}

- (void) end
{
	// FIXME: we should release a lock here
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
    SET_STRING_ELT(cv, 0, mkChar([str cString]));    
    pstr=R_ParseVector(cv, count, &ps);
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
	if (!active) return nil;
    ps=[self parse: str];
    if (ps==nil) return nil;
	if([ps type]==NILSXP) { [ps release]; return nil; }
    xr=[self evaluateExpressions: ps];
	[ps release];
	return xr;
}

- (RSEXP*) evaluateString: (NSString*) str withParts: (int) count
{
    RSEXP *ps, *xr;
	if (!active) return nil;
    ps=[self parse: str withParts: count];
    if (ps==nil) return nil;
	if([ps type]==NILSXP) { [ps release]; return nil; }
    xr=[self evaluateExpressions: ps];
	[ps release];
	return xr;
}

- (BOOL) executeString: (NSString*) str
{
    RSEXP *ps, *xr;
	BOOL success=NO;
	if (!active) return NO;
    ps=[self parse: str];
    if (ps==nil) return NO;
    xr=[self evaluateExpressions: ps];
	[ps release];
	if (xr!=nil) success=YES;
	if (xr) [xr release];
	return success;
}

@end
