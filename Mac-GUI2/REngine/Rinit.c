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
 *                     Copyright (C) 2002-2005   The R Foundation
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


/* PrivateHeaders */
#include "config.h"
#include "Defn.h"
#include "Print.h"
#include "IOStuff.h"

#include <Rversion.h>
#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#endif

#include <R.h>
#include <Rinternals.h>
#include "Rinit.h"
#include "Rcallbacks.h"
#include <R_ext/Parse.h>

#include <Rdevices.h>
#include <R_ext/GraphicsDevice.h>

#include <Rdefines.h>

/* This constant defines the maximal length of single ReadConsole input, which usually corresponds to the maximal length of a single line. The buffer is allocated dynamically, so an arbitrary size is fine. */
#ifndef MAX_R_LINE_SIZE
#define MAX_R_LINE_SIZE 32767
#endif

/*----------- more recent R versions have Rinterface ---------*/
#if (R_VERSION < R_Version(2,4,0))
#error "Your R version is too old, 2.4.0 or higher is required"
#endif

#define CSTACK_DEFNS
#include <Rinterface.h>

/* and SaveAction is not officially exported */
extern SA_TYPE SaveAction;

/* --- those are not part of the Rinterface --- */
extern SEXP (*ptr_do_packagemanger)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_datamanger)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_browsepkgs)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_wsbrowser)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_hsbrowser)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_dataentry)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_selectlist)(SEXP, SEXP, SEXP, SEXP);
extern void (*ptr_do_flushconsole)();
extern void (*ptr_R_ProcessEvents)();

extern int  (*ptr_R_EditFile)(const char *); /* in r-devel ptr_Raqua_Edit is no longer used*/
extern SEXP (*ptr_do_selectlist)(SEXP, SEXP, SEXP, SEXP);
extern int  (*ptr_Raqua_CustomPrint)(const char *, SEXP); /* custom print proxy for help/search/pkg-info */
extern int  (*ptr_CocoaSystem)(const char *);

int end_Rmainloop(void);    /* from src/main.c */
int Rf_initialize_R(int ac, char **av); /* from src/unix/system.c */

/*--- note: the REPL code was modified, R_ReplState is not the same as used internally in R */

typedef struct {
  ParseStatus    status;
  int            prompt_type;
  int            browselevel;
  int            buflen;
  unsigned char *buf;
  unsigned char *bufp;
} R_ReplState;

static int  RGUI_ReplIteration(SEXP rho, int savestack, int browselevel, R_ReplState *state);
static void RGUI_ReplConsole(SEXP rho, int savestack, int browselevel);

/*---------- implementation -----------*/

static R_ReplState state = {0, 1, 0, (MAX_R_LINE_SIZE+1), NULL, NULL};

char *lastInitRError = 0;

/* Note: R_SignalHandlers are evaluated in setup_Rmainloop which is called inside initR */
int initR(int argc, char **argv, int save_action) {
    if (!getenv("R_HOME")) {
        lastInitRError = "R_HOME is not set. Please set all required environment variables before running this program.";
        return -1;
    }
    
    int stat=Rf_initialize_R(argc, argv);
    if (stat<0) {
        lastInitRError = "Failed to initialize R!";;
        return -2;
    }

	if (state.buflen<128) state.buflen=1024;
	state.buf=(unsigned char*) malloc(state.buflen);
	
   // printf("R primary initialization done. Setting up parameters.\n");

    R_Outputfile = NULL;
    R_Consolefile = NULL;
    R_Interactive = 1;
    SaveAction = (save_action==Rinit_save_yes)?SA_SAVE:((save_action==Rinit_save_no)?SA_NOSAVE:SA_SAVEASK);

	R_CStackLimit = -1;

    /* ptr_R_Suicide = Re_Suicide; */
    /* ptr_R_CleanUp = Re_CleanUp; */
    ptr_R_ShowMessage = Re_ShowMessage;
    ptr_R_ReadConsole =  Re_ReadConsole;
    ptr_R_WriteConsole = Re_WriteConsole;
    ptr_R_ResetConsole = Re_ResetConsole;
    ptr_do_flushconsole = Re_FlushConsole;
    ptr_R_ClearerrConsole = Re_ClearerrConsole;
    ptr_R_Busy = Re_RBusy;
    ptr_R_ProcessEvents =  Re_ProcessEvents;

	ptr_R_EditFile = Re_Edit;
	ptr_Raqua_CustomPrint = Re_CustomPrint;
	
    ptr_R_ShowFiles = Re_ShowFiles;
	ptr_R_EditFiles = Re_EditFiles;
    ptr_R_ChooseFile = Re_ChooseFile;
	
	
	ptr_do_packagemanger = Re_packagemanger;
	ptr_do_datamanger = Re_datamanger;
	ptr_do_browsepkgs = Re_browsepkgs;
	ptr_do_wsbrowser = Re_do_wsbrowser;
	ptr_do_hsbrowser = Re_do_hsbrowser;
	
	ptr_CocoaSystem = Re_system;
	ptr_do_dataentry = Re_dataentry;
	ptr_do_selectlist = Re_do_selectlist;

	setup_Rmainloop();

    return 0;
}

static int firstRun=1;

void setRSignalHandlers(int val) {
	R_SignalHandlers = val;
}

void run_REngineRmainloop(int delayed)
{
    /* Here is the real R read-eval-loop. */
    /* We handle the console until end-of-file. */

	firstRun=delayed;
	
    R_IoBufferInit(&R_ConsoleIob);
    SETJMP(R_Toplevel.cjmpbuf);
    R_GlobalContext = R_ToplevelContext = &R_Toplevel;
#ifdef REINSTALL_SIGNAL_HANDLERS
    signal(SIGINT, handleInterrupt);
    signal(SIGUSR1,onsigusr1);
    signal(SIGUSR2,onsigusr2);
#ifdef Unix
    signal(SIGPIPE, onpipe);
#endif
#endif
	if (firstRun) {
		firstRun=0;
		return;
	}
	
    RGUI_ReplConsole(R_GlobalEnv, 0, 0);
	end_Rmainloop(); /* must go here */
}

extern void Re_WritePrompt(char *prompt);
extern char *R_PromptString(int browselevel, int type); /* from main.c */
extern void Rf_callToplevelHandlers(SEXP expr, SEXP value, Rboolean succeeded, 
			Rboolean visible);  /* from main.c */

#if (R_VERSION < R_Version(2,3,0))
extern int Rf_ParseBrowser(SEXP, SEXP);
#else
									/* FIXME: Rf_ParseBrowser is no longer exported
									we will need to update it manually - ugly!
									This is copy/paste from R 2.3.0 release
									*/
static void printwhere(void)
{
	RCNTXT *cptr;
	int lct = 1;
	
	for (cptr = R_GlobalContext; cptr; cptr = cptr->nextcontext) {
		if ((cptr->callflag & (CTXT_FUNCTION | CTXT_BUILTIN)) &&
			(TYPEOF(cptr->call) == LANGSXP)) {
				Rprintf("where %d: ", lct++);
				PrintValue(cptr->call);
		}
	}
	Rprintf("\n");
}
									
static int Rf_ParseBrowser(SEXP CExpr, SEXP rho)
{
	int rval = 0;
	if (isSymbol(CExpr)) {
		if (!strcmp(CHAR(PRINTNAME(CExpr)), "n")) {
			SET_DEBUG(rho, 1);
			rval = 1;
		}
		if (!strcmp(CHAR(PRINTNAME(CExpr)), "c")) {
			rval = 1;
			SET_DEBUG(rho, 0);
		}
		if (!strcmp(CHAR(PRINTNAME(CExpr)), "cont")) {
			rval = 1;
			SET_DEBUG(rho, 0);
		}
		if (!strcmp(CHAR(PRINTNAME(CExpr)), "Q")) {
												
			/* Run onexit/cend code for everything above the target.
			   The browser context is still on the stack, so any error
			   will drop us back to the current browser.  Not clear
			   this is a good thing.  Also not clear this should still
			   be here now that jump_to_toplevel is used for the
			   jump. */
			R_run_onexits(R_ToplevelContext);
												
			/* this is really dynamic state that should be managed as such */
			R_BrowseLevel = 0;
			SET_DEBUG(rho, 0); /*PR#1721*/
							
			jump_to_toplevel();
		}
		if (!strcmp(CHAR(PRINTNAME(CExpr)),"where")) {
			printwhere();
			/* SET_DEBUG(rho, 1); */
			rval = 2;
		}
	}
	return rval;
}
#endif


static void RGUI_ReplConsole(SEXP rho, int savestack, int browselevel)
{
    int status;

    R_IoBufferWriteReset(&R_ConsoleIob);
    state.buf[0] = '\0';
    state.buf[state.buflen-1] = '\0'; /* stopgap measure */
    state.bufp = state.buf;
    if(R_Verbose)
		REprintf(" >R_ReplConsole(): before \"for(;;)\" {Rinit.c}\n");
    for(;;) {
		status = RGUI_ReplIteration(rho, savestack, browselevel, &state);
		if(status < 0)
			return;
    }
}

static int cleanExit=1;

static int RGUI_ReplIteration(SEXP rho, int savestack, int browselevel, R_ReplState *state)
{
    int c, browsevalue;
    SEXP value;
    Rboolean wasDisplayed = FALSE;

	if (!cleanExit) // if previous eval was aborted, we may need to clear the busy flag
		R_Busy(0);
	
    if(!*state->bufp) {
	    if (Re_ReadConsole(R_PromptString(browselevel, state->prompt_type),
						   state->buf, state->buflen-1, 1) == 0) return(-1);
	    state->bufp = state->buf;
    }
#ifdef SHELL_ESCAPE
    if (*state->bufp == '!') {
#ifdef HAVE_SYSTEM
	    R_system(&(state->buf[1]));
#else
	    Rprintf("error: system commands are not supported in this version of R.\n");
#endif /* HAVE_SYSTEM */
	    state->buf[0] = '\0';
	    return(0);
    }
#endif /* SHELL_ESCAPE */
    while((c = *state->bufp++)) {
	    R_IoBufferPutc(c, &R_ConsoleIob);
	    if(c == ';' || c == '\n') break;
    }

    R_PPStackTop = savestack;
    R_CurrentExpr = R_Parse1Buffer(&R_ConsoleIob, 0, &state->status);
    switch(state->status) {
    case PARSE_NULL:

	    if (browselevel)
		    return(-1);
	    R_IoBufferWriteReset(&R_ConsoleIob);
	    state->prompt_type = 1;
	    return(1);

    case PARSE_OK:

 	    R_IoBufferReadReset(&R_ConsoleIob);
	    R_CurrentExpr = R_Parse1Buffer(&R_ConsoleIob, 1, &state->status);
	    if (browselevel) {
		    browsevalue = Rf_ParseBrowser(R_CurrentExpr, rho);
		    if(browsevalue == 1 )
			    return(-1);
		    if(browsevalue == 2 ) {
			    R_IoBufferWriteReset(&R_ConsoleIob);
			    return(0);
		    }
	    }
	    R_Visible = 0;
	    R_EvalDepth = 0;
	    PROTECT(R_CurrentExpr);
		cleanExit=0;
	    R_Busy(1);
	    value = eval(R_CurrentExpr, rho);
	    SET_SYMVALUE(R_LastvalueSymbol, value);
	    wasDisplayed = R_Visible;
	    if (R_Visible)
		    PrintValueEnv(value, rho);
	    if (R_CollectWarnings)
			PrintWarnings();
	    Rf_callToplevelHandlers(R_CurrentExpr, value, TRUE, wasDisplayed);
	    R_CurrentExpr = value; /* Necessary? Doubt it. */
	    R_Busy(0);
		cleanExit=1;
	    UNPROTECT(1);
	    R_IoBufferWriteReset(&R_ConsoleIob);
	    state->prompt_type = 1;
	    return(1);

    case PARSE_ERROR:

	    state->prompt_type = 1;
	    error("syntax error");
	    R_IoBufferWriteReset(&R_ConsoleIob);
	    return(1);

    case PARSE_INCOMPLETE:

	    R_IoBufferReadReset(&R_ConsoleIob);
	    state->prompt_type = 2;
	    return(2);

    case PARSE_EOF:

	    return(-1);
	    break;
    }

    return(0);
}

