#include <gtk/gtk.h>
#include <R.h>
#ifndef WIN32
#include "R_ext/eventloop.h"
#include <gdk/gdkx.h>
#else
#include <sys/types.h>
#endif


void
R_gtk_eventHandler(void *userData)
{
    while (gtk_events_pending())
	gtk_main_iteration();  
}

void
R_gtk_setEventHandler()
{
    static InputHandler *h = NULL;
    if(!h)
	h = addInputHandler(R_InputHandlers, ConnectionNumber(GDK_DISPLAY()),
			    R_gtk_eventHandler, -1);
}
