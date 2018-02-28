#include <R.h>
#include <Rinternals.h>
#include <Rgraphics.h>
#include <Rdevices.h>
#include <R_ext/GraphicsDevice.h>
#include <R_ext/GraphicsEngine.h>

#define BEGIN_SUSPEND_INTERRUPTS
#define END_SUSPEND_INTERRUPTS

/* #include <gtk/gtk.h>a */
#include "devGNOME.h"

static void createGnomeDevice(char *display, double width, double height,
			      double ps, int aa)
{
    NewDevDesc *dev;
    GEDevDesc *dd;
    
    R_CheckDeviceAvailable();
    BEGIN_SUSPEND_INTERRUPTS {
        /* Allocate and initialize the device driver data */
        if (!(dev = (NewDevDesc *) calloc(1, sizeof(NewDevDesc))))
            return;
        /* Do this for early redraw attempts */
        dev->displayList = R_NilValue;
        /* Make sure that this is initialised before a GC can occur.
         * This (and displayList) get protected during GC
         */
        dev->savedSnapshot = R_NilValue;
        if (!GnomeDeviceDriver((DevDesc*)dev, display, width, height, ps)){
            free(dev);
            error("Unable to start GNOME device");
        }
        gsetVar(install(".Device"), mkString("gnome"), R_NilValue);
        dd = GEcreateDevDesc(dev);
        dd->newDevStruct = 1;	
       addDevice((DevDesc*) dd);
        GEinitDisplayList(dd);
    } END_SUSPEND_INTERRUPTS;

    gdk_flush();
}

void do_gnome(char **dpy, double *width, double *height, double *ps, int *aa)
{
    char *display, *vmax;
    vmax = vmaxget();

    display = R_alloc(strlen(dpy[0]) + 1, sizeof(char));
    strcpy (display, dpy[0]);
    if (*width <= 0 || *height <= 0)
        error("invalid width or height for GNOME device");

    createGnomeDevice(display, *width, *height, *ps, 0); 
}
