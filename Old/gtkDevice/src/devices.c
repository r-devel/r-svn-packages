#include <gtk/gtk.h>
#include "devGTK.h"

#ifndef BEGIN_SUSPEND_INTERRUPTS
# define BEGIN_SUSPEND_INTERRUPTS
# define END_SUSPEND_INTERRUPTS
#endif

typedef Rboolean 
(*GtkDeviceCreateFun)(pDevDesc, char *display, double width, 
		      double height, double pointsize);

static pGEDevDesc
createGtkDevice(char *display, double width, double height, double ps, 
		GtkDeviceCreateFun init_fun)
{
    pGEDevDesc dd;
    pDevDesc dev;

    R_GE_checkVersionOrDie(R_GE_version);
    R_CheckDeviceAvailable();
    BEGIN_SUSPEND_INTERRUPTS {
	/* Allocate and initialize the device driver data */
	if (!(dev = (pDevDesc) calloc(1, sizeof(NewDevDesc))))
	    return NULL;
#if R_VERSION < R_Version(2, 7, 0)
	/* Do this for early redraw attempts */
	dev->displayList = R_NilValue;
#endif
	if (! init_fun (dev, display, width, height, ps)) {
	    free(dev);
	    PROBLEM  "unable to start device gtk" ERROR;
	}
#if R_VERSION < R_Version(2, 7, 0)
	gsetVar(install(".Device"), mkString("GTK"), R_NilValue);
	dd = GEcreateDevDesc(dev);
	Rf_addDevice((DevDesc*) dd);
	GEinitDisplayList(dd);
#else
	dd = GEcreateDevDesc(dev);
	GEaddDevice2(dd, "GTK");
#endif
    } END_SUSPEND_INTERRUPTS;

    gdk_flush();

    return(dd);
}


void
do_GTK(char **dpy, double *in_width, double *in_height, double *in_pointsize)
{
    char *display, *vmax;
    double height, width, ps;

/*    gcall = call; */
    vmax = vmaxget();
    display = R_alloc(strlen(dpy[0]) + 1, sizeof(char));
    strcpy(display, dpy[0]);
    width = *in_width;
    height = *in_height;
    if (width <= 0 || height <= 0) {
	PROBLEM "Gtk device: invalid width or height" ERROR;
    }
    ps = *in_pointsize;
 
    createGtkDevice(display, width, height, ps, GTKDeviceDriver);

    vmaxset(vmax);
    /*   return R_NilValue; */
}


SEXP
do_asGTKDevice(SEXP widget, SEXP w, SEXP h, SEXP pointsize)
{
    GtkWidget *drawing_widget = (GtkWidget*) R_ExternalPtrAddr(widget);
    double width, height, ps;
    SEXP ans = Rf_allocVector(LGLSXP, 1);

    width = asReal(w);
    height = asReal(h);
    ps = asReal(pointsize);

    LOGICAL(ans)[0] = (createGtkDevice((char *) drawing_widget, width, height, ps, GTKDeviceFromWidget) != NULL);

    return(ans);
}

