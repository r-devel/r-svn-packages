#ifndef R_DEV_GTK_H
#define R_DEV_GTK_H

#include <gtk/gtk.h>

#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#if R_VERSION < R_Version(2, 7, 0)
# include <Rgraphics.h>
# include <Rdevices.h>
# include <R_ext/GraphicsDevice.h>
typedef NewDevDesc* pDevDesc;
typedef GEDevDesc* pGEDevDesc;
typedef R_GE_gcontext* pGEcontext;
#endif
#include <R_ext/GraphicsEngine.h>

typedef struct {
    /* R Graphics Parameters */
    /* Local device copy so that we can detect */
    /* when parameter changes. */

    double cex;				/* Character expansion */
    double srt;				/* String rotation */

    /* gint bg; */                      /* Background */
    int fill;
    int col;

    gint lty, lwd;                      /* line params */

    /* GTK Driver Specific */

    int windowWidth;			/* Window width (pixels) */
    int windowHeight;			/* Window height (pixels) */
    gboolean resize;			/* Window resized */
    GtkWidget *window;			/* Graphics frame */
    GtkWidget *drawing;                 /* Drawable window */

    GdkPixmap *pixmap;                  /* Backing store */

    GdkGC *wgc;
    GdkColor gcol_bg;
    GdkRectangle clip;
    GdkCursor *gcursor;

    int fontface;			/* Typeface */
    int fontsize;			/* Size in points */
    gboolean usefixed;

    PangoFont *font;
    PangoFontDescription *fontdesc;
    PangoContext *context; /* NOT USED YET */

} gtkDesc;


Rboolean GTKDeviceDriver(pDevDesc dd, char *display, double width, 
			 double height, double pointsize);

Rboolean GTKDeviceFromWidget(pDevDesc dd, char *w, double width, 
			     double height, double pointsize);

#endif /* ifndef R_DEV_GTK_H */
