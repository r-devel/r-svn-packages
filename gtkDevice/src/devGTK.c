/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2001   Lyndon Drake
 *                            and the R Development Core Team
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef GTK2

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_RINT
#define R_rint(x) rint(x)
#else
#define R_rint(x) ((int) x + 0.5)
#endif

#include <gtk/gtk.h>
#include <R.h>
#include <Rinternals.h>
#include <Rgraphics.h>
#include <Rdevices.h>
#include <R_ext/GraphicsDevice.h>
#include <R_ext/GraphicsEngine.h>

#include "devGTK.h"
#include "gdkrotated.h"

/*
#define DEBUG_GTK
*/

#define CURSOR		GDK_CROSSHAIR		/* Default cursor */
#define MM_PER_INCH	25.4			/* mm -> inch conversion */

#define IS_100DPI ((int) (1./pixelHeight() + 0.5) == 100)

static unsigned int adobe_sizes = 0x0403165D;
#define ADOBE_SIZE(I) ((I) > 7 && (I) < 35 && (adobe_sizes & (1<<((I)-8))))

/* routines from here */

/* Device driver actions */
static void GTK_Activate(NewDevDesc *dd);
static void GTK_Circle(double x, double y, double r,
		       int col, int fill, double gamma, int lty, double lwd,
		       NewDevDesc *dd);
static void GTK_Clip(double x0, double x1, double y0, double y1, 
		     NewDevDesc *dd);
static void GTK_Close(NewDevDesc *dd);
static void GTK_Deactivate(NewDevDesc *dd);
static void GTK_Hold(NewDevDesc *dd);
static Rboolean GTK_Locator(double *x, double *y, NewDevDesc *dd);
static void GTK_Line(double x1, double y1, double x2, double y2,
		     int col, double gamma, int lty, double lwd,
		     NewDevDesc *dd);
static void GTK_MetricInfo(int c, int font, double cex, double ps,
			      double* ascent, double* descent,
			      double* width, NewDevDesc *dd);
static void GTK_Mode(int mode, NewDevDesc *dd);
static void GTK_NewPage(int fill, double gamma, NewDevDesc *dd);
static void GTK_Polygon(int n, double *x, double *y, 
			int col, int fill, double gamma, int lty, double lwd,
			NewDevDesc *dd);
static void GTK_Polyline(int n, double *x, double *y, 
			    int col, double gamma, int lty, double lwd,
			    NewDevDesc *dd);
static void GTK_Rect(double x0, double y0, double x1, double y1,
		     int col, int fill, double gamma, int lty, double lwd,
		     NewDevDesc *dd);
static void GTK_Size(double *left, double *right,
		     double *bottom, double *top,
		     NewDevDesc *dd);
static double GTK_StrWidth(char *str, int font,
			      double cex, double ps, NewDevDesc *dd);
static void GTK_Text(double x, double y, char *str, 
		     double rot, double hadj, 
		     int col, double gamma, int font, double cex, double ps,
		     NewDevDesc *dd);
static Rboolean GTK_Open(NewDevDesc*, gtkDesc*, char*, double, double);

static gint initialize(NewDevDesc *dd);

/* Pixel Dimensions (Inches) */

static double pixelWidth(void)
{
    double width, widthMM;
    width = gdk_screen_width();
    widthMM = gdk_screen_width_mm();
    return ((double)widthMM / (double)width) / MM_PER_INCH;
}

static double pixelHeight(void)
{
    double height, heightMM;
    height = gdk_screen_height();
    heightMM = gdk_screen_height_mm();
    return ((double)heightMM / (double)height) / MM_PER_INCH;
}

/* font stuff */

static char *fontname_R6 = "-*-helvetica-%s-%s-*-*-%d-*-*-*-*-*-*-*";
static char *symbolname	 = "-adobe-symbol-*-*-*-*-%d-*-*-*-*-*-*-*";

static char *slant[] = {"r", "o"};
static char *weight[] = {"medium", "bold"};

static GHashTable *font_htab = NULL;

struct _FontMetricCache {
    gint ascent[255];
    gint descent[255];
    gint width[255];
    gint font_ascent;
    gint font_descent;
    gint max_width;
};

#define SMALLEST 2

static GdkFont* RGTKLoadFont (gint face, gint size)
{
    gchar *fontname;
    GdkFont *tmp_font;
    gint pixelsize;

    if (face < 1 || face > 5)
	face = 1;

    if (size < SMALLEST) size = SMALLEST;

    /* Here's a 1st class fudge: make sure that the Adobe design sizes
       8, 10, 11, 12, 14, 17, 18, 20, 24, 25, 34 can be obtained via
       an integer "size" at 100 dpi, namely 6, 7, 8, 9, 10, 12, 13,
       14, 17, 18, 24 points. It's almost y = x * 100/72, but not
       quite. The constants were found using lm(). --pd */
    if (IS_100DPI) size = R_rint (size * 1.43 - 0.4);

    /* 'size' is the requested size, 'pixelsize' the size of the
       actually allocated font*/
    pixelsize = size;

    if(face == 5)
	fontname = g_strdup_printf(symbolname, pixelsize);
    else
	fontname = g_strdup_printf(fontname_R6,
				   weight[(face-1)%2],
				   slant[((face-1)/2)%2],
				   pixelsize);
#ifdef DEBUG_GTK
    Rprintf("loading:\n%s\n", fontname);
#endif
    tmp_font = gdk_font_load(fontname);
#ifdef DEBUG_GTK
    if (tmp_font) Rprintf("success\n"); else Rprintf("failure\n");
#endif

    if (!tmp_font) {
	static int near[]=
	{14,14,14,17,17,18,20,20,20,20,24,24,24,25,25,25,25};
	/* 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29  */
	
	/* If ADOBE_SIZE(pixelsize) is true at this point then
	   the user's system does not have the standard ADOBE font set
	   so we just have to use a "fixed" font.
	   If we can't find a "fixed" font then something is seriously
	   wrong */
	if (ADOBE_SIZE (pixelsize)) {
	    tmp_font = gdk_font_load ("fixed");
	    if (!tmp_font)
		error("Could not find any X11 fonts\nCheck that the Font Path is correct.");
	}
	
	if (pixelsize < 8)
	    pixelsize = 8;
	else if (pixelsize == 9)
	    pixelsize = 8;
	else if (pixelsize >= 13 && pixelsize < 30) 
	    pixelsize = near[size-13];
	else
	    pixelsize = 34;
	
	g_free(fontname);
	if(face == 5)
	    fontname = g_strdup_printf(symbolname, pixelsize);
	else
	    fontname = g_strdup_printf(fontname_R6,
				       weight[(face-1)%2],
				       slant[((face-1)/2)%2],
				       pixelsize);	    
#ifdef DEBUG_GTK
	Rprintf("loading:\n%s\n", fontname);
#endif
	tmp_font = gdk_font_load (fontname);
#ifdef DEBUG_GTK
	if (tmp_font) Rprintf("success\n"); else Rprintf("failure\n");
#endif
    }
    
    if(tmp_font) {
	g_hash_table_insert(font_htab, (gpointer) g_strdup(fontname),
			    (gpointer) tmp_font);
	if (fabs( (pixelsize - size)/(double)size ) > 0.1)
	    warning("GTK used font size %d when %d was requested",
		    pixelsize, size);
    }

    g_free(fontname);
    return tmp_font;
}

static void SetFont(NewDevDesc *dd, gint face, gint size)
{
    GdkFont *tmp_font;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    if (face < 1 || face > 5) face = 1;

    if (!gtkd->usefixed &&
	(size != gtkd->fontsize	|| face != gtkd->fontface)) {
	tmp_font = RGTKLoadFont(face, size);
	if(tmp_font) {
	    gtkd->font = tmp_font;
	    gtkd->fontface = face;
	    gtkd->fontsize = size;
	    gdk_gc_set_font(gtkd->wgc, tmp_font);
	} else 
	    error("X11 font at size %d could not be loaded", size);
    }
}

static gint SetBaseFont(gtkDesc *gtkd)
{
    gtkd->fontface = 1;
    gtkd->fontsize = 12;
    gtkd->usefixed = 0;

    if(font_htab == NULL) {
	font_htab = g_hash_table_new (g_str_hash, g_str_equal);
    }

    gtkd->font = RGTKLoadFont (gtkd->fontface, gtkd->fontsize);

    if(gtkd->font != NULL) {
	gdk_gc_set_font(gtkd->wgc, gtkd->font);
	return 1;
    }

    gtkd->usefixed = 1;
    gtkd->font = gdk_font_load ("fixed");

    if(gtkd->font != NULL) {
	gdk_gc_set_font(gtkd->wgc, gtkd->font);
	return 1;
    }

    return 0;
}


/* set the r, g, b, and pixel values of gcol to color */
static void SetColor(GdkColor *gcol, int color)
{
    int red, green, blue;

    red = R_RED(color);
    green = R_GREEN(color);
    blue = R_BLUE(color);
    gcol->red = 0;
    gcol->green = 0;
    gcol->blue = 0;
    gcol->pixel = gdk_rgb_xpixel_from_rgb((red << 16)|(green << 8)|(blue));
}

/* set the line type */
static void SetLineType(NewDevDesc *dd, int newlty, double nlwd)
{
    static gint8 dashlist[8];
    gint i, j, newlwd;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    newlwd = nlwd;
    if(newlty != gtkd->lty || newlwd != gtkd->lwd) {
	gtkd->lty = newlty;
	gtkd->lwd = newlwd;

	if(newlty == 0) {
	    if(newlwd <= 1)
		newlwd = 0;

	    gdk_gc_set_line_attributes(gtkd->wgc, newlwd,
				       GDK_LINE_SOLID,
				       GDK_CAP_BUTT,
				       GDK_JOIN_ROUND);
	}
	else {
	    if(newlwd < 1)
		newlwd = 1;

	    for(i = 0; (i < 8) && (newlty != 0); i++) {
		j = newlty & 15;
		if(j == 0) j = 1;
		j = j * newlwd;
		if(j > 255) j = 255;
		dashlist[i] = j;
		newlty = newlty >> 4;
	    }

	    /* set dashes */
	    gdk_gc_set_dashes(gtkd->wgc, 0, dashlist, i);
	    gdk_gc_set_line_attributes(gtkd->wgc, newlwd,
				       GDK_LINE_ON_OFF_DASH,
				       GDK_CAP_BUTT,
				       GDK_JOIN_ROUND);
	}
    }
}

/* signal functions */

static gint realize_event(GtkWidget *widget, gpointer data)
{
    NewDevDesc *dd;

    dd = (NewDevDesc *) data;
    g_return_val_if_fail(dd != NULL, FALSE);

    return(initialize(dd));
}


static gint initialize(NewDevDesc *dd)
{
    gtkDesc *gtkd;
    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    /* create gc */
    gtkd->wgc = gdk_gc_new(gtkd->drawing->window);

    /* set the cursor */
    gtkd->gcursor = gdk_cursor_new(GDK_CROSSHAIR);
    gdk_window_set_cursor(gtkd->drawing->window, gtkd->gcursor);

    /* set window bg */
    gdk_window_set_background(gtkd->drawing->window, &gtkd->gcol_bg);

    if(gtkd->wgc)
	gdk_gc_set_foreground(gtkd->wgc, &gtkd->gcol_bg);

    /* create offscreen drawable */
    if(gtkd->windowWidth > 0 && gtkd->windowHeight > 0) {
	gtkd->pixmap = gdk_pixmap_new(gtkd->drawing->window,
				      gtkd->windowWidth, gtkd->windowHeight,
				      -1);
	gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
			   gtkd->windowWidth, gtkd->windowHeight);
    }

    return FALSE;
}

static gint configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{
    NewDevDesc *dd;
    gtkDesc *gtkd;

    dd = (NewDevDesc *) data;
    g_return_val_if_fail(dd != NULL, FALSE);

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    /* check for resize */
    if((GTK_WIDGET_REALIZED(gtkd->drawing)) && ((gtkd->windowWidth != event->width) || (gtkd->windowHeight != event->height))) {
	gtkd->windowWidth = event->width;
	gtkd->windowHeight = event->height;

	gtkd->resize = TRUE;
    }

    return FALSE;
}

static void GTK_resize(NewDevDesc *dd);

static gint expose_event(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
    NewDevDesc *dd;
    gtkDesc *gtkd;

    dd = (NewDevDesc *) data;
    g_return_val_if_fail(dd != NULL, FALSE);

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    if(gtkd->wgc == NULL)
	realize_event(gtkd->drawing, dd);
    else if(gtkd->resize != 0) {
	GTK_resize(dd); 
    }
    else if (gtkd->pixmap) {
	gdk_draw_pixmap(gtkd->drawing->window, gtkd->wgc, gtkd->pixmap,
			event->area.x, event->area.y, event->area.x, event->area.y,
			event->area.width, event->area.height);
    }
    else {
	GEplayDisplayList((GEDevDesc*) Rf_GetDevice(Rf_devNumber((DevDesc*)dd)));
    }

    return FALSE;
}

static gint delete_event(GtkWidget *widget, GdkEvent *event, gpointer data)
{
    NewDevDesc *dd;

    dd = (NewDevDesc *) data;
    g_return_val_if_fail(dd != NULL, FALSE);

    Rf_KillDevice((DevDesc*) Rf_GetDevice(Rf_devNumber((DevDesc*) dd)));

    return TRUE;
}

/*
static void tb_activate_cb(GtkWidget *widget, gpointer data)
{
    NewDevDesc *dd;

    dd = (NewDevDesc *) data;
    g_return_if_fail(dd != NULL);

    selectDevice(Rf_devNumber((DevDesc*)dd));
}

static void tb_close_cb(GtkWidget *widget, gpointer data)
{
    NewDevDesc *dd;

    dd = (NewDevDesc *) data;
    g_return_if_fail(dd != NULL);

    Rf_KillDevice((DevDesc*) Rf_GetDevice(Rf_devNumber((DevDesc*) dd)));
}
*/

/* create window etc */
static Rboolean GTK_Open(NewDevDesc *dd, gtkDesc *gtkd, char *dsp, double w, double h)
{
    gint iw, ih;

    /* initialise pointers */
    gtkd->drawing = NULL;
    gtkd->wgc = NULL;
    gtkd->gcursor = NULL;

    /* initialise colour */
    gdk_rgb_init();
    gtk_widget_push_visual(gdk_rgb_get_visual());
    gtk_widget_push_colormap(gdk_rgb_get_cmap());

    /* create window etc */
    gtkd->windowWidth = iw = w / pixelWidth();
    gtkd->windowHeight = ih = h / pixelHeight();

    gtkd->window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    gtk_window_set_policy(GTK_WINDOW(gtkd->window), TRUE, TRUE, FALSE);
    gtk_widget_realize(gtkd->window);

    /* create drawingarea */
    gtkd->drawing = gtk_drawing_area_new();
    gtk_widget_set_events(gtkd->drawing,
			  GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK);

    /* connect to signal handlers, etc */
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "realize",
		       (GtkSignalFunc) realize_event, (gpointer) dd);

    /* drawingarea properties */
    gtk_widget_set_usize(gtkd->drawing, iw, ih);

    /* setup background color */
    /* gtkd->bg = dd->bg = R_RGB(255, 255, 255); */
    SetColor(&gtkd->gcol_bg, R_RGB(255, 255, 255)); /* FIXME canvas color */

    /* place and realize the drawing area */
    gtk_container_add(GTK_CONTAINER(gtkd->window), gtkd->drawing);
    gtk_widget_realize(gtkd->drawing);

    /* connect to signal handlers, etc */
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "configure_event",
		       (GtkSignalFunc) configure_event, (gpointer) dd);
    gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "expose_event",
		       (GtkSignalFunc) expose_event, (gpointer) dd);
    gtk_signal_connect(GTK_OBJECT(gtkd->window), "delete_event",
		       (GtkSignalFunc) delete_event, (gpointer) dd);

    /* show everything */
    gtk_widget_show_all(gtkd->window);
  
    /* initialise line params */
    gtkd->lty = -1;
    gtkd->lwd = -1;

    /* create offscreen drawable */
    gtkd->pixmap = gdk_pixmap_new(gtkd->drawing->window,
				  gtkd->windowWidth, gtkd->windowHeight,
				  -1);
    gdk_gc_set_foreground(gtkd->wgc, &gtkd->gcol_bg);
    gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
		       gtkd->windowWidth, gtkd->windowHeight);


    /* let other widgets use the default colour settings */
    gtk_widget_pop_visual();
    gtk_widget_pop_colormap();

    /* Set base font */
    if(!SetBaseFont(gtkd)) {
	Rprintf("can't find X11 font\n");
	return FALSE;
    }

    /* we made it! */
    return TRUE;
}

static double GTK_StrWidth(char *str, int font,
			      double cex, double ps, NewDevDesc *dd)
{
    int size;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    size = cex * ps + 0.5;
    SetFont(dd, font, size);

    return (double) gdk_string_width(gtkd->font, str);
}

static void GTK_MetricInfo(int c, int font, double cex, double ps,
			      double* ascent, double* descent,
			      double* width, NewDevDesc *dd)
{
    gint size;
    gint lbearing, rbearing, iascent, idescent, iwidth;
    gint maxwidth;
    gchar tmp[2];
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    size = cex * ps + 0.5;
    SetFont(dd, font, size);

    if(c == 0) {
	maxwidth = 0;

	for(c = 0; c <= 255; c++) {
	    g_snprintf(tmp, 2, "%c", (gchar) c);
	    iwidth = gdk_string_width(gtkd->font, tmp);
	    if (iwidth > maxwidth)
		maxwidth = iwidth;
	}

	*ascent = (double) gtkd->font->ascent;
	*descent = (double) gtkd->font->descent;
	*width = (double) maxwidth;
    }
    else {
	g_snprintf(tmp, 2, "%c", (gchar) c);
	gdk_string_extents(gtkd->font, tmp,
			   &lbearing, &rbearing,
			   &iwidth, &iascent, &idescent);

	*ascent = (double) iascent;
	*descent = (double) idescent;
	*width = (double) iwidth;
	/* This formula does NOT work for spaces
	*width = (double) (lbearing+rbearing);
	*/
    }
}

/* set clipping */
static void GTK_Clip(double x0, double x1, double y0, double y1, NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    gtkd->clip.x = dd->clipLeft = (int) MIN(x0, x1);
    gtkd->clip.width = abs( (int) x0 - (int) x1) + 1;
    dd->clipRight = dd->clipLeft + gtkd->clip.width;

    gtkd->clip.y = dd->clipBottom = (int) MIN(y0, y1);
    gtkd->clip.height = abs( (int) y0 - (int) y1) + 1;
    dd->clipTop = dd->clipBottom + gtkd->clip.height;

    /* Setting the clipping rectangle works when drawing to a window
       but not to the backing pixmap. This is a GTK+ bug that is
       unlikely to be fixed in this version (9 Jul 2002) - MTP
    */
    /* gdk_gc_set_clip_rectangle(gtkd->wgc, &gtkd->clip); */
}

static void GTK_Size(double *left, double *right,
		     double *bottom, double *top,
		     NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    *left = 0.0;
    *right =  gtkd->windowWidth;
    *bottom = gtkd->windowHeight;
    *top = 0.0;
}

static void GTK_resize(NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    if (gtkd->resize != 0) {
	dd->left = 0.0;
	dd->right = gtkd->windowWidth;
	dd->bottom = gtkd->windowHeight;
	dd->top = 0.0;
	gtkd->resize = 0;

	if(gtkd->pixmap)
	    gdk_pixmap_unref(gtkd->pixmap);
	if(gtkd->windowWidth > 0 && gtkd->windowHeight > 0) {
	    gtkd->pixmap = gdk_pixmap_new(gtkd->drawing->window,
					  gtkd->windowWidth, gtkd->windowHeight,
					  -1);
	    if(gtkd->wgc) {
		gdk_gc_set_foreground(gtkd->wgc, &gtkd->gcol_bg);
		gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
				   gtkd->windowWidth, gtkd->windowHeight);
	    }
	}
	GEplayDisplayList((GEDevDesc*) Rf_GetDevice(Rf_devNumber((DevDesc*)dd)));
    }
}

/* clear the drawing area */
static void GTK_NewPage(int fill, double gamma, NewDevDesc *dd)
{
    gtkDesc *gtkd;

    g_return_if_fail(dd != NULL);

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);
    g_return_if_fail(gtkd->drawing != NULL);
    g_return_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing));

    if(gtkd->drawing->window == NULL)
	return;
    if(gtkd->fill != fill && R_OPAQUE(fill)) {
	SetColor(&gtkd->gcol_bg, fill);
	gtkd->fill = fill;
	gdk_window_set_background(gtkd->drawing->window, &gtkd->gcol_bg);
    }

    gdk_window_clear(gtkd->drawing->window);

    if(gtkd->wgc) {
	gdk_gc_set_foreground(gtkd->wgc, &gtkd->gcol_bg);
	gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
			   gtkd->windowWidth, gtkd->windowHeight);
    }
}

/* kill off the window etc */
static void GTK_Close(NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    if(gtkd->window)
	gtk_widget_destroy(gtkd->window);

    if(gtkd->pixmap)
	gdk_pixmap_unref(gtkd->pixmap);

    gdk_flush();
    free(gtkd);
}

#define title_text_inactive "R graphics device %d"
#define title_text_active "R graphics device %d - Active"

static void GTK_Activate(NewDevDesc *dd)
{
    gtkDesc *gtkd;
    gint devnum;
    gchar *title_text;

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);
    if(!gtkd->window)
	return;

    devnum = Rf_devNumber((DevDesc*)dd) + 1;

    title_text = g_strdup_printf(title_text_active, devnum);

    gtk_window_set_title(GTK_WINDOW(gtkd->window), title_text);

    g_free(title_text);
}

static void GTK_Deactivate(NewDevDesc *dd)
{
    gtkDesc *gtkd;
    gint devnum;
    gchar *title_text;

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_if_fail(gtkd != NULL);
    if(!gtkd->window)
	return;

    devnum = Rf_devNumber((DevDesc*)dd) + 1;

    title_text = g_strdup_printf(title_text_inactive, devnum);

    gtk_window_set_title(GTK_WINDOW(gtkd->window), title_text);

    g_free(title_text);
}

/* drawing stuff */

static void GTK_Rect(double x0, double y0, double x1, double y1,
		     int col, int fill, double gamma, int lty, double lwd,
		     NewDevDesc *dd)
{
    double tmp;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill, gcol_outline;


    if(!gtkd->drawing->window)
	return;

    if(x0 > x1) {
	tmp = x0;
	x0 = x1;
	x1 = tmp;
    }
    if(y0 > y1) {
	tmp = y0;
	y0 = y1;
	y1 = tmp;
    }


    if (R_OPAQUE(fill)) {
	SetColor(&gcol_fill, fill);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	SetLineType(dd, lty, lwd);

	gdk_draw_rectangle(gtkd->drawing->window,
			   gtkd->wgc, TRUE,
			   (gint) x0, (gint) y0,
			   (gint) x1 - (gint) x0,
			   (gint) y1 - (gint) y0);
	gdk_draw_rectangle(gtkd->pixmap,
			   gtkd->wgc, TRUE,
			   (gint) x0, (gint) y0,
			   (gint) x1 - (gint) x0,
			   (gint) y1 - (gint) y0);
    }
    if (R_OPAQUE(col)) {
	SetColor(&gcol_outline, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

	SetLineType(dd, lty, lwd);

	gdk_draw_rectangle(gtkd->drawing->window,
			   gtkd->wgc, FALSE,
			   (gint) x0, (gint) y0,
			   (gint) x1 - (gint) x0,
			   (gint) y1 - (gint) y0);
	gdk_draw_rectangle(gtkd->pixmap,
			   gtkd->wgc, FALSE,
			   (gint) x0, (gint) y0,
			   (gint) x1 - (gint) x0,
			   (gint) y1 - (gint) y0);
    }
}

static void GTK_Circle(double x, double y, double r,
		       int col, int fill, double gamma, int lty, double lwd,
		       NewDevDesc *dd)
{
    GdkColor gcol_fill, gcol_outline;
    gint ix, iy, ir;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    if(!gtkd->drawing->window)
	return;

    ix = x - r;
    iy = y - r;
    ir = 2 * floor(r + 0.5);

    if (R_OPAQUE(fill)) {
	SetColor(&gcol_fill, fill);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	gdk_draw_arc(gtkd->drawing->window,
		     gtkd->wgc, TRUE,
		     ix, iy, ir, ir,
		     0, 23040);
	gdk_draw_arc(gtkd->pixmap,
		     gtkd->wgc, TRUE,
		     ix, iy, ir, ir,
		     0, 23040);
    }
    if (R_OPAQUE(col)) {
	SetColor(&gcol_outline, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

	SetLineType(dd, lty, lwd);

	gdk_draw_arc(gtkd->drawing->window,
		     gtkd->wgc, FALSE,
		     ix, iy, ir, ir,
		     0, 23040);
	gdk_draw_arc(gtkd->pixmap,
		     gtkd->wgc, FALSE,
		     ix, iy, ir, ir,
		     0, 23040);
    }
}

static void GTK_Line(double x1, double y1, double x2, double y2,
		     int col, double gamma, int lty, double lwd,
		     NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill;
    gint ix1, iy1, ix2, iy2;

    if(!gtkd->drawing->window)
	return;

    ix1 = (gint) x1;  iy1 = (gint) y1;
    ix2 = (gint) x2;  iy2 = (gint) y2;

    if (R_OPAQUE(col)) {
	SetColor(&gcol_fill, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	SetLineType(dd, lty, lwd);

	gdk_draw_line(gtkd->drawing->window,
		      gtkd->wgc, ix1, iy1, ix2, iy2);
	gdk_draw_line(gtkd->pixmap,
		      gtkd->wgc, ix1, iy1, ix2, iy2);
    }
}

static void GTK_Polyline(int n, double *x, double *y, 
			    int col, double gamma, int lty, double lwd,
			    NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill;
    GdkPoint *points;
    int i;

    if(!gtkd->drawing->window)
	return;

    points = g_new0(GdkPoint, n);

    for(i = 0; i < n; i++) {
	points[i].x = (gint16) x[i];
	points[i].y = (gint16) y[i];
    }

    if (R_OPAQUE(col)) {
	SetColor(&gcol_fill, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	SetLineType(dd, lty, lwd);

	gdk_draw_lines(gtkd->drawing->window,
		       gtkd->wgc, points, n);
	gdk_draw_lines(gtkd->pixmap,
		       gtkd->wgc, points, n);
    }

    g_free(points);
}

static void GTK_Polygon(int n, double *x, double *y, 
			int col, int fill, double gamma, int lty, double lwd,
			NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill, gcol_outline;
    GdkPoint *points;
    int i;
 
    if(!gtkd->drawing->window)
	return;

    points = g_new0(GdkPoint, n + 1);

    for(i = 0; i < n; i++) {
	points[i].x = (gint16) x[i];
	points[i].y = (gint16) y[i];
    }

    if (R_OPAQUE(fill)) {
	SetColor(&gcol_fill, fill);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	gdk_draw_polygon(gtkd->drawing->window,
			 gtkd->wgc, TRUE, points, n);
	gdk_draw_polygon(gtkd->pixmap,
			 gtkd->wgc, TRUE, points, n);
    }
    if (R_OPAQUE(col)) {
	SetColor(&gcol_outline, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_outline);

	SetLineType(dd, lty, lwd);

	gdk_draw_polygon(gtkd->drawing->window,
			 gtkd->wgc, FALSE, points, n);
	gdk_draw_polygon(gtkd->pixmap,
			 gtkd->wgc, FALSE, points, n);
    }

    g_free(points);
}

static void GTK_Text(double x, double y, char *str, 
		     double rot, double hadj, 
		     int col, double gamma, int font, double cex, double ps,
		     NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill;
    gint size;
    double rrot = DEG2RAD * rot;

    if(!gtkd->drawing->window)
	return;
    size = cex * ps + 0.5;
    SetFont(dd, font, size);
    gdk_gc_set_font(gtkd->wgc, gtkd->font);

    if (R_OPAQUE(col)) {
	SetColor(&gcol_fill, col);
	gdk_gc_set_foreground(gtkd->wgc, &gcol_fill);

	gdk_draw_text_rot(gtkd->drawing->window,
			  gtkd->font, gtkd->wgc,
			  (int) x, (int) y,
			  gtkd->windowWidth, gtkd->windowHeight,
			  str, strlen(str), rrot);
	gdk_draw_text_rot(gtkd->pixmap,
			  gtkd->font, gtkd->wgc, 
			  (int) x, (int) y,
			  gtkd->windowWidth, gtkd->windowHeight,
			  str, strlen(str), rrot);
    }
}


typedef struct _GTK_locator_info GTK_locator_info;

struct _GTK_locator_info {
    guint x;
    guint y;
    gboolean button1;
};

static void locator_button_press(GtkWidget *widget,
				 GdkEventButton *event,
				 gpointer user_data)
{
    GTK_locator_info *info;

    info = (GTK_locator_info *) user_data;

    info->x = event->x;
    info->y = event->y;
    if(event->button == 1)
	info->button1 = TRUE;
    else
	info->button1 = FALSE;

    gtk_main_quit();
}

static Rboolean GTK_Locator(double *x, double *y, NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GTK_locator_info *info;
    guint handler_id;
    gboolean button1;

    info = g_new(GTK_locator_info, 1);

    /* Flush any pending events */
    while(gtk_events_pending())
	gtk_main_iteration();

    /* connect signal */
    handler_id = gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "button-press-event",
				    (GtkSignalFunc) locator_button_press, (gpointer) info);

    /* run the handler */
    gtk_main();

    *x = (double) info->x;
    *y = (double) info->y;
    button1 = info->button1;

    /* clean up */
    gtk_signal_disconnect(GTK_OBJECT(gtkd->drawing), handler_id);
    g_free(info);

    if(button1)
	return TRUE;

    return FALSE;
}

static void GTK_Mode(gint mode, NewDevDesc *dd)
{
#ifdef XSYNC
    if(mode == 0)
	gdk_flush();
#else
    gdk_flush();
#endif
}

static void GTK_Hold(NewDevDesc *dd)
{
}


/* Device driver entry point */
Rboolean
GTKDeviceDriver(DevDesc *odd, char *display, double width, 
		double height, double pointsize)
{
    NewDevDesc *dd;
    int ps;
    gchar tmp[2];
    gint c, rbearing, lbearing;
    double max_rbearing, min_lbearing;
    gtkDesc *gtkd;

    dd = (NewDevDesc*) odd;

    if(!(gtkd = (gtkDesc *) malloc(sizeof(gtkDesc))))
	return FALSE;

    dd->deviceSpecific = (void *) gtkd;

    /* font loading */
    ps = pointsize;
    if(ps < 6 || ps > 24) ps = 12;
    ps = 2 * (ps / 2);
    gtkd->fontface = -1;
    gtkd->fontsize = -1;
    dd->startfont = 1; 
    dd->startps = ps;
    dd->startcol = 0;
    dd->startfill = NA_INTEGER;
    dd->startlty = LTY_SOLID; 
    dd->startgamma = 1;

    /* device driver start */
    if(!GTK_Open(dd, gtkd, display, width, height)) {
	free(gtkd);
	return FALSE;
    }

    dd->newDevStruct = 1;

    /* setup data structure */
    dd->open = GTK_Open;
    dd->close = GTK_Close;
    dd->activate = GTK_Activate;
    dd->deactivate = GTK_Deactivate;
    dd->size = GTK_Size;
    dd->newPage = GTK_NewPage;
    dd->clip = GTK_Clip;
    dd->strWidth = GTK_StrWidth;
    dd->text = GTK_Text;
    dd->rect = GTK_Rect;
    dd->circle = GTK_Circle;
    dd->line = GTK_Line;
    dd->polyline = GTK_Polyline;
    dd->polygon = GTK_Polygon;
    dd->locator = GTK_Locator;
    dd->mode = GTK_Mode;
    dd->hold = GTK_Hold;
    dd->metricInfo = GTK_MetricInfo;

    dd->left = 0;
    dd->right = gtkd->windowWidth;
    dd->bottom = gtkd->windowHeight;
    dd->top = 0;

    /* nominal character sizes */
    max_rbearing = 0;
    min_lbearing = 10000; /* just a random big number */
    for(c = 0; c <= 255; c++) {
	g_snprintf(tmp, 2, "%c", (gchar) c);
	gdk_string_extents(gtkd->font, tmp,
			   &lbearing, &rbearing,
			   NULL, NULL, NULL);
	if(lbearing < min_lbearing || c == 0)
	    min_lbearing = lbearing;
	if(rbearing > max_rbearing)
	    max_rbearing = rbearing;
    }

    dd->cra[0] = max_rbearing - min_lbearing;
    dd->cra[1] = (double) gtkd->font->ascent + (double) gtkd->font->descent;

    /* character addressing offsets */
    dd->xCharOffset = 0.4900;
    dd->yCharOffset = 0.3333;
    dd->yLineBias = 0.1;

    /* inches per raster unit */
    dd->ipr[0] = pixelWidth();
    dd->ipr[1] = pixelHeight();

    /* device capabilities */
    dd->canResizePlot= TRUE;
    dd->canChangeFont= FALSE;
    dd->canRotateText= TRUE;
    dd->canResizeText= TRUE;
    dd->canClip = FALSE; /* See comment in GTK_Clip */
    dd->canHAdj = 0;/* not better? {0, 0.5, 1} */
    dd->canChangeGamma = FALSE;

    /* gtk device description stuff */
    gtkd->cex = 1.0;
    gtkd->srt = 0.0;
    gtkd->resize = FALSE;

    dd->displayListOn = TRUE;

    /* finish */
    return TRUE;
}




Rboolean
GTKDeviceFromWidget(DevDesc *odd, char *widget, double width, double height, double pointsize)
{
    NewDevDesc *dd = (NewDevDesc *) odd;
    double ps = pointsize;
    gint iw, ih, w, h;
    GtkWidget *drawing = (GtkWidget *) widget;
    gtkDesc *gtkd;

    gchar tmp[2];
    gint  c, rbearing, lbearing;
    double max_rbearing, min_lbearing;

    GTK_DRAWING_AREA(drawing);

    if(!(gtkd = (gtkDesc *) malloc(sizeof(gtkDesc))))
	return FALSE;

    w = width; h = height;
    gtkd->window = NULL;
    gtkd->pixmap = NULL;
    gtkd->drawing = NULL;
    gtkd->wgc = NULL;
    gtkd->gcursor = NULL;
    gtkd->resize = 1;

    /* font loading */
    ps = pointsize;
    if(ps < 6 || ps > 24) ps = 12;
    ps = 2 * (ps / 2);
    gtkd->fontface = -1;
    gtkd->fontsize = -1;
    dd->startfont = 1; 
    dd->startps = ps;
    dd->startcol = 0;
    dd->startfill = NA_INTEGER;
    dd->startlty = LTY_SOLID; 
    dd->startgamma = 1;

    /* device driver start */
    {
#if 0
	GtkArg args[2];
	args[0].name = "GtkWidget::width";  
	args[1].name = "GtkWidget::height";  
	gtk_object_getv(GTK_OBJECT(drawing), 2, args);
	w = GTK_VALUE_INT(args[0]);
	if(w < 0)
	    w = 0;
	h = GTK_VALUE_INT(args[1]);
	if(h < 0)
	    h = 0;
#endif

	gtkd->drawing = drawing;
	if(GTK_WIDGET_REALIZED(gtkd->drawing))
	    gtk_widget_add_events(gtkd->drawing,
				  GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK);
	else
	    gtk_widget_set_events(gtkd->drawing,
				  GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK);
	gdk_rgb_init();
	gtk_widget_push_visual(gdk_rgb_get_visual());
	gtk_widget_push_colormap(gdk_rgb_get_cmap());

	gtkd->windowWidth = iw = w / pixelWidth();
	gtkd->windowHeight = ih = h / pixelHeight();
	/* connect to signal handlers, etc */
	gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "realize",
			   (GtkSignalFunc) realize_event, (gpointer) dd);


	SetColor(&gtkd->gcol_bg, R_RGB(255, 255, 255)); /*FIXME canvas color*/

	/* connect to signal handlers, etc */
	gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "configure_event",
			   (GtkSignalFunc) configure_event, (gpointer) dd);
	gtk_signal_connect(GTK_OBJECT(gtkd->drawing), "expose_event",
			   (GtkSignalFunc) expose_event, (gpointer) dd);

	dd->deviceSpecific = (void *) gtkd;

	/* initialise line params */
	gtkd->lty = -1;
	gtkd->lwd = -1;


	/* let other widgets use the default colour settings */
	gtk_widget_pop_visual();
	gtk_widget_pop_colormap();

	/* Set base font */
	if(!SetBaseFont(gtkd)) {
	    Rprintf("can't find X11 font\n");
	    return FALSE;
	}
    }

    dd->newDevStruct = 1;

    /* setup data structure */
    dd->open = GTK_Open;
    dd->close = GTK_Close;
    dd->activate = GTK_Activate;
    dd->deactivate = GTK_Deactivate;
    dd->size = GTK_Size;
    dd->newPage = GTK_NewPage;
    dd->clip = GTK_Clip;
    dd->strWidth = GTK_StrWidth;
    dd->text = GTK_Text;
    dd->rect = GTK_Rect;
    dd->circle = GTK_Circle;
    dd->line = GTK_Line;
    dd->polyline = GTK_Polyline;
    dd->polygon = GTK_Polygon;
    dd->locator = GTK_Locator;
    dd->mode = GTK_Mode;
    dd->hold = GTK_Hold;
    dd->metricInfo = GTK_MetricInfo;

    dd->left = 0;
    dd->right = gtkd->windowWidth;
    dd->bottom = gtkd->windowHeight;
    dd->top = 0;

    /* nominal character sizes */
    max_rbearing = 0;
    min_lbearing = 10000; /* just a random big number */
    for(c = 0; c <= 255; c++) {
	g_snprintf(tmp, 2, "%c", (gchar) c);
	gdk_string_extents(gtkd->font, tmp,
			   &lbearing, &rbearing,
			   NULL, NULL, NULL);
	if(lbearing < min_lbearing || c == 0)
	    min_lbearing = lbearing;
	if(rbearing > max_rbearing)
	    max_rbearing = rbearing;
    }

    dd->cra[0] = max_rbearing - min_lbearing;
    dd->cra[1] = (double) gtkd->font->ascent + (double) gtkd->font->descent;

    /* character addressing offsets */
    dd->xCharOffset = 0.4900;
    dd->yCharOffset = 0.3333;
    dd->yLineBias = 0.1;

    /* inches per raster unit */
    dd->ipr[0] = pixelWidth();
    dd->ipr[1] = pixelHeight();

    /* device capabilities */
    dd->canResizePlot= TRUE;
    dd->canChangeFont= FALSE;
    dd->canRotateText= TRUE;
    dd->canResizeText= TRUE;
    dd->canClip = FALSE; /* See comment in GTK_Clip */
    dd->canHAdj = 0;/* not better? {0, 0.5, 1} */
    dd->canChangeGamma = FALSE;

    /* gtk device description stuff */
    gtkd->cex = 1.0;
    gtkd->srt = 0.0;
    gtkd->resize = TRUE;

    dd->displayListOn = TRUE;

    /* finish */
    return TRUE;
}

#endif /* ifndef GTK2 */
