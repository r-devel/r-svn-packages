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

#ifdef GTK2

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <locale.h>

#include <R.h>
#include <Rinternals.h>
#include <Rgraphics.h>
#include <R_ext/GraphicsDevice.h>
#include <R_ext/GraphicsEngine.h>

#include "devGTK.h"
#include "gdkrotated.h"

#define CURSOR		GDK_CROSSHAIR		/* Default cursor */
#define MM_PER_INCH	25.4			/* mm -> inch conversion */


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

/*
static char *fontname_R6 = "-adobe-helvetica-%s-%s-*-*-*-%d-*-*-*-*-*-*";
static char *symbolname = "-adobe-symbol-*-*-*-*-*-%d-*-*-*-*-*-*";
static char *fixedname = "fixed";
*/

static PangoStyle  slant[] = {PANGO_STYLE_NORMAL, PANGO_STYLE_OBLIQUE};
static PangoWeight weight[] = {PANGO_WEIGHT_NORMAL, PANGO_WEIGHT_BOLD};
static PangoFontDescription *fontdesc;

struct _FontMetricCache {
    gint ascent[255];
    gint descent[255];
    gint width[255];
    gint font_ascent;
    gint font_descent;
    gint max_width;
};

static void
text_extents (PangoFont   *font,
	      PangoContext *context,
	      const gchar *text,
	      gint         text_length,
	      gint        *lbearing,
	      gint        *rbearing,
	      gint        *width,
	      gint        *ascent,
	      gint        *descent)
{
    PangoLayout *layout;
    PangoLayoutLine *line;
    PangoRectangle rect, lrect;
    char *utf8;
    gsize bytes_read, bytes_written;

    layout = pango_layout_new(context);
    utf8 = g_locale_to_utf8(text, text_length, &bytes_read, &bytes_written,
			    NULL);
    pango_layout_set_text(layout, utf8, bytes_written);
    pango_layout_get_pixel_extents(layout, NULL, &rect);

    line = pango_layout_get_line(layout, 0);
    pango_layout_line_get_pixel_extents(line, NULL, &lrect);

    if(ascent)
	*ascent = PANGO_ASCENT(lrect);
    if(descent)
	*descent = rect.height - PANGO_ASCENT(lrect);
    if(width)
	*width = rect.width;
    if(lbearing)
	*lbearing = PANGO_LBEARING(rect);
    if(rbearing)
	*rbearing = PANGO_RBEARING(rect);

    g_free(utf8);
    g_object_unref(layout);
}


static PangoFont *RGTKLoadFont(PangoFontDescription *fontdesc,
			       gtkDesc *gtkd)
{
    /* static GHashTable *font_htab = NULL; */
    PangoFont *tmp_font;
    PangoContext *context;

    /* FIXME: DISABLING FONT HASHING
    if(!font_htab) {
	font_htab = g_hash_table_new((GHashFunc) pango_font_description_hash,
				     (GEqualFunc)pango_font_description_equal);
    }
    
    tmp_font = g_hash_table_lookup(font_htab, (gconstpointer) fontdesc);

    if(!tmp_font) {
	PangoContext *context;
	context = gtk_widget_get_pango_context(GTK_WIDGET(gtkd->drawing));
	tmp_font = pango_context_load_font(context, fontdesc);
	if (!tmp_font) {
	    g_hash_table_insert(font_htab, (gpointer) fontdesc,
				(gpointer) tmp_font);
        }
    }
    */

    context = gtk_widget_get_pango_context(GTK_WIDGET(gtkd->drawing));
    tmp_font = pango_context_load_font(context, fontdesc);
    
    return tmp_font;
}

static gint SetBaseFont(gtkDesc *gtkd)
{
    PangoFontDescription *fontdesc;

    gtkd->fontface = 1; /* Role of fontface ? Not used here */
    gtkd->fontsize = 12;
    gtkd->usefixed = 0;

    fontdesc = pango_font_description_new();
    pango_font_description_set_family(fontdesc, "helvetica");
    pango_font_description_set_style(fontdesc, slant[0]);
    pango_font_description_set_weight(fontdesc, weight[0]);
    pango_font_description_set_size(fontdesc, gtkd->fontsize * PANGO_SCALE);

    gtkd->font = RGTKLoadFont(fontdesc, gtkd);
    pango_font_description_free(fontdesc);

    if(gtkd->font)
	return TRUE;

    gtkd->usefixed = 1;
    fontdesc = pango_font_description_from_string("fixed");
    gtkd->font = RGTKLoadFont(fontdesc, gtkd);
    pango_font_description_free(fontdesc);

    if(gtkd->font)
	return TRUE;

    return FALSE;
}

#define SMALLEST 8
#define LARGEST 24


static void SetFont(NewDevDesc *dd, gint face, gint size)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    PangoFont *tmp_font;
    PangoFontDescription *fontdesc;
    PangoContext *context;

    if (gtkd->usefixed)
	return;
    
    if (face < 1 || face > 5)
	face = 1;

    size = CLAMP (size, SMALLEST, LARGEST);

    gtkd->fontface = face;
    gtkd->fontsize = size;

    fontdesc = pango_font_description_new();
    if (face == 5) {
	pango_font_description_set_family(fontdesc, "symbol");
	pango_font_description_set_size(fontdesc, PANGO_SCALE * size);
    }
    else {
	pango_font_description_set_family(fontdesc, "helvetica");
	pango_font_description_set_weight(fontdesc, weight[(face-1)%2]);
	pango_font_description_set_style(fontdesc, slant[((face-1)/2)%2]);
	pango_font_description_set_size(fontdesc, PANGO_SCALE * size);
    }
    context = gtk_widget_get_pango_context (gtkd->drawing);
    pango_context_set_font_description (context, fontdesc);
    pango_font_description_free(fontdesc);

    if ((tmp_font = RGTKLoadFont(fontdesc, gtkd))) {
	gtkd->font = tmp_font;
    }
    /* FIXME: Should indicate failure */
}

static void SetColor(GdkColor *gcol, int color)
/* 
   Convert an R color to a GDK color
   NB.  The pixel value is not set, so in some instances
   you will need to call gdk_rgb_get_color() with a suitable
   colormap after calling this function
*/
{
    gcol->red = R_RED(color) << 8;
    gcol->green = R_GREEN(color) << 8;
    gcol->blue = R_BLUE(color) << 8;
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

static gboolean realize_event(GtkWidget *widget, NewDevDesc *dd)
{
    g_return_val_if_fail(dd != NULL, FALSE);
    
    return(initialize(dd));
}


static gboolean initialize(NewDevDesc *dd)
{
    gtkDesc *gtkd;
    GdkGCClass *gcclass;

    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_val_if_fail(gtkd != NULL, FALSE);
    g_return_val_if_fail(gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail(GTK_IS_DRAWING_AREA(gtkd->drawing), FALSE);

    /* create gc */
    gtkd->wgc = gdk_gc_new(gtkd->drawing->window);
    gcclass = GDK_GC_GET_CLASS(gtkd->wgc);

    /* set the cursor */
    gtkd->gcursor = gdk_cursor_new(GDK_CROSSHAIR);
    gdk_window_set_cursor(gtkd->drawing->window, gtkd->gcursor);

    /* set window bg */
    gdk_rgb_find_color(gtk_widget_get_colormap(gtkd->drawing),
		       &gtkd->gcol_bg);
    gdk_window_set_background(gtkd->drawing->window, &gtkd->gcol_bg);

    if(gtkd->wgc)
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gtkd->gcol_bg);
    
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

static gint configure_event(GtkWidget *widget,
			    GdkEventConfigure *event,
			    gpointer data)
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
    if((GTK_WIDGET_REALIZED(gtkd->drawing)) && 
       ((gtkd->windowWidth != event->width) || 
	(gtkd->windowHeight != event->height)))
    {
	gtkd->windowWidth = event->width;
	gtkd->windowHeight = event->height;

	gtkd->resize = TRUE;
    }

    return FALSE;
}

static void GTK_resize (NewDevDesc *dd);

static gint expose_event (GtkWidget *widget,
			  GdkEventExpose *event,
			  NewDevDesc *dd)
{
    gtkDesc *gtkd;
    
    g_return_val_if_fail (dd != NULL, FALSE);
    gtkd = (gtkDesc *) dd->deviceSpecific;
    g_return_val_if_fail (gtkd != NULL, FALSE);
    g_return_val_if_fail (gtkd->drawing != NULL, FALSE);
    g_return_val_if_fail (GTK_IS_DRAWING_AREA (gtkd->drawing), FALSE);

    if(gtkd->wgc == NULL) {
	realize_event (gtkd->drawing, dd);
    }
    else if(gtkd->resize != 0) {
	GTK_resize (dd); 
    }
    else if(gtkd->pixmap) {
	gdk_draw_pixmap (gtkd->drawing->window, gtkd->wgc, gtkd->pixmap,
			 event->area.x, event->area.y,
			 event->area.x, event->area.y,
			 event->area.width, event->area.height);
    }
    else {
	GEplayDisplayList((GEDevDesc*) Rf_GetDevice(devNumber((DevDesc*)dd)));
    }
    
    return FALSE;
}

static gint delete_event(GtkWidget *widget, GdkEvent *event, NewDevDesc *dd)
{
    g_return_val_if_fail (dd != NULL, FALSE);

    Rf_KillDevice ((DevDesc*) Rf_GetDevice (devNumber ((DevDesc*) dd)));
    
    return TRUE;
}

/*
static void tb_activate_cb(GtkWidget *widget, NewDevDesc *dd)
{
    g_return_if_fail (dd != NULL);
    
    selectDevice (devNumber ((DevDesc*)dd));
}

static void tb_close_cb(GtkWidget *widget, NewDevDesc *dd)
{
    g_return_if_fail(dd != NULL);
    
    Rf_KillDevice ((DevDesc*) Rf_GetDevice (devNumber ((DevDesc*) dd)));
}
*/

/* create window etc */
static Rboolean GTK_Open(NewDevDesc *dd, gtkDesc *gtkd, char *dsp, double w,
			 double h)
{
    gint iw, ih;

    /* initialise pointers */
    gtkd->drawing = NULL;
    gtkd->wgc = NULL;
    gtkd->gcursor = NULL;

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
    /* FIXME: There should be a canvas argument to gtk() which is
       used here */
    SetColor(&gtkd->gcol_bg, R_RGB(255, 255, 255));

    /* place and realize the drawing area */
    gtk_container_add(GTK_CONTAINER(gtkd->window), gtkd->drawing);
    gtk_widget_realize(gtkd->drawing);

    /* connect to signal handlers, etc */
    g_signal_connect(G_OBJECT (gtkd->drawing), "configure_event",
		     G_CALLBACK (configure_event), dd);
    g_signal_connect(G_OBJECT (gtkd->drawing), "expose_event",
		     G_CALLBACK (expose_event), dd);
    g_signal_connect(G_OBJECT(gtkd->window), "delete_event",
		     G_CALLBACK (delete_event), dd);
    
    /* show everything */
    gtk_widget_show_all(gtkd->window);
    
    /* initialise line params */
    gtkd->lty = -1;
    gtkd->lwd = -1;

    /* create offscreen drawable */
    gtkd->pixmap = gdk_pixmap_new(gtkd->drawing->window,
				  gtkd->windowWidth, gtkd->windowHeight,
				  -1);
    gdk_gc_set_rgb_fg_color(gtkd->wgc, &gtkd->gcol_bg);
    gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
		       gtkd->windowWidth, gtkd->windowHeight);

    /* Set base font */
    if(!SetBaseFont(gtkd)) {
	Rprintf("can't find X11 font\n");
	return FALSE;
    }

    /* we made it! */
    return TRUE;
}


static double GTK_StrWidth (char *str, int font,
			    double cex, double ps, NewDevDesc *dd)
{
    int size, width;
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    size = cex * ps + 0.5;
    SetFont(dd, font, size);
    
    text_extents(gtkd->font, gtk_widget_get_pango_context(gtkd->drawing),
		 str, strlen(str),
		 0, 0, &width, 0, 0);

    return (double) width;
}


static void GTK_MetricInfo (int c, int font, double cex, double ps,
			    double* ascent, double* descent,
			    double* width, NewDevDesc *dd)
{
    gint size;
    gint lbearing, rbearing, iascent, idescent, iwidth;
    gint maxwidth;
    gchar text[2];
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    PangoFontMetrics *font_metrics;
    gchar *lang;

    /* 
       FIXME: lbearing rbearing unused
       FIXME: Get language right
    */
    lang = setlocale(LC_MESSAGES, NULL);
    font_metrics = pango_font_get_metrics(gtkd->font, 
					  pango_language_from_string(lang));

    size = cex * ps + 0.5;
    SetFont(dd, font, size);

    if(c == 0) {
	maxwidth = 0;

	for(c = 0; c <= 255; c++) {
	    g_snprintf(text, 2, "%c", (gchar) c);
	    text_extents(gtkd->font, 
			 gtk_widget_get_pango_context(gtkd->drawing),
			 text, strlen(text), 0, 0, &iwidth, 0, 0);
	    if (iwidth > maxwidth)
		maxwidth = iwidth;
	}

	*ascent = (double) pango_font_metrics_get_ascent(font_metrics);
	*descent = (double) pango_font_metrics_get_descent(font_metrics);
	*width = maxwidth;
    }
    else 
    {
	g_snprintf(text, 2, "%c", (gchar) c);
	text_extents(gtkd->font, 
		     gtk_widget_get_pango_context(gtkd->drawing),
		     text, strlen(text),
		     0, 0, &iwidth, &iascent, &idescent);
	*ascent = iascent;
	*descent = idescent;
	*width = iwidth;
    }
}

/* set clipping */
static void GTK_Clip (double x0, double x1, double y0, double y1,
		      NewDevDesc *dd)
{
    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;

    gtkd->clip.x = MIN(x0, x1);
    gtkd->clip.width = abs(x0 - x1);

    gtkd->clip.y = MIN(y0, y1);
    gtkd->clip.height = abs(y0 - y1);

    gdk_gc_set_clip_rectangle(gtkd->wgc, &gtkd->clip);
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
		gdk_gc_set_rgb_fg_color(gtkd->wgc, &gtkd->gcol_bg);
		gdk_draw_rectangle(gtkd->pixmap, gtkd->wgc, TRUE, 0, 0,
				   gtkd->windowWidth, gtkd->windowHeight);
	    }
	}
    }

    GEplayDisplayList ((GEDevDesc*) Rf_GetDevice(devNumber((DevDesc*)dd)));
}

/* clear the drawing area */
static void GTK_NewPage (int fill, double gamma, NewDevDesc *dd)
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
	gdk_rgb_find_color(gtk_widget_get_colormap(gtkd->drawing),
			   &gtkd->gcol_bg);
	gdk_window_set_background(gtkd->drawing->window, &gtkd->gcol_bg);
	gtkd->fill = fill;
    }

    gdk_window_clear(gtkd->drawing->window);

    if(gtkd->wgc) {
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gtkd->gcol_bg);
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
	g_object_unref(gtkd->pixmap);
    
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

    devnum = devNumber((DevDesc*)dd) + 1;

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

    devnum = devNumber((DevDesc*)dd) + 1;

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_outline);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_outline);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);

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
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);

	gdk_draw_polygon(gtkd->drawing->window,
			 gtkd->wgc, TRUE, points, n);
	gdk_draw_polygon(gtkd->pixmap,
			 gtkd->wgc, TRUE, points, n);
    }
    if (R_OPAQUE(col)) {
	SetColor(&gcol_outline, col);
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_outline);

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
    PangoContext *context;
    PangoLayout *layout;

    gtkDesc *gtkd = (gtkDesc *) dd->deviceSpecific;
    GdkColor gcol_fill;
    gint size;
    double rrot = DEG2RAD * rot;

    if(!gtkd->drawing->window)
	return;
    size = cex * ps + 0.5;
    SetFont(dd, font, size);

    context = gtk_widget_get_pango_context(gtkd->drawing);
    layout = pango_layout_new(context);
    /* FIXME String should be UTF-8 */
    pango_layout_set_text(layout, str, -1);

    /* gdk_gc_set_font(gtkd->wgc, gtkd->font); */

    if (R_OPAQUE(col)) {
	SetColor(&gcol_fill, col);
	gdk_gc_set_rgb_fg_color(gtkd->wgc, &gcol_fill);
	
	gdk_draw_layout_rot(gtkd->drawing->window,
			    gtkd->wgc,
			    (int) x,
			    (int) y,
			    layout, 
			    rrot);
	gdk_draw_layout_rot(gtkd->pixmap,
			    gtkd->wgc,
			    (int) x,
			    (int) y,
			    layout, 
			    rrot);
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
    gint cumwidth, c, rbearing, lbearing, ascent, descent;
    /* fixme lbearing rbearing unused */
    double max_rbearing, min_lbearing;
    PangoFontMetrics *font_metrics;
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
    cumwidth = 0;
    max_rbearing = 0;
    min_lbearing = 10000; /* just a random big number */
    for(c = 0; c <= 255; c++) {
	gint ilbearing, irbearing;
	g_snprintf(tmp, 2, "%c", (gchar) c);
	text_extents(gtkd->font, 
		     gtk_widget_get_pango_context(gtkd->drawing),
		     tmp, strlen(tmp),
		     &ilbearing, &irbearing, 0, 0, 0);
	if(ilbearing < min_lbearing || c == 0)
	    min_lbearing = ilbearing;
	if(irbearing > max_rbearing || c == 0)
	    max_rbearing = irbearing;
    }
    font_metrics = pango_font_get_metrics(gtkd->font, /* FIXME: language */
					  pango_language_from_string("en"));
    ascent = pango_font_metrics_get_ascent(font_metrics);
    descent = pango_font_metrics_get_descent(font_metrics);
    dd->cra[1] = PANGO_PIXELS(ascent) + PANGO_PIXELS(descent);
    pango_font_metrics_unref(font_metrics);
    /* dd->cra[0] ??? */

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
    dd->canClip = FALSE;/* FIXME: really? */
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
    /* FIXME: lbearing, rbearing unused */
    gint cumwidth, c, rbearing, lbearing, ascent, descent;
    double max_rbearing, min_lbearing;
    PangoFontMetrics *font_metrics;

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
    cumwidth = 0;
    max_rbearing = 0;
    min_lbearing = 10000; /* just a random big number */
    for(c = 0; c <= 255; c++) {
	gint ilbearing, irbearing;
	g_snprintf(tmp, 2, "%c", (gchar) c);
	text_extents(gtkd->font, 
		     gtk_widget_get_pango_context(gtkd->drawing),
		     tmp, strlen(tmp),
		     &ilbearing, &irbearing, 0, 0, 0);
	if(ilbearing < min_lbearing || c == 0)
	    min_lbearing = ilbearing;
	if(irbearing > max_rbearing || c == 0)
	    max_rbearing = irbearing;
    }
    font_metrics = pango_font_get_metrics(gtkd->font, /* FIXME: language */
					  pango_language_from_string("en"));
    ascent = pango_font_metrics_get_ascent(font_metrics);
    descent = pango_font_metrics_get_descent(font_metrics);
    dd->cra[1] = PANGO_PIXELS(ascent) + PANGO_PIXELS(descent);
    pango_font_metrics_unref(font_metrics);
    /* FIXME cra[0]; */

    /*
    max_rbearing = 0;
    min_lbearing = 10000;
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
    */

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
    dd->canClip = FALSE;/* FIXME: really? */
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

#endif /* ifdef GTK2 */














