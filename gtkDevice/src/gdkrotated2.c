/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1998-2002   Lyndon Drake
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

#include <math.h>
#include "gdkrotated.h"

static GdkPixmap* 
draw_layout_to_bitmap (PangoLayout *layout)
    /* 
       Render layout to a newly allocated bitmap
    */
{
    GdkPixmap *pixmap;
    PangoRectangle rect;
    GdkGC *gc;
    
    pango_layout_get_pixel_extents (layout, NULL, &rect);
    pixmap = gdk_pixmap_new (NULL, rect.width, rect.height, 1);
    gc = gdk_gc_new(pixmap);

    /* Clear the background */
    gdk_gc_set_function(gc, GDK_CLEAR);
    gdk_draw_rectangle(pixmap, gc, TRUE, 0, 0,
		       rect.width, rect.height);

    /* Render the layout onto the bitmap */
    gdk_gc_set_function(gc, GDK_SET);
    gdk_draw_layout(pixmap, gc, 0, 0, layout);

    g_object_unref(gc);
    return pixmap;
}

static GdkPixbuf*
gdk_pixbuf_rotate(GdkPixbuf *source, double angle)
    /* 
       Returns a newly allocated pixbuf that contains the source
       pixbuf, rotated through "angle" radians.
    */
{
    GdkPixbuf *rpixbuf;
    int width, height, rowstride, rwidth, rheight, rrowstride;
    double costheta, sintheta;
    int i, j, k;
    guchar *pixels, *rpixels;
    int nchan, rnchan;

    g_return_val_if_fail(gdk_pixbuf_get_bits_per_sample(source)==8, NULL);
    g_return_val_if_fail(gdk_pixbuf_get_colorspace(source)==GDK_COLORSPACE_RGB,
			 NULL);

    width = gdk_pixbuf_get_width(source);
    height = gdk_pixbuf_get_height(source);

    costheta = floor(cos(angle) * 1000.0 + 0.5) / 1000.0;
    sintheta = floor(sin(angle) * 1000.0 + 0.5) / 1000.0;
    
    /* Size of new pixbuf */
    rwidth = (int) (height * fabs(sintheta) + width * fabs(costheta));
    rheight = (int) (width * fabs(sintheta) + height * fabs(costheta));

    /* Create a new pixbuf of the required size and fill it
       with transparent black */
    rpixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, TRUE, 8, rwidth, rheight);
    gdk_pixbuf_fill(rpixbuf, 0x00000000);

    /* 
       Pixel-wise copying
    */
    rowstride = gdk_pixbuf_get_rowstride(source);
    pixels = gdk_pixbuf_get_pixels(source);
    nchan = gdk_pixbuf_get_n_channels(source);
    
    rrowstride = gdk_pixbuf_get_rowstride(rpixbuf);
    rpixels = gdk_pixbuf_get_pixels(rpixbuf);
    rnchan = gdk_pixbuf_get_n_channels(rpixbuf);

    for (i = 0; i < rwidth; i++)
    {
	double x, y;
	
	x = (i + 0.5) - rwidth/2;
	for (j = 0; j < rheight; j++)
	{
	    double u, v;
	    int p, q;
	    
	    y = rheight/2 - (j + 0.5);

	    u = x * costheta + y * sintheta;
	    v = - x * sintheta + y * costheta;

	    p = (int) (u + width/2);
	    q = (int) (height/2 - v);

	    if (p >= 0 && p < width && q >= 0 && q < height) {
		for (k = 0; k < nchan; k++) {
		    rpixels[j * rrowstride + rnchan * i + k]
			= pixels[q * rowstride + nchan * p + k];
		}
	    }
	}
    }

    return(rpixbuf);
}

static void repaint_pixbuf (GdkPixbuf *pixbuf, guchar red, guchar green,
			    guchar blue, GdkGC *gc)
    /* Repaints the pixels in pixbuf matching red, green and blue with
       the foreground cover of the graphics context gc */
{
    int i, j;
    GdkGCValues values, values2;
    GdkColormap *colormap;
    GdkColor fgcolor;
    guchar *pixels;
    int width, height, rowstride, nchan, index;

    gdk_gc_get_values(gc, &values);
    /* 
       Under X11, gdk_gc_get_values only fills in the pixel field
       for the foreground colour. We need to finish the job.  This
       is a bug in gdk.
    */
    colormap = gdk_gc_get_colormap(gc);
    gdk_colormap_query_color(colormap, values.foreground.pixel,
			     &(values.foreground));

    pixels = gdk_pixbuf_get_pixels(pixbuf);
    width = gdk_pixbuf_get_width(pixbuf);
    height = gdk_pixbuf_get_height(pixbuf);
    rowstride = gdk_pixbuf_get_rowstride(pixbuf);
    nchan = gdk_pixbuf_get_n_channels(pixbuf);

    for (i = 0; i < width; i ++) {
	for (j = 0; j < height; j++) {
	    index = j * rowstride + i * nchan;
	    if (pixels[index] == red &&	
		pixels[index + 1] == green && 
		pixels[index + 2] == blue) 
	    {
		pixels[index] =     values.foreground.red >> 8;
		pixels[index + 1] = values.foreground.green >> 8;
		pixels[index + 2] = values.foreground.blue >> 8;
	    }
	}
    }


}

void gdk_draw_layout_rot(GdkDrawable *drawable,
			 GdkGC *gc,
			 int x,
			 int y,
			 PangoLayout *layout,
			 double angle)
/* 
   GDK does not handle rotated text, so we have to render the text
   onto a pixmap, rotate it, and then render the result back onto the
   drawable.

   Since the rotation involves pixel-wise manipulation of the pixmap,
   we copy it over to a pixbuf - which is a client side object
   in X terminology - for greater efficiency.
*/
{
    PangoRectangle lrect;	
    PangoLayoutLine *line;

    g_return_if_fail(GDK_IS_GC(gc));
    g_return_if_fail(GDK_IS_DRAWABLE(drawable));
    
    line = pango_layout_get_line(layout, 0);
    pango_layout_line_get_pixel_extents(line, NULL, &lrect);    
    
    if(angle == 0.0) {
	gdk_draw_layout(drawable, gc,
			x - PANGO_LBEARING(lrect),
			y - PANGO_ASCENT(lrect), layout);
    }
    else {
	PangoRectangle rect;
	GdkPixmap *pixmap;
	GdkPixbuf *pixbuf, *tpixbuf, *rpixbuf;
	double xpivot, ypivot;

	/* Render layout to a bitmap and then copy it to a pixbuf */
	pango_layout_get_pixel_extents(layout, NULL, &rect);
	pixmap = draw_layout_to_bitmap(layout);
	pixbuf = gdk_pixbuf_get_from_drawable(NULL, pixmap, NULL,
					      0, 0, 0, 0, 
					      rect.width, rect.height);
	g_object_unref(pixmap);
	if (!pixbuf)
	    return;

	/* 
	   Make a copy of pixbuf with an alpha channel, making
	   the black pixels transparent and repainting the white
	   pixels with the foreground colour
	*/
	tpixbuf = gdk_pixbuf_add_alpha(pixbuf, TRUE, 0, 0, 0);
	repaint_pixbuf(tpixbuf, 255, 255, 255, gc);
	g_object_unref(pixbuf);

	rpixbuf = gdk_pixbuf_rotate (tpixbuf, angle);
	if (rpixbuf) {
	    double u, v;
	    int xnew, ynew;

	    /* 
	       u, v are coordinates of the centre of the text,
	       with respect to the pivot point 
	    */
	    u = (rect.width/2 - PANGO_LBEARING(lrect)) * cos (angle) - 
		(rect.height/2 - PANGO_ASCENT(lrect)) * sin (angle); 
	    v = (rect.width/2 - PANGO_LBEARING(lrect)) * sin (angle) + 
		(rect.height/2 - PANGO_RBEARING(lrect)) * cos (angle);

	    /* 
	       xnew, ynew are coordinates of the top left hand corner
	       of rpixbuf
	    */
	    xnew = x + u - 3 * gdk_pixbuf_get_width(rpixbuf)/2;
	    ynew = y + v - 3 * gdk_pixbuf_get_height(rpixbuf)/2;

	    gdk_pixbuf_render_to_drawable (rpixbuf, drawable, gc, 
					   0, 0, xnew, ynew, -1, -1, 
					   GDK_RGB_DITHER_NONE, 0, 0);
	    g_object_unref(rpixbuf);
	}
	g_object_unref(tpixbuf);
    }
}

#endif /* ifdef GTK2 */









