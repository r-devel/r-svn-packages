#include <gtk/gtk.h>


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rgraphics.h>
#include <Rdevices.h>
#include <R_ext/GraphicsDevice.h>
#include <R_ext/GraphicsEngine.h>

#include "devGTK.h"
/*
 The idea is quite common. We want to allow for rubber-banding on the screen.
 
 It is initiated from within R by calling the routine 
   R_initRubberBand()
 which returns a reference to the data structure which manages the rubber band on the screen.

 As the mouse is moved (with the button down), the rubber band and its information are updated.
 The R "program" can listen for events also and query the rubber band for its current location.

 See drawMotion.S in the inst/examples/ directory.

 This currently is only for gtk-1.2.*. This is because it relies on RGtk.
 Support for gtk2 will come when we release RGtk2.
*/

typedef enum  {MOTION, RELEASE, PRESS, MAX_EVENTS} EventTypes;

typedef struct {
  gboolean active;

  int x1, y1;
  int x2, y2;

  int isDrawn;

  GtkWidget *drawing_area;  
  GdkGC *gc; 

  guint eventId[3];
} RubberBandInfo;




RubberBandInfo *R_getRubberBandInfo(SEXP ref);
static gint button_press_event(GtkWidget *w, GdkEventButton *event, RubberBandInfo *info);
static SEXP makeRubberBandReference(RubberBandInfo *info);
static void initializeRubberBand(RubberBandInfo *info, GtkWidget *drawing_area);
void draw_rubber_band(RubberBandInfo *info, int x, int y);
void ConnectToEvents(GtkWidget *drawing_area, RubberBandInfo *info);
void setGC(GdkGC *gc, SEXP lwd, SEXP lty, SEXP color);




SEXP
R_getRubberBandCoordinates(SEXP ref)
{
   RubberBandInfo *info;
   SEXP ans = R_NilValue;

   info = R_getRubberBandInfo(ref);

   if(info) {
      ans = allocVector(INTSXP, 4);
      INTEGER(ans)[0] = info->x1; 
      INTEGER(ans)[1] = info->y1; 
      INTEGER(ans)[2] = info->x2; 
      INTEGER(ans)[3] = info->y2; 
   }

   return(ans);
}


static void
connect(RubberBandInfo *info)
{
  if(info->eventId[PRESS] == 0)
      info->eventId[PRESS] = gtk_signal_connect (GTK_OBJECT (info->drawing_area), "button_press_event",
                                                   (GtkSignalFunc) button_press_event, info);
}


static void
disconnect(RubberBandInfo *info)
{
    int i;
    for(i = 0; i < 2; i++){
       if(info->eventId[i] > 0) {
          gtk_signal_disconnect(GTK_OBJECT(info->drawing_area), info->eventId[i]);
          info->eventId[i] = 0;
       }
    }
}

SEXP
R_connectRubberBand(SEXP ref, SEXP con)
{
   RubberBandInfo *info;
   info = R_getRubberBandInfo(ref);
   if(!LOGICAL(con)[0])
      disconnect(info);
   else
      connect(info);


   return(ScalarLogical(TRUE));
}

RubberBandInfo *
R_getRubberBandInfo(SEXP ref)
{
   RubberBandInfo *info;
   if(TYPEOF(ref) != EXTPTRSXP) {
      PROBLEM "expecting an external pointer object for the rubber band information"
      ERROR;
   }

   if(R_ExternalPtrTag(ref) != Rf_install("RubberBandInfo")) {
      PROBLEM "expecting an external pointer object with tag RubberBandInfo"
      ERROR;
   }

   info = R_ExternalPtrAddr(ref);

   if(!info) {
      PROBLEM "NULL external pointer passed for rubber band info. This is probably left over from a previous session!"
      ERROR;
   }

   return(info);
}


SEXP
R_initRubberBand(SEXP s_w, SEXP lwd, SEXP lty, SEXP color, SEXP rconnect)
{
  GtkWidget *w = (GtkWidget *) R_ExternalPtrAddr(s_w);
  RubberBandInfo *info;

  info =  (RubberBandInfo *) calloc(1, sizeof(RubberBandInfo));
  if(!info) {
    PROBLEM "Cannot allocate %d bytes for a RubberBandInfo object", (int)sizeof(RubberBandInfo)
    ERROR;
  }

  initializeRubberBand(info, w);
  setGC(info->gc, lwd, lty, color);
  if(LOGICAL(rconnect)[0])
     connect(info);

  return(makeRubberBandReference(info));
}



void
releaseRubberBandInfo(SEXP ref)
{
  RubberBandInfo *info;
  info = (RubberBandInfo *) R_ExternalPtrAddr(ref);
  if(!info)
    return;

  free(info);
}


static SEXP
makeRubberBandReference(RubberBandInfo *info)
{
   SEXP ans;
   PROTECT(ans = R_MakeExternalPtr(info, Rf_install("RubberBandInfo"), R_NilValue));

   R_RegisterCFinalizer(ans, releaseRubberBandInfo);
   SET_CLASS(ans, mkString("RubberBandInfo"));

   UNPROTECT(1);
   return(ans);
}

static void
initializeRubberBand(RubberBandInfo *info, GtkWidget *drawing_area)
{
    int i;

    info->isDrawn = 0;
    info->active = TRUE;
    info->drawing_area = drawing_area;

    info->x1 = info->x2 = info->y1 = info->y2 = 0;

    for(i = 0; i < MAX_EVENTS; i++)
        info->eventId[i] = 0;

    info->gc = gdk_gc_new(drawing_area->window); 
    gdk_gc_copy( info->gc, drawing_area->style->black_gc); 
    gdk_gc_set_function(info->gc, GDK_XOR);

/*
    gdk_window_set_cursor(drawing_area->window, );
*/
}



static gint
motion_notify_event(GtkWidget *w, GdkEventMotion *event, RubberBandInfo *info)
{
  int x, y;
  GdkModifierType state;


  if(info->active == -1)
      return(TRUE);

  if (event->is_hint)
    gdk_window_get_pointer (event->window, &x, &y, &state);
  else
  {
    x = event->x;
    y = event->y;
    state = event->state;
  }

  if(GDK_BUTTON1_MASK) {
     draw_rubber_band(info, x, y);
  }

  return TRUE;
}


static gint
button_release_event(GtkWidget *w, GdkEventButton *event, RubberBandInfo *info)
{
   disconnect(info);

   if(info->isDrawn) {
      info->isDrawn = 0;
      draw_rubber_band(info, info->x2, info->y2);
      info->active = FALSE;
   }

    return TRUE;

}
static gint
button_press_event(GtkWidget *w, GdkEventButton *event, RubberBandInfo *info)
{
  if (event->button == 1) {

     info->eventId[MOTION] = gtk_signal_connect (GTK_OBJECT (w), "motion_notify_event",
                                                 (GtkSignalFunc) motion_notify_event, info);

     info->eventId[RELEASE] = gtk_signal_connect (GTK_OBJECT (w), "button_release_event",
                                                 (GtkSignalFunc) button_release_event, info);

     info->x1 = event->x;
     info->y1 = event->y;

     info->isDrawn = 0;
     draw_rubber_band(info, event->x, event->y);

  }
  return TRUE;
}



void
draw_rubber_band(RubberBandInfo *info, int x, int y)
{
    int w, h;

    if(info->isDrawn) {
       w = info->x2 - info->x1;
       h = info->y2 - info->y1;

       if(w != 0 && y != 0)
           gdk_draw_rectangle(info->drawing_area->window, info->gc, FALSE, info->x1, info->y1, w, h);
    }

    w = x - info->x1;
    h = y - info->y1;
    gdk_draw_rectangle(info->drawing_area->window, info->gc, FALSE, info->x1, info->y1, w, h);

    info->x2 = x; info->y2 = y;
    info->isDrawn = 1;
}


void
setGC(GdkGC *gc, SEXP lwd, SEXP lty, SEXP scolor)
{
   GdkColor color = {0, 0, 0};

   color.red = INTEGER(scolor)[0];
   color.green = INTEGER(scolor)[1];
   color.blue = INTEGER(scolor)[2];

   color.pixel = gdk_rgb_xpixel_from_rgb(color.red | color.green | color.blue);

   gdk_gc_set_foreground(gc, &color);
   gdk_gc_set_line_attributes(gc, INTEGER(lwd)[0], INTEGER(lty)[0], GDK_CAP_ROUND, GDK_JOIN_ROUND);
}
