gtk <- function(display = "", width = 7, height = 7, pointsize = 12)
{
    .C("do_GTK", as.character(display), as.numeric(width),
       as.numeric(height), as.numeric(pointsize))
    return(invisible(TRUE))
}

GTK <- gtk

asGtkDevice <- function(widget, width = 300, height = 300, pointsize = 12)
{
    require("RGtk")
    if(!inherits(widget, "GtkDrawingArea")) {
        stop("Widget being used as a Gtk Device must be or extend the GtkDrawingWidget class")
    }

    .C("do_asGTKDevice", widget, as.numeric(width), as.numeric(height),
       as.numeric(pointsize))
}



