gnome <- function(display = "", width = 7, height = 7, pointsize = 12)
{
    .C("do_gnome", as.character(display), as.numeric(width),
       as.numeric(height), as.numeric(pointsize), PACKAGE = "gnomeDevice")
    return(invisible(TRUE))
}

GNOME <- gnome



