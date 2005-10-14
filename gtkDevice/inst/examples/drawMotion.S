library(RGtk)

d = gtkDrawingArea()
d$SetEvents(as.integer(sum(2^c(2, 3, 4, 8, 9))))
w = gtkWindow()
w$SetUsize(600, 600)

w$Add(d)

library(gtkDevice)

dev = asGtkDevice(d)

if(TRUE) {
  plot(1:10, lty = 2, lwd = 3, cex = 3, col = "red", main  = "Testing fonts", font.main = 3)
} else {
    # Delay the call to plot until the widget is actually displayed.
  d$AddCallback("expose_event", function(...) { plot(1:10)})
}


if(FALSE) {
  i = rubberBand(d)
  # move the rubber band.
  getRubberBandCoordinates(i)
}

