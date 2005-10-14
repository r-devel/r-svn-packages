library(gtkDevice)
library(RGtk)

win = gtkWindow()
win$SetUsize(300, 300)
d = gtkDrawingArea()

win$Add(d)

dev = asGtkDevice(d)
plot(1:10)

# band = .Call("R_initRubberBand", d)

