.First.lib <- function(libname, pkgname)
{
    library.dynam("gtkDevice", pkgname, libname)
    .C("loadGTK")
    .C("R_gtk_setEventHandler") 
}

.Last.lib <- function(libname, pkgname)
{
    devices <- dev.list()
    gtk.devices <- devices[names(devices)=="GTK"]
    if(length(gtk.devices) > 0) {
        dev.off(gtk.devices)
    }
}
