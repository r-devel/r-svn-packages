.First.lib <- function(libname, pkgname)
{
    library.dynam("gnomeDevice", pkgname, libname)
    .C("loadGNOME", PACKAGE = "gnomeDevice")
    .C("R_gtk_setEventHandler",  PACKAGE = "gnomeDevice") 
}

.Last.lib <- function(libname, pkgname)
{
    devices <- dev.list()
    gnome.devices <- devices[names(devices)=="GNOME"]
    if(length(gnome.devices) > 0) {
        dev.off(gnome.devices)
    }
}
