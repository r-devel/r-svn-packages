sopen <- function(filename) .Call("sopen", filename)

"sqlite.data.frame" <- function(x, name=NULL) 
    .Call("sdf_create_sdf", as.data.frame(x), name)

workspace <- NULL;
create.workspace <- function() workspace <- .Call("sdf_initialize");


.onUnload <- function(libpath) library.dynam.unload("SQLiteDF", libpath)
