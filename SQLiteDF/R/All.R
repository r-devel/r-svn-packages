sopen <- function(filename) .Call("sopen", filename)

"sqlite.data.frame" <- function(x, name=NULL) 
    .Call("sdf_create_sdf", as.data.frame(x), name)

workspace <- NULL;
create.workspace <- function() workspace <- .Call("sdf_initialize");

.onLoad <- function(libname, pkgname) .Call("sdf_init_workspace");

.onUnload <- function(libpath) {
    .Call("sdf_finalize_workspace")
    library.dynam.unload("SQLiteDF", libpath)
}

length <- function(x) UseMethod("length");
length.default <- .Primitive("length");

# sdf functions
nrow.sqlite.data.frame <- function(x) .Call("sdf_get_rows", x);

# workspace functions
ls.sdf <- function(pattern=NULL) .Call("sdf_list_sdfs", pattern);
get.sdf <- function(name) .Call("sdf_get_sdf", name);

# S3 methods for sqlite.data.frame
names.sqlite.data.frame <- function(x) .Call("sdf_get_names", x);
length.sqlite.data.frame <- function(x) .Call("sdf_get_length", x);
dim.sqlite.data.frame <- function(x) 
    c(nrow.sqlite.data.frame(x), length.sqlite.data.frame(x))
