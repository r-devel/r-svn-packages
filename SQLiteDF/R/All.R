"sqlite.data.frame" <- function(x, name=NULL) {
    if (inherits(x, "sqlite.data.frame")) x 
    else .Call("sdf_create_sdf", as.data.frame(x), name)
}

.onLoad <- function(libname, pkgname) .Call("sdf_init_workspace");

.onUnload <- function(libpath) {
    .Call("sdf_finalize_workspace")
    library.dynam.unload("SQLiteDF", libpath)
}

length <- function(x) UseMethod("length");
length.default <- .Primitive("length");

# workspace functions
lsSdf <- function(pattern=NULL) .Call("sdf_list_sdfs", pattern);
getSdf <- function(name) .Call("sdf_get_sdf", name);

readTableSdf <- function(filename, sep=",", quote="\"'", name=NULL, 
                            rownames, colnames)
    .Call("sdf_import_table", filename, name, sep, quote, rownames, colnames);

attachSdf <- function(sdf_filename, sdf_iname=NULL) 
    .Call("sdf_attach_sdf", sdf_filename, sdf_iname);
detachSdf <- function(iname) .Call("sdf_detach_sdf", iname);

# -------------------------------------------------------------------------
# S3 methods for sqlite.data.frame
# -------------------------------------------------------------------------
names.sqlite.data.frame <- function(x) .Call("sdf_get_names", x);
length.sqlite.data.frame <- function(x) .Call("sdf_get_length", x);
nrow.sqlite.data.frame <- function(x) .Call("sdf_get_row_count", x);
dim.sqlite.data.frame <- function(x) 
    c(nrow.sqlite.data.frame(x), length.sqlite.data.frame(x));
dimnames.sqlite.data.frame <- function(x) list(row.names(x), names(x))
"$.sqlite.data.frame" <- function(x, name) .Call("sdf_get_variable", x, name);
"[[.sqlite.data.frame" <- function(x, idx) {
    if (length(idx) != 1) stop("index must be a 1-element vector.");
    if (is.character(idx)) .Call("sdf_get_variable", x, idx)
    else if (is.numeric(idx)) {
        if (idx > length(x)) stop("subscript out of bounds")
        else .Call("sdf_get_variable", x, names(x)[idx])
    } else stop("don't know how to handle index.");
}
"[.sqlite.data.frame" <- function(x, row=NULL, col=NULL) {
#    if (row == NULL || col == NULL) { data.frame(NULL) }
#    if (missing(row)) row = NULL;
#    if (missing(col)) col = NULL; 
    .Call("sdf_get_index", x, row, col); 
}

# -------------------------------------------------------------------------
# S3 methods for sqlite.vector
# -------------------------------------------------------------------------
"[.sqlite.vector" <- function(x, idx) {
    # temporary, better to be in C because assumption is length(x) is large
    if (is.numeric(idx) && all(idx <= 0)) idx <- (1:length(x))[idx];
    .Call("sdf_get_variable_index", x, idx);
}
length.sqlite.vector <- function(x) .Call("sdf_get_variable_length", x);
# methods to "coerce" to ordinary vectors
as.double.sqlite.vector <- function(x) as.double(x[1:length(x)])
as.character.sqlite.vector <- function(x) as.character(x[1:length(x)])
as.logical.sqlite.vector <- function(x) as.logical(x[1:length(x)])
as.integer.sqlite.vector <- function(x) as.integer(x[1:length(x)])
