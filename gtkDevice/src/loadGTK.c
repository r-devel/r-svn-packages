#include "devGTK.h"
#include <R_ext/Rdynload.h>

/**
  Ensure that gtk is loaded. We may want to leave this to 
  another package, e.g.  RGtk, so that people can add their own 
  options.
 */
void loadGTK()
{
    char **argv; 
    int argc = 1;
    argv = (char **) g_malloc(argc * sizeof(char *));
    argv[0] = g_strdup("R");
    gtk_init(&argc, &argv);
}
