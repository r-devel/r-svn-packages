#include <gnome.h>
#include <Rversion.h>

/*
  Ensure that gnome is loaded. We may want to leave this to
  another package, so that people can add their own options.
*/


void loadGNOME()
{
	static gboolean gnome_initialized = FALSE;

	char **argv; 
	int argc = 1;
	argv = (char **) g_malloc(argc * sizeof(char *));
	argv[0] = g_strdup("R");
	if (!gnome_initialized) {
		gnome_init("R",
			   g_strdup_printf("%s.%s %s (%s-%s-%s)", 
					   R_MAJOR, R_MINOR,
					   R_STATUS, R_YEAR, R_MONTH, R_DAY),
			   argc, argv);
		printf("Initializing GNOME!\n");
		gnome_initialized=TRUE;
	}
}
