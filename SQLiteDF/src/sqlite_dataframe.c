#include <stdio.h>
#include <string.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

SEXP sdf_create_sdf(SEXP df, SEXP name) {
    SEXP ret, tmp;
    char *basename, *filename; 
    int file_idx = 0;

    /* check if name exists */
    if (IS_CHARACTER(name)) {
        int baselen = strlen(CHAR(STRING_ELT(name,0)));
        basename = R_alloc(baselen, sizeof(char));
        strcpy(basename, CHAR(STRING_ELT(name,0)));
        R_alloc(baselen + 9, sizeof(char)); /* <name>10000.db\0 */
        sprintf(filename, "%s.db", basename);
    } else if (name == R_NilValue) {
        basename = "data";
        filename = R_alloc(13, sizeof(char)); /* data10000.db\0 */
        sprintf(filename, "%s.db", basename);
        file_idx = 1;
    } else {
        /* how to error??? */
        return R_NilValue;
    }

    FILE *f;
    do {
        f = fopen(filename, "rb+");
        if (f == NULL) break;
        fclose(f);
        sprintf(filename, "%s%d.db", file_idx, basename);
    } while (file_idx < 10000);
    
    /* found a free number */
    if (file_idx < 10000) {
        /* create a sdf db file */
            /* attach to workspace.db */
            /* create tables */
            /* add rows to workspace */
        /* register */
    }
        
    PROTECT(ret = NEW_LOGICAL(1));
    LOGICAL(ret)[0] = (name == R_NilValue);
    UNPROTECT(1);
    return ret;
}


SEXP sopen(SEXP name) {
    char *filename;
    
    if (IS_CHARACTER(name)) {
        filename = CHAR(STRING_ELT(name,0));
        Rprintf("%s\n", filename);
    }
    /* sqlite3 *db;
    int res = sqlite3_open(filename, &db); */
    SEXP ret;
    PROTECT(ret = NEW_LOGICAL(1));
    LOGICAL(ret)[0] = IS_CHARACTER(name);
    /* sqlite3_close(db); */ 
    UNPROTECT(1);
    return ret;
}
