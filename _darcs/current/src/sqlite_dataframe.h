#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

#ifndef __SQLITE_DATAFRAME__
#define __SQLITE_DATAFRAME__

#include "sqlite3.h"

#define WORKSPACE_COLUMNS 3

int _is_r_sym(char *sym);
int _file_exists(char *filename);
int _empty_callback(void *data, int ncols, char **row, char **cols);
char *_r2iname(char *internal_name, char *filename);
char *_fixname(char *rname);
char *_get_full_pathname2(char *relpath);
int _sqlite_error(int res);
SEXP _getListElement(SEXP list, char *varname);
int _expand_buf(int i, int size);
const char *_get_column_type(const char *class, int type);
int _add_sdf1(char *filename, char *internal_name);
SEXP _create_sdf_sexp(char *iname);
int _get_factor_levels1(const char *iname, const char *varname, SEXP var);
SEXP _shrink_vector(SEXP vec, int len);

#define _sqlite_exec(sql) sqlite3_exec(g_workspace, sql, _empty_callback, NULL, NULL)

/* R object accessors shortcuts */
#define CHAR_ELT(str, i) CHAR(STRING_ELT(str,i))

/* SDF object accessors shortcuts */
#define SDF_INAME(sdf) CHAR(STRING_ELT(_getListElement(sdf, "iname"),0))
#define SVEC_VARNAME(sdf) CHAR(STRING_ELT(_getListElement(sdf, "varname"),0))

/* # of sql buffers */
#define NBUFS 4

/* possible var types when stored in sqlite as integer */
#define VAR_INTEGER 0
#define VAR_FACTOR  1
#define VAR_ORDERED 2


#ifndef __SQLITE_WORKSPACE__
extern sqlite3 *g_workspace;
extern char *g_sql_buf[NBUFS];
extern int g_sql_buf_sz[NBUFS];
#endif

#endif
