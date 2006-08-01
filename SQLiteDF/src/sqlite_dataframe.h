#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

#ifndef __SQLITE_DATAFRAME__
#define __SQLITE_DATAFRAME__

#define WORKSPACE_COLUMNS 5
#ifndef WIN32
#define MAX_ATTACHED 30     /* 31 including workspace.db */
#else
#define MAX_ATTACHED 10     /* set to 10 until I have recompiled it to 30 */
#endif

/* utilities for checking characteristics of arg */
int _is_r_sym(char *sym);
int _file_exists(char *filename);
int _sdf_exists2(char *iname);

/* sdf utilities */
int USE_SDF1(const char *iname, int exists);  /* call this before doing anything on an SDF */
SEXP _create_sdf_sexp(const char *iname);  /* create a SEXP for an SDF */
int _add_sdf1(char *filename, char *internal_name); /* add SDF to workspace */
void _delete_sdf2(const char *iname); /* remove SDF from workspace */
int _get_factor_levels1(const char *iname, const char *varname, SEXP var);
int _get_row_count2(const char *table, int quote);
SEXP _get_rownames(const char *sdf_iname);
char *_get_full_pathname2(char *relpath); /* get full path given relpath, used in workspace mgmt */

/* utilities for creating SDF's */
char *_create_sdf_skeleton1(SEXP name, int *o_namelen);
int _copy_factor_levels2(const char *factor_type, const char *iname_src,
        const char *colname_src, const char *iname_dst, const char *colname_dst);

/* R utilities */
SEXP _getListElement(SEXP list, char *varname);
SEXP _shrink_vector(SEXP vec, int len); /* shrink vector size */

/* sqlite utilities */
int _empty_callback(void *data, int ncols, char **row, char **cols);
int _sqlite_error(int res);
const char *_get_column_type(const char *class, int type); /* get sqlite type corresponding to R class & type */

/* global buffer (g_sql_buf) utilities */
int _expand_buf(int i, int size);  /* expand ith buf if size > buf[i].size */


/* sqlite.vector utilities */
SEXP sdf_get_variable(SEXP sdf, SEXP name);
SEXP sdf_detach_sdf(SEXP internal_name);

/* misc utilities */
char *_r2iname(char *internal_name, char *filename);
char *_fixname(char *rname);

/* register functions to sqlite */
void __register_vector_math();

#define _sqlite_exec(sql) sqlite3_exec(g_workspace, sql, _empty_callback, NULL, NULL)
#ifndef SET_ROWNAMES
#define SET_ROWNAMES(x,n) setAttrib(x, R_RowNamesSymbol, n)
#endif

/* override R's, which does a GetRowNames which is different I believe */
#undef GET_ROWNAMES
#define GET_ROWNAMES(x) getAttrib(x, R_RowNamesSymbol)

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

/* detail constants (see _get_sdf_detail2 in sqlite_workspace.c) */
#define SDF_DETAIL_EXISTS 0
#define SDF_DETAIL_FULLFILENAME 1

/* R SXP type constants */
#define FACTORSXP 11
#define ORDEREDSXP 12

#ifndef __SQLITE_WORKSPACE__
extern sqlite3 *g_workspace;
extern char *g_sql_buf[NBUFS];
extern int g_sql_buf_sz[NBUFS];
#endif

#endif
