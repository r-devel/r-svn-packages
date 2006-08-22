#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

#ifndef __SQLITE_DATAFRAME__
#define __SQLITE_DATAFRAME__

/* # of sql buffers */
#define NBUFS 4

#ifdef __SQLITE_WORKSPACE__
SEXP SDF_RowNamesSymbol;
SEXP SDF_VectorTypeSymbol;
SEXP SDF_DimSymbol;
SEXP SDF_DimNamesSymbol;
#else 
extern SEXP SDF_RowNamesSymbol;
extern SEXP SDF_VectorTypeSymbol;
extern SEXP SDF_DimSymbol;
extern SEXP SDF_DimNamesSymbol;

extern sqlite3 *g_workspace;
extern char *g_sql_buf[NBUFS];
extern int g_sql_buf_sz[NBUFS];
#endif

/* sdf type attributes */
#define GET_SDFVECTORTYPE(x) getAttrib(x, SDF_VectorTypeSymbol) 
#define SET_SDFVECTORTYPE(x,t) setAttrib(x, SDF_VectorTypeSymbol, t) 
#define TEST_SDFVECTORTYPE(x, t) (strcmp(CHAR(asChar(GET_SDFVECTORTYPE(x))), t) == 0)
#define GET_SDFROWNAMES(x) getAttrib(x, SDF_RowNamesSymbol)
#define SET_SDFROWNAMES(x, r) setAttrib(x, SDF_RowNamesSymbol, r)

#ifndef SET_ROWNAMES
#define SET_ROWNAMES(x, r) setAttrib(x, R_RowNamesSymbol, r)
#endif

#define WORKSPACE_COLUMNS 6
#define MAX_ATTACHED 30     /* 31 including workspace.db */

/* utilities for checking characteristics of arg */
int _is_r_sym(char *sym);
int _file_exists(char *filename);
int _sdf_exists2(char *iname);

/* sdf utilities */
int USE_SDF1(const char *iname, int exists, int protect);  /* call this before doing anything on an SDF */
int UNUSE_SDF2(const char *iname); /* somewhat like UNPROTECT */
SEXP _create_sdf_sexp(const char *iname);  /* create a SEXP for an SDF */
int _add_sdf1(char *filename, char *internal_name); /* add SDF to workspace */
void _delete_sdf2(const char *iname); /* remove SDF from workspace */
int _get_factor_levels1(const char *iname, const char *varname, SEXP var);
int _get_row_count2(const char *table, int quote);
SEXP _get_rownames(const char *sdf_iname);
char *_get_full_pathname2(char *relpath); /* get full path given relpath, used in workspace mgmt */
int _is_factor2(const char *iname, const char *factor_type, const char *colname);
SEXP _get_rownames2(const char *sdf_iname);

/* utilities for creating SDF's */
char *_create_sdf_skeleton1(SEXP name, int *o_namelen, int protect);
int _copy_factor_levels2(const char *factor_type, const char *iname_src,
        const char *colname_src, const char *iname_dst, const char *colname_dst);
int _create_factor_table2(const char *iname, const char *factor_type, 
        const char *colname);
char *_create_svector1(SEXP name, const char *type, int * _namelen, int protect);

/* R utilities */
SEXP _getListElement(SEXP list, char *varname);
SEXP _shrink_vector(SEXP vec, int len); /* shrink vector size */

/* sqlite utilities */
int _empty_callback(void *data, int ncols, char **row, char **cols);
int _sqlite_error_check(int res, const char *file, int line);
const char *_get_column_type(const char *class, int type); /* get sqlite type corresponding to R class & type */
sqlite3* _is_sqlitedb(char *filename);

/* global buffer (g_sql_buf) utilities */
int _expand_buf(int i, int size);  /* expand ith buf if size > buf[i].size */


/* workspace utilities */
int _prepare_attach2();  /* prepare workspace before attaching a sqlite db */

/* misc utilities */
char *_r2iname(char *internal_name, char *filename);
char *_fixname(char *rname);
char *_str_tolower(char *out, const char *ref);

/* register functions to sqlite */
void __register_vector_math();

#define _sqlite_exec(sql) sqlite3_exec(g_workspace, sql, _empty_callback, NULL, NULL)
#define _sqlite_error(res) _sqlite_error_check((res), __FILE__, __LINE__)

#ifdef __SQLITE_DEBUG__
#define _sqlite_begin  { _sqlite_error(_sqlite_exec("begin")); Rprintf("begin at "  __FILE__  " line %d\n",  __LINE__); }
#define _sqlite_commit  { _sqlite_error(_sqlite_exec("commit")); Rprintf("commit at "  __FILE__  " line %d\n",  __LINE__); }
#else
#define _sqlite_begin  _sqlite_error(_sqlite_exec("begin")) 
#define _sqlite_commit _sqlite_error(_sqlite_exec("commit"))
#endif

/* R object accessors shortcuts */
#define CHAR_ELT(str, i) CHAR(STRING_ELT(str,i))

/* SDF object accessors shortcuts */
#define SDF_INAME(sdf) CHAR(STRING_ELT(_getListElement(sdf, "iname"),0))
#define SVEC_VARNAME(sdf) CHAR(STRING_ELT(_getListElement(sdf, "varname"),0))

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



/* top level functions */
/* sqlite_workspace.c */
SEXP sdf_init_workspace();
SEXP sdf_finalize_workspace();
SEXP sdf_list_sdfs(SEXP pattern);
SEXP sdf_get_sdf(SEXP name);
SEXP sdf_attach_sdf(SEXP filename, SEXP internal_name);
SEXP sdf_detach_sdf(SEXP internal_name);
SEXP sdf_rename_sdf(SEXP sdf, SEXP name);

/* sqlite_dataframe.c */
SEXP sdf_create_sdf(SEXP df, SEXP name);
SEXP sdf_get_names(SEXP sdf);
SEXP sdf_get_length(SEXP sdf);
SEXP sdf_get_row_count(SEXP sdf);
SEXP sdf_import_table(SEXP _filename, SEXP _name, SEXP _sep, SEXP _quote, 
        SEXP _rownames, SEXP _colnames);
SEXP sdf_get_index(SEXP sdf, SEXP row, SEXP col);
SEXP sdf_rbind(SEXP sdf, SEXP data);
SEXP sdf_get_iname(SEXP sdf);

/* sqlite_vector.c */
SEXP sdf_get_variable(SEXP sdf, SEXP name);
SEXP sdf_get_variable_length(SEXP svec);
SEXP sdf_get_variable_index(SEXP svec, SEXP idx);
/* SEXP sdf_set_variable_index(SEXP svec, SEXP idx, SEXP value); */
SEXP sdf_variable_summary(SEXP svec, SEXP maxsum);
SEXP sdf_do_variable_math(SEXP func, SEXP vector, SEXP other_args);
SEXP sdf_do_variable_op(SEXP func, SEXP vector, SEXP op2, SEXP arg_reversed);
SEXP sdf_do_variable_summary(SEXP func, SEXP vector, SEXP na_rm);
SEXP sdf_sort_variable(SEXP svec, SEXP decreasing);

/* sqlite_external.c */
SEXP sdf_import_sqlite_table(SEXP _dbfilename, SEXP _tblname, SEXP _sdfiname);

/* sqlite_matrix.c */
SEXP sdf_as_matrix(SEXP sdf, SEXP name);
#endif
