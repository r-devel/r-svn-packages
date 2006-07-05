#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"

#define __UTIL__
#include "sqlite_dataframe.h"

int _is_r_sym(char *sym) {
    int i, len = strlen(sym);

    if (isalpha(sym[0]) || sym[0] == '_') i = 1;
    else if (sym[0] == '.') { if (isdigit(sym[1])) return FALSE; i = 2; }
    else return FALSE;

    for (; i < len; i++) {
        if (!(isalnum(sym[i]) || sym[i] == '.' || sym[i] == '_')) return FALSE;
    }

    return TRUE;
}

int _file_exists(char *filename) {
    FILE *f; int ret = FALSE;
    if ((f = fopen(filename, "rb")) != NULL) { fclose(f); ret = TRUE; }
    return ret;
}

int _empty_callback(void *data, int ncols, char **rows, char **cols) {
    return 0;
}

char *_fixname(char *rname) {
    char *tmp = rname;
    while (*tmp) { if (*tmp == '.') *tmp = '_'; tmp++; }
    return rname;
}


char *_r2iname(char *rname, char *iname) {
    /* basically s/[.!@#$%^&*()-+=]/_/g */
    strcpy(iname, rname);
    char *tmp = iname;
    while ((tmp = strpbrk(tmp, ".!@#$%^&*()-+=")) == NULL) *tmp = '_';
    return iname;
}

/*
int _sql_exec_noret(char *sql, ...) {
    int len;
    _check_sql_buf();

    len = sprintf(g_sql_buf, sql);
}*/

int _check_sql_buf(int i) {
    if (i >= NBUFS) return FALSE;
    int ret = strlen(g_sql_buf[i]);
    if ((ret*1.0/g_sql_buf_sz[i]) > 0.6) {
        g_sql_buf_sz[i] *= 2;
        g_sql_buf[i] = Realloc(g_sql_buf[i], g_sql_buf_sz[i], char);
    } else ret = 0;
    return ret;
} 

int _expand_buf(int i, int size) {
    int expanded = FALSE;
    if (size >= g_sql_buf_sz[i]) {
        g_sql_buf_sz[i] *= 2;
        g_sql_buf[i] = Realloc(g_sql_buf[i], g_sql_buf_sz[i], char);
        expanded = TRUE;
    }
    return expanded;
}

int _sqlite_error(int res) {
    int ret = FALSE;
    if (res != SQLITE_OK) { 
        Rprintf("ERROR: %s\n", sqlite3_errmsg(g_workspace));
        ret = TRUE;
    }
    return ret;
}

const char *_get_column_type(const char *class, int type) {
    if (type == INTSXP) return "int";
    else if (type == REALSXP) return "double";
    else if (type == CHARSXP) return "text";
    else if (type == LGLSXP) return "bit";
    else if (strcmp(class, "factor") == 0) return "int";
    else if (strcmp(class, "ordered") == 0) return "int";
    
    return NULL;
}

const char *_get_r_class(const char *db, const char *type) {
    return NULL;
}

int __count_callback(void *data, int ncols, char **rows, char **cols) {
    *((int*) data) = atoi(rows[0]);
    return 0;
}

int _get_row_count2(const char *table) {
    sprintf(g_sql_buf[2], "select count(*) from %s", table);
   
    int ret, res;
    res = sqlite3_exec(g_workspace, g_sql_buf[2], __count_callback, &ret, NULL);
    if (res != SQLITE_OK) return -1;
    return ret;
}

/* TODO: windows version */
char *_get_full_pathname2(char *relpath) {
    char *tmp1, *tmp2, tmp3, tmp4;
    int buflen, relpathlen;

    relpathlen = strlen(relpath);
    if (relpath[0] == '/') {
        _expand_buf(2, relpathlen);
        strcpy(g_sql_buf[2], "/");
        buflen = 1;
    } else {
        while (TRUE) {
            tmp1 = getcwd(g_sql_buf[2], g_sql_buf_sz[2]);
            if (tmp1 == NULL) _expand_buf(2, g_sql_buf_sz[2]+1);
            else { 
                buflen = strlen(tmp1); 
                strcpy(g_sql_buf[2]+buflen,"/");
                buflen += 1;
                break;
            }
        }
    }

    /* we'll go along the relpath string "normalizing" relative paths */
    tmp1 = relpath;
    while (tmp1[0]) {
        tmp2 = tmp1;
        while(!(*tmp2 == '/' || *tmp2 == 0)) tmp2++;
        
        tmp3 = *tmp2;
        *tmp2 = 0; /* temporarily "end" string at that point */
        if (strcmp(tmp1, ".") == 0) {
            /* nothing to do */
        } else if (strcmp(tmp1, "..") == 0) {
            /* remove top dir */
            if (buflen > 1) {
                buflen--;
                do buflen--; 
                while (g_sql_buf[2][buflen] != '/' && buflen > 1);
                g_sql_buf[2][++buflen] = 0;
            }
        } else { 
            /* non-relative path part, append to buf */
            tmp4 = 0;
            if (tmp3 == '/') { *tmp2 = '/'; tmp4 = *(tmp2+1); *(tmp2+1) = 0; }
            _expand_buf(2, buflen+strlen(tmp1));
            strcpy(g_sql_buf[2] + buflen, tmp1);
            buflen += strlen(tmp1);
            if (tmp3 == '/') { *(tmp2+1) = tmp4; }
        }

        if (tmp3 == 0) break;
        else { *tmp2 = '/'; tmp1 = tmp2 + 1; }
    }

    return g_sql_buf[2];
}

/* based on p70 of R-exts.pdf */
SEXP _getListElement(SEXP list, char *varname) {
    SEXP ret = R_NilValue, names = GET_NAMES(list);
    int i;

    for (i = 0; i < LENGTH(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), varname) == 0) {
            ret = VECTOR_ELT(list, i);
            break;
        }
    }

    return ret;
}

SEXP _get_rownames(const char *sdf_iname) {
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[0], "select [row name] from [%s].sdf_data", sdf_iname);
    int res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);

    sqlite3_finalize(stmt);
    if (_sqlite_error(res)) return R_NilValue; 

    SEXP ret, value;
    int nprotected = 1;
    PROTECT(ret = NEW_LIST(2)); 

    /* set list names */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(value, 0, mkChar("iname"));
    SET_STRING_ELT(value, 1, mkChar("varname"));
    SET_NAMES(ret, value);

    /* set list values */
    SET_VECTOR_ELT(ret, 0, mkString(sdf_iname));
    SET_VECTOR_ELT(ret, 1, mkString("row name"));

    /* set class */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(value, 0, mkChar("sqlite.vector"));
    SET_STRING_ELT(value, 1, mkChar("character"));
    SET_CLASS(ret, value);

    UNPROTECT(nprotected);
    return ret;
}

SEXP _create_sdf_sexp(const char *iname) {
    SEXP names, class, variable, ret;
    int nprotected = 0;
    PROTECT(ret = NEW_LIST(1)); nprotected++;
    PROTECT(names = mkString("iname")); nprotected++;
    SET_NAMES(ret, names);

    PROTECT(variable = mkString(iname)); nprotected++;
    SET_VECTOR_ELT(ret, 0, variable);

    /* set class */
    PROTECT(class = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(class, 0, mkChar("sqlite.data.frame"));
    SET_STRING_ELT(class, 1, mkChar("data.frame"));
    SET_CLASS(ret, class);
    SET_ROWNAMES(ret, _get_rownames(iname));
    
    UNPROTECT(nprotected);
    return ret;
}

static void __attach_levels2(char *table, SEXP var, int len) {
    SEXP levels;
    int idx = 0, res;
    sqlite3_stmt *stmt;

    PROTECT(levels = NEW_CHARACTER(len));
    sprintf(g_sql_buf[2], "select level, label from %s order by level asc",
            table);
    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
    _sqlite_error(res);
    while (sqlite3_step(stmt) == SQLITE_ROW) { 
        SET_STRING_ELT(levels, idx, mkChar((char *)sqlite3_column_text(stmt, 1)));
        idx++;
    }
    SET_LEVELS(var, levels);
    UNPROTECT(1);
}
    
int _get_factor_levels1(const char *iname, const char *varname, SEXP var) {
    int ret = VAR_INTEGER;
    
    sprintf(g_sql_buf[1], "[%s].[factor %s]", iname, varname);
    int res = _get_row_count2(g_sql_buf[1]);
    if (res > 0) { /* res is exptected to be {-1} \union I+ */
        __attach_levels2(g_sql_buf[1], var, res);
        ret = VAR_FACTOR;
    }

    sprintf(g_sql_buf[1], "[%s].[ordered %s]", iname, varname);
    res = _get_row_count2(g_sql_buf[1]);
    if (res > 0) {
        __attach_levels2(g_sql_buf[1], var, res);
        ret = VAR_ORDERED;
    }

    return ret;
}

SEXP _shrink_vector(SEXP vec, int len) {
    int origlen = LENGTH(vec);
    SEXP ret = vec;

    if (vec == R_NilValue) return vec;
    else if (origlen > len) {
        int type = TYPEOF(vec), i;
        if (type == CHARSXP) {
            PROTECT(ret = NEW_CHARACTER(len));
            for (i = 0; i < len; i++) {
                SET_STRING_ELT(ret, i, STRING_ELT(vec, i));
            }
        } else if (type == INTSXP) {
            /* for non-strings, memcpy is more efficient but ... */
            PROTECT(ret = NEW_INTEGER(len));
            for (i = 0; i < len; i++) INTEGER(ret)[i] = INTEGER(vec)[i];
        } else if (type == REALSXP) {
            PROTECT(ret = NEW_NUMERIC(len));
            for (i = 0; i < len; i++) REAL(ret)[i] = REAL(vec)[i];
        } else if (type == LGLSXP) {
            PROTECT(ret = NEW_LOGICAL(len));
            for (i = 0; i < len; i++) LOGICAL(ret)[i] = LOGICAL(vec)[i];
        } else return ret;
        UNPROTECT(1);
    }

    return ret;
}

