#include "sqlite_dataframe.h"
#include <math.h>
#include "Rmath.h"

SEXP sdf_get_variable(SEXP sdf, SEXP name) {
    if (!IS_CHARACTER(name)) {
        Rprintf("ERROR: argument is not a string.\n");
        return R_NilValue;
    }

    char *iname = SDF_INAME(sdf);
    char *varname = CHAR_ELT(name, 0);

    /* check if sdf & varname w/in that sdf exists */
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[0], "select [%s] from [%s].sdf_data", varname, iname);

    int res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);

    if (_sqlite_error(res)) return R_NilValue;

    const char *coltype = sqlite3_column_decltype(stmt, 0);
    sqlite3_finalize(stmt);


    SEXP ret, value, class = R_NilValue; int nprotected = 0;
    PROTECT(ret = NEW_LIST(2)); nprotected++;

    /* set list names */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(value, 0, mkChar("iname"));
    SET_STRING_ELT(value, 1, mkChar("varname"));
    SET_NAMES(ret, value);

    /* set list values */
    SET_VECTOR_ELT(ret, 0, mkString(iname));
    SET_VECTOR_ELT(ret, 1, mkString(varname));

    /* set class */
    int type = -1;
    if (strcmp(coltype, "text") == 0) class = mkChar("character");
    else if (strcmp(coltype, "double") == 0) class = mkChar("numeric");
    else if (strcmp(coltype, "bit") == 0) class = mkChar("logical");
    else if (strcmp(coltype, "integer") == 0 || strcmp(coltype, "int") == 0) {
        /* determine if int, factor or ordered */
        type = _get_factor_levels1(iname, varname, ret);
        switch(type) {
            case VAR_INTEGER: class = mkChar("integer"); break;
            case VAR_FACTOR: class = mkChar("factor"); break;
            case VAR_ORDERED: class = mkChar("ordered");
        }

    }

    if (type != VAR_ORDERED) {
        PROTECT(value = NEW_CHARACTER(2)); nprotected++;
        SET_STRING_ELT(value, 0, mkChar("sqlite.vector"));
        SET_STRING_ELT(value, 1, class);
    } else {
        PROTECT(value = NEW_CHARACTER(3)); nprotected++;
        SET_STRING_ELT(value, 0, mkChar("sqlite.vector"));
        SET_STRING_ELT(value, 1, class);
        SET_STRING_ELT(value, 2, mkChar("factor"));
    }
    SET_CLASS(ret, value);

    UNPROTECT(nprotected);
    return ret;

}

int _get_vector_index_typed_result(sqlite3_stmt *stmt, SEXP *ret, int idx_or_len) {
    int added = 1;
    if (*ret == NULL || *ret == R_NilValue) {
        const char *coltype = sqlite3_column_decltype(stmt, 0);
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL) {
            added = 0;
        }
        
        if (strcmp(coltype, "text") == 0) {
            PROTECT(*ret = NEW_CHARACTER(idx_or_len));
            if (added) 
                SET_STRING_ELT(*ret, 0, mkChar((char *)sqlite3_column_text(stmt, 0)));
        } else if (strcmp(coltype, "double") == 0) {
            PROTECT(*ret = NEW_NUMERIC(idx_or_len));
            if (added) REAL(*ret)[0] = sqlite3_column_double(stmt, 0);
        } else if (strcmp(coltype, "bit") == 0) {
            PROTECT(*ret = NEW_LOGICAL(idx_or_len));
            if (added) INTEGER(*ret)[0] = sqlite3_column_int(stmt, 0);
        } else if (strcmp(coltype, "integer") == 0 || 
                   strcmp(coltype, "int") == 0) {
            /* caller should just copy off the vars level attr for factors */
            PROTECT(*ret = NEW_INTEGER(idx_or_len));
            if (added) INTEGER(*ret)[0] = sqlite3_column_int(stmt, 0);
        } else added = 0;

        UNPROTECT(1);
    } else {
        const char *coltype = sqlite3_column_decltype(stmt, 0);
        if (sqlite3_column_type(stmt, 0) == SQLITE_NULL) {
            added = 0;
        } else if (strcmp(coltype, "text") == 0) {
            SET_STRING_ELT(*ret, idx_or_len, 
                    mkChar((char *)sqlite3_column_text(stmt, 0)));
        } else if (strcmp(coltype, "double") == 0) {
            REAL(*ret)[idx_or_len] = sqlite3_column_double(stmt, 0);
        } else if (strcmp(coltype, "bit") == 0) {
            INTEGER(*ret)[idx_or_len] = sqlite3_column_int(stmt, 0);
        } else if (strcmp(coltype, "integer") == 0 ||
                   strcmp(coltype, "int") == 0) {
            /* caller should just copy off the vars level attr for factors */
            INTEGER(*ret)[idx_or_len] = sqlite3_column_int(stmt, 0);
        } else added = 0;
    }
        
    return added;
}

SEXP sdf_get_variable_length(SEXP svec) {
    char *iname = SDF_INAME(svec);
    sprintf(g_sql_buf[0], "[%s].sdf_data", iname);

    SEXP ret;
    PROTECT(ret = NEW_INTEGER(1));
    INTEGER(ret)[0] = _get_row_count2(g_sql_buf[0]);
    UNPROTECT(1);
    return ret;
}

    
SEXP sdf_get_variable_index(SEXP svec, SEXP idx) {
    SEXP ret = R_NilValue, tmp;
    char *iname = SDF_INAME(svec), *varname = SVEC_VARNAME(svec);
    int index, idxlen, i, retlen=0, res;

    /* check if sdf exists */
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[0], "select [%s] from [%s].sdf_data", varname, iname);
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
    sqlite3_finalize(stmt);
    if (_sqlite_error(res)) { return ret; }

    idxlen = LENGTH(idx);
    if (idxlen < 1) return ret;

    sprintf(g_sql_buf[0], "select [%s] from [%s].sdf_data limit ?,1",
            varname, iname);
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);

    /* get data based on index */
    if (IS_NUMERIC(idx)) {
        index = ((int) REAL(idx)[0]) - 1;
        if (index < 0 && idxlen == 1) return ret;

        if (index >= 0) {
            sqlite3_bind_int(stmt, 1, index);
            res = sqlite3_step(stmt);
            if (res == SQLITE_ROW) { 
                retlen = _get_vector_index_typed_result(stmt, &ret, idxlen);
            } 
        } 
        
        if (index < 0 || res != SQLITE_ROW) {
            /* something wrong w/ 1st idx, and it is quietly ignored. we make 
             * a "dummy" call to setup the SEXP */
            sqlite3_bind_int(stmt, 1, 0);
            _get_vector_index_typed_result(stmt, &ret, idxlen - 1);
            retlen = 0;
        }

        if (idxlen > 1) {
            for (i = 1; i < idxlen; i++) {
                index = ((int) REAL(idx)[i]) - 1;
                if (index < 0) continue;
                sqlite3_reset(stmt);
                sqlite3_bind_int(stmt, 1, index);
                res = sqlite3_step(stmt);
                if (res == SQLITE_ROW)
                    retlen += _get_vector_index_typed_result(stmt, &ret, retlen); 
            }
        }
    } else if (IS_INTEGER(idx)) {
        /* similar to REAL (IS_NUMERIC) above, except that we don't have
         * to cast idx to int. can't refactor this out, sucks. */
        index = INTEGER(idx)[0] - 1;
        if (index < 0 && idxlen == 1) return ret;

        if (index >= 0) {
            sqlite3_bind_int(stmt, 1, index);
            res = sqlite3_step(stmt);
            if (res == SQLITE_ROW) { 
                retlen = _get_vector_index_typed_result(stmt, &ret, idxlen);
            } 
        } 
        
        if (index < 0 || res != SQLITE_ROW) {
            /* something wrong w/ 1st idx, and it is quietly ignored. we make 
             * a "dummy" call to setup the SEXP */
            sqlite3_bind_int(stmt, 1, 0);
            _get_vector_index_typed_result(stmt, &ret, idxlen - 1);
            retlen = 0;
        }

        if (idxlen > 1) {
            for (i = 1; i < idxlen; i++) {
                index = INTEGER(idx)[i] - 1;
                if (index < 0) continue;
                sqlite3_reset(stmt);
                sqlite3_bind_int(stmt, 1, index);
                sqlite3_step(stmt);
                retlen += _get_vector_index_typed_result(stmt, &ret, retlen); 
            }
        }

    } else if (IS_LOGICAL(idx)) {
        /* have to deal with recycling */
        sprintf(g_sql_buf[0], "[%s].sdf_data", iname);
        int veclen = _get_row_count2(g_sql_buf[0]);

        /* find if there is any TRUE element in the vector */
        for (i = 0; i < idxlen && i < veclen; i++) {
            if (LOGICAL(idx)[i]) {
                sqlite3_bind_int(stmt, 1, i);
                sqlite3_step(stmt);
                /* there are at least (idxlen-i) TRUE per cycle of the LOGICAL
                 * index. there are at least (veclen/idxlen) cycles (int div).
                 * at the last cycle, if (veclen%idxlen > 0), there will be
                 * at least (veclen%idxlen - i) if veclen%idxlen > i */
                retlen = (idxlen-i) * (veclen/idxlen);
                if (veclen%idxlen > i) retlen += (veclen%idxlen - i);

                /* create the vector */
                retlen = _get_vector_index_typed_result(stmt, &ret, retlen);
                break;
            }
        }

        if (i < idxlen && i < veclen) {
            for (i++; i < veclen; i++) {
                if (LOGICAL(idx)[i%idxlen]) {
                    sqlite3_reset(stmt);
                    sqlite3_bind_int(stmt, 1, i);
                    sqlite3_step(stmt);
                    retlen += _get_vector_index_typed_result(stmt, &ret, retlen);
                }
            }
        }
    }

    sqlite3_finalize(stmt);

    if (ret != R_NilValue) {
        ret = _shrink_vector(ret, retlen);
        tmp = GET_LEVELS(svec);
        if (tmp != R_NilValue) {
            SET_LEVELS(ret, duplicate(tmp));
            if (LENGTH(GET_CLASS(svec)) == 2) {
                SET_CLASS(ret, mkString("factor"));
            } else {
                PROTECT(tmp = NEW_CHARACTER(2));
                SET_STRING_ELT(tmp, 0, mkChar("ordered"));
                SET_STRING_ELT(tmp, 1, mkChar("factor"));
                UNPROTECT(1);
            }
        }
    }

    return ret;
}


SEXP sdf_do_variable_math(SEXP func, SEXP vector, SEXP extra_args, SEXP _nargs) {
    char *iname, *iname_src, *varname_src, *funcname;
    int namelen, res, nargs;
    sqlite3_stmt *stmt;

    /* get data from arguments (function name and sqlite.vector stuffs) */
    funcname = CHAR_ELT(func, 0);
    iname_src = SDF_INAME(vector);
    varname_src = SVEC_VARNAME(vector);
    nargs = INTEGER(_nargs)[0];

    /* check nargs */
    if (nargs > 2) {
        Rprintf("Error: Don't know how to handle Math functions w/ more than 2 args\n");
        return R_NilValue;
    }

    /* create a new sdf, with 1 column named V1 */
    iname = _create_sdf_skeleton2(R_NilValue, &namelen);
    if (iname == NULL) return R_NilValue;

    sprintf(g_sql_buf[0], "create table [%s].sdf_data ([row name] text, "
            "V1 double)", iname);
    res = _sqlite_exec(g_sql_buf[0]);
    _sqlite_error(res);

    /* insert into <newsdf>.col, row.names select func(col), rownames */
    sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
            "select [row name], %s([%s]) from [%s].sdf_data", iname, funcname,
            varname_src, iname_src);

    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0); 
    if (_sqlite_error(res)) {
        sprintf(g_sql_buf[0], "detach %s", iname);
        _sqlite_exec(g_sql_buf[0]);

        /* we will return a string with the file name, and do file.remove
         * at R */
        iname[namelen] = '.';
        return mkString(iname);
    }

    sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    /* add to workspace */
    strcpy(g_sql_buf[0], iname);
    iname[namelen] = '.';
    _add_sdf1(iname,  g_sql_buf[0]);
    iname[namelen] = 0;

    /* create sqlite.vector sexp */
    SEXP ret, value; int nprotected = 0;
    PROTECT(ret = NEW_LIST(2)); nprotected++;

    /* set list names */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(value, 0, mkChar("iname"));
    SET_STRING_ELT(value, 1, mkChar("varname"));
    SET_NAMES(ret, value);

    /* set list values */
    SET_VECTOR_ELT(ret, 0, mkString(iname));
    SET_VECTOR_ELT(ret, 1, mkString("V1"));

    /* set sexp class */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_VECTOR_ELT(value, 0, mkChar("sqlite.vector"));
    SET_VECTOR_ELT(value, 1, mkChar("numeric"));
    SET_CLASS(ret, value);

    UNPROTECT(nprotected);

    return ret;

}



/****************************************************************************
 * VECTOR MATH/OPS/GROUP OPERATIONS
 ****************************************************************************/
int __vecmath_checkarg(sqlite3_context *ctx, sqlite3_value *arg, double *value) {
    int ret = 1;
    if (sqlite3_value_type(arg) == SQLITE_NULL) { 
        sqlite3_result_null(ctx); 
        ret = 0;
    } else {
        if (sqlite3_value_type(arg) == SQLITE_INTEGER) 
            *value = sqlite3_value_int(arg); 
        else *value = sqlite3_value_double(arg); 
    }
    return ret;
}

#define SQLITE_MATH_FUNC1(name, func) static void __vecmath_ ## name(\
        sqlite3_context *ctx, int argc, sqlite3_value **argv) { \
    double value; \
    if (__vecmath_checkarg(ctx, argv[0], &value)) { \
        sqlite3_result_double(ctx, func(value)); \
    }  \
}

/* SQLITE_MATH_FUNC1(abs, abs)   in SQLite */
SQLITE_MATH_FUNC1(sign, sign)   /* in R */
SQLITE_MATH_FUNC1(sqrt, sqrt)
SQLITE_MATH_FUNC1(floor, floor)
SQLITE_MATH_FUNC1(ceiling, ceil)
SQLITE_MATH_FUNC1(trunc, ftrunc) /* in R */
/* SQLITE_MATH_FUNC1(round, )   in SQLite */
/* SQLITE_MATH_FUNC1(signif, ) 2 arg */
SQLITE_MATH_FUNC1(exp, exp)
/* SQLITE_MATH_FUNC1(log, ) 2 arg */
SQLITE_MATH_FUNC1(cos, cos)
SQLITE_MATH_FUNC1(sin, sin)
SQLITE_MATH_FUNC1(tan, tan)
SQLITE_MATH_FUNC1(acos, acos)
SQLITE_MATH_FUNC1(asin, asin)
SQLITE_MATH_FUNC1(atan, atan)
SQLITE_MATH_FUNC1(cosh, cosh)
SQLITE_MATH_FUNC1(sinh, sinh)
SQLITE_MATH_FUNC1(tanh, tanh)
SQLITE_MATH_FUNC1(acosh, acosh)  /* nowhere in include?? */
SQLITE_MATH_FUNC1(asinh, asinh)  /* nowhere in include?? */
SQLITE_MATH_FUNC1(atanh, atanh)  /* nowhere in include?? */
SQLITE_MATH_FUNC1(lgamma, lgammafn) /* in R */
SQLITE_MATH_FUNC1(gamma, gammafn) /* in R */
/* SQLITE_MATH_FUNC1(gammaCody, gammaCody)   * in R ?? */
SQLITE_MATH_FUNC1(digamma, digamma) /* in R */    
SQLITE_MATH_FUNC1(trigamma, trigamma) /* in R */

#define VMENTRY1(func)  {#func, __vecmath_ ## func}
void __register_vector_math() {
    int i, res;
    static const struct {
        char *name;
        void (*func)(sqlite3_context*, int, sqlite3_value**);
    } arr_func1[] = {
        VMENTRY1(sign),
        VMENTRY1(sqrt),
        VMENTRY1(floor),
        VMENTRY1(ceiling),
        VMENTRY1(trunc),
        VMENTRY1(exp),
        VMENTRY1(cos),
        VMENTRY1(sin),
        VMENTRY1(tan),
        VMENTRY1(acos),
        VMENTRY1(asin),
        VMENTRY1(atan),
        VMENTRY1(cosh),
        VMENTRY1(sinh),
        VMENTRY1(tanh),
        VMENTRY1(acosh),
        VMENTRY1(asinh),
        VMENTRY1(atanh),
        VMENTRY1(lgamma),
        VMENTRY1(gamma),
        VMENTRY1(digamma),
        VMENTRY1(trigamma)
    };

    int func1_len = sizeof(arr_func1) / sizeof(arr_func1[0]);

    for (i = 0; i < func1_len; i++) {
        res = sqlite3_create_function(g_workspace, arr_func1[i].name, 1, 
                SQLITE_ANY, NULL, arr_func1[i].func, NULL, NULL);
        _sqlite_error(res);
    }
}
