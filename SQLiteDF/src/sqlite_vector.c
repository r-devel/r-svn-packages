#include "sqlite_dataframe.h"
#include <math.h>
#include "Rmath.h"

/****************************************************************************
 * UTILITY FUNCTIONS
 ****************************************************************************/
char *_create_svector1(SEXP name, const char *type, int * _namelen, int protect) {
    int namelen, res;
    char *iname = _create_sdf_skeleton1(name, &namelen, protect);

    if (iname == NULL) return NULL;

    sprintf(g_sql_buf[2], "create table [%s].sdf_data ([row name] text, "
            "V1 %s, primary key ([row name]))", iname, type);
    res = _sqlite_exec(g_sql_buf[2]);
    _sqlite_error(res);

    if (_namelen != NULL) *_namelen = namelen;
    return iname;
}

static SEXP _create_svector_sexp(const char *iname, const char *varname, 
        const char *type) {
    SEXP ret, value; int nprotected = 0;
    PROTECT(ret = NEW_LIST(2)); nprotected++;

    /* set list names */
    PROTECT(value = NEW_CHARACTER(2)); nprotected++;
    SET_STRING_ELT(value, 0, mkChar("iname"));
    SET_STRING_ELT(value, 1, mkChar("varname"));
    SET_NAMES(ret, value);

    /* set list values */
    SET_VECTOR_ELT(ret, 0, mkString(iname));
    SET_VECTOR_ELT(ret, 1, mkString(varname));

    /* set sexp class */
    SET_CLASS(ret, mkString("sqlite.vector"));

    /* set sdf.vector.type */
    SET_SDFVECTORTYPE(ret, mkString(type));

    UNPROTECT(nprotected);

    return ret;
}

/* if ret == NULL, 3rd arg is the length of the vector created. otherwise
 * it is the index in the vector ret where we will put the result extracted
 * from ret */
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

/****************************************************************************
 * SVEC FUNCTIONS
 ****************************************************************************/
SEXP sdf_get_variable(SEXP sdf, SEXP name) {
    char *iname, *varname, *svec_type = NULL;
    const char *coltype;
    int type = -1, res, nprotected = 0;
    SEXP ret, value;

    if (!IS_CHARACTER(name)) {
        Rprintf("ERROR: argument is not a string.\n");
        return R_NilValue;
    }

    iname = SDF_INAME(sdf);
    varname = CHAR_ELT(name, 0);

    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

    /* check if sdf & varname w/in that sdf exists */
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[0], "select [%s] from [%s].sdf_data", varname, iname);

    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);

    if (_sqlite_error(res)) return R_NilValue;

    coltype = sqlite3_column_decltype(stmt, 0);
    sqlite3_finalize(stmt);
    
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
    if (strcmp(coltype, "text") == 0) svec_type = "character";
    else if (strcmp(coltype, "double") == 0) svec_type = "numeric";
    else if (strcmp(coltype, "bit") == 0) svec_type = "logical";
    else if (strcmp(coltype, "integer") == 0 || strcmp(coltype, "int") == 0) {
        /* determine if int, factor or ordered */
        type = _get_factor_levels1(iname, varname, ret);
        switch(type) {
            case VAR_INTEGER: svec_type = "integer"; break;
            case VAR_FACTOR: svec_type = "factor"; break;
            case VAR_ORDERED: svec_type = "ordered";
        }

    }

    SET_CLASS(ret, mkString("sqlite.vector"));
    SET_SDFVECTORTYPE(ret, mkString(svec_type));

    UNPROTECT(nprotected);
    return ret;

}


SEXP sdf_get_variable_length(SEXP svec) {
    char *iname = SDF_INAME(svec);
    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

    return ScalarInteger(_get_row_count2(iname, 1));
}

    
SEXP sdf_get_variable_index(SEXP svec, SEXP idx) {
    SEXP ret = R_NilValue, tmp;
    char *iname = SDF_INAME(svec), *varname = SVEC_VARNAME(svec);
    int index, idxlen, i, retlen=0, res;
    sqlite3_stmt *stmt;

    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

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
            sqlite3_reset(stmt);
            sqlite3_bind_int(stmt, 1, 0);
            res = sqlite3_step(stmt);
            if (res == SQLITE_ROW) {
                _get_vector_index_typed_result(stmt, &ret, idxlen - 1);
            }
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
            sqlite3_reset(stmt);
            sqlite3_bind_int(stmt, 1, 0);
            res = sqlite3_step(stmt);
            if (res == SQLITE_ROW) {
                _get_vector_index_typed_result(stmt, &ret, idxlen - 1);
            }
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
        int veclen = _get_row_count2(iname, 1);

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
            if (TEST_SDFVECTORTYPE(svec, "factor")) {
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

/* sqlite.vector.[<- */
SEXP sdf_set_variable_index(SEXP svec, SEXP idx, SEXP value) {
    int idx_len, val_len;

    /* match levels if factor or ordered */
    if (inherits(value, "ordered")) {
    } else if (inherits(value, "factor")) {
    }

    idx_len = LENGTH(idx);
    val_len = LENGTH(value);

    if (IS_NUMERIC(idx)) {
    } else if (IS_INTEGER(idx)) {
    } else if (IS_LOGICAL(idx)) {
    }

    return R_NilValue;
}

SEXP sdf_variable_summary(SEXP svec, SEXP maxsum) {
    char *iname, *varname, *type;
    sqlite3_stmt *stmt;
    int nprotected = 0;
    SEXP ret, names;

    iname = SDF_INAME(svec);
    varname = SVEC_VARNAME(svec);
    USE_SDF1(iname, TRUE, FALSE);

    if ((TEST_SDFVECTORTYPE(svec, "ordered") && ((type = "ordered"))) ||
            (TEST_SDFVECTORTYPE(svec, "factor") && ((type = "factor")))) {
        int nrows, i, max_rows = INTEGER(maxsum)[0];


        sprintf(g_sql_buf[0], "[%s].[%s %s]", iname, type, varname);
        nrows = _get_row_count2(g_sql_buf[0], FALSE);
        if (nrows <= max_rows) max_rows = nrows;
        else { nrows = max_rows; max_rows--; }

        PROTECT(ret = NEW_INTEGER(nrows)); nprotected = 1;
        PROTECT(names = NEW_CHARACTER(nrows)); nprotected++;

        sprintf(g_sql_buf[0], "select [%s].[%s %s].label, count(*) from "
                "[%s].sdf_data join [%s].[%s %s] on [%s].sdf_data.[%s]=[%s].[%s %s].level "
                "group by [%s].sdf_data.[%s], [%s].[%s %s].level order by count(*) desc",
                iname, type, varname, iname, iname, type, varname, iname, varname,
                iname, type, varname, iname, varname, iname, type, varname);
        sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);

        for (i = 0; i < max_rows; i++) {
            sqlite3_step(stmt);
            SET_STRING_ELT(names, i, mkChar((char *)sqlite3_column_text(stmt, 0)));
            INTEGER(ret)[i] = sqlite3_column_int(stmt, 1);
        }

        if (nrows > max_rows) {
            int others_sum = 0;
            SET_STRING_ELT(names, nrows-1, mkChar("(Others)"));
            while (sqlite3_step(stmt) == SQLITE_ROW) {
                others_sum += sqlite3_column_int(stmt, 1);
            }
            INTEGER(ret)[nrows-1] = others_sum;
        }
    } else if (TEST_SDFVECTORTYPE(svec, "logical")) {
        sprintf(g_sql_buf[0], "select count(*) from "
                "[%s].sdf_data group by [%s] order by [%s]", iname, varname, varname);

        PROTECT(names = NEW_CHARACTER(3)); nprotected = 1;
        PROTECT(ret = NEW_CHARACTER(3)); nprotected++;

        SET_STRING_ELT(names, 0, mkChar("Mode"));
        SET_STRING_ELT(names, 1, mkChar("FALSE"));
        SET_STRING_ELT(names, 2, mkChar("TRUE"));

        sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
        SET_STRING_ELT(ret, 0, mkChar("logical"));
        sqlite3_step(stmt);
        SET_STRING_ELT(ret, 1, mkChar((const char *)sqlite3_column_text(stmt, 0)));
        sqlite3_step(stmt);
        SET_STRING_ELT(ret, 2, mkChar((const char *)sqlite3_column_text(stmt, 0)));

    } else return R_NilValue;

    sqlite3_finalize(stmt);
    SET_NAMES(ret, names);
    SET_CLASS(ret, mkString("table"));
    UNPROTECT(nprotected);
    return ret;
}
            

/* the global accumulator should be safe if we only do 1 cummulative or
 * aggregate at any time. it won't work for stuffs like "select max(col)-min(col)
 * from sdf_data". */
static long double g_accumulator = 0.0; /* accumulator var for cumsum, cumprod, etc. */
static int g_start = 0;            /* flag for start of cummulation */
static int g_narm = 0;             /* for Summary group */

SEXP sdf_do_variable_math(SEXP func, SEXP vector, SEXP other_args) {
    char *iname, *iname_src, *varname_src, *funcname;
    int namelen, res;
    sqlite3_stmt *stmt;

    /* get data from arguments (function name and sqlite.vector stuffs) */
    funcname = CHAR_ELT(func, 0);
    iname_src = SDF_INAME(vector);
    varname_src = SVEC_VARNAME(vector);

    if (!USE_SDF1(iname_src, TRUE, TRUE)) return R_NilValue;

    /* create a new sdf, with 1 column named V1 */
    iname = _create_svector1(R_NilValue, "double", &namelen, TRUE);

    /* insert into <newsdf>.col, row.names select func(col), rownames */
    if (strcmp(funcname, "round") == 0 || strcmp(funcname, "signif") == 0) {
        double digits = REAL(_getListElement(other_args, "digits"))[0];
        sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
                "select [row name], %s([%s],?) from [%s].sdf_data", iname, funcname,
                varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0); 
        if (_sqlite_error(res)) goto vecmath_prepare_error;
        res = sqlite3_bind_double(stmt, 1, digits);
    } else if (strcmp(funcname, "log") == 0) {
        double base = REAL(_getListElement(other_args, "base"))[0];
        sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
                "select [row name], %s([%s],?) from [%s].sdf_data", iname, funcname,
                varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0); 
        if (_sqlite_error(res)) goto vecmath_prepare_error;
        res = sqlite3_bind_double(stmt, 1, base);
    } else {
        sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
                "select [row name], %s([%s]) from [%s].sdf_data", iname, funcname,
                varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0); 
    }

    if (_sqlite_error(res)) {
vecmath_prepare_error:
        sprintf(g_sql_buf[0], "detach %s", iname);
        _sqlite_exec(g_sql_buf[0]);

        /* we will return a string with the file name, and do file.remove
         * at R */
        iname[namelen] = '.';
        return mkString(iname);
    }

    g_accumulator = 0.0;   /* initialize accumulator */
    g_start = 1;           /* flag that we are at start of accumulating */
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    UNUSE_SDF2(iname);
    UNUSE_SDF2(iname_src);

    return _create_svector_sexp(iname, "V1", "numeric");
}

SEXP sdf_do_variable_op(SEXP func, SEXP vector, SEXP op2, SEXP arg_reversed) {
    char *iname = NULL, *iname_src, *varname_src, *funcname;
    int res, functype = -1, op2_len, svec_len, i, reversed;
    sqlite3_stmt *stmt, *stmt2;

    /* get data from arguments (function name and sqlite.vector stuffs) */
    funcname = CHAR_ELT(func, 0);
    iname_src = SDF_INAME(vector);
    varname_src = SVEC_VARNAME(vector);
    reversed = LOGICAL(arg_reversed)[0];

    if (!USE_SDF1(iname_src, TRUE, TRUE)) return R_NilValue;
    switch(funcname[0]) {
        case '+' :
        case '-' :
        case '*' :
        case '/' :
        case '^' :
        case '%' : /* %% and %/% */ 
            functype = 0; break;  /* output is REAL */
        case '&' :
        case '|' :
        case '!' :  /* ! and != */
            if (funcname[1] == 0) { functype = 1; break; }
        case '=' :  /* == */
        case '<' :  /* < and <= */
        case '>' :  /* > and >= */
            functype = 2; 
    }

    svec_len = _get_row_count2(iname_src, 1);
    if (functype == 0 || functype == 2 || (functype == 1 && funcname[0] != '!')) {
        char *insert_fmt_string1c, *insert_fmt_string1s;
        char *insert_fmt_string2c, *insert_fmt_string2s;
        int vec_idx, sdf_idx;

        iname = _create_svector1(R_NilValue, (functype == 0) ? "double" : "bit", NULL, TRUE);

        /* insert_fmt_string prefixes:
         * 1 - op2 is an ordinary vector, op2_len = 1, straightforward insert-select
         * 2 - op2 is an ordinary vector, op2_len > 1, may need recycling, loop both sides
         * s - operator will be printed as string
         * c - operator is either INTDIV or RMOD (int division, R modulo). initially, I
         *     thought I could cast both op to int then use / and % of sqlite, which is just
         *     the 2nd character of the funcname. however, R's INTDIV and casts to int 
         *     after division (e.g. 3.5 %/% 1.5 == 2). RMOD operates on doubles too,
         *     which is the remainder of the largest int multiple of "divisor"
         *     (e.g. 3.5 %% 1.5 = 0.5)
         */ 
        if (functype == 1) { /* boolean binary operators */
            if (!reversed) {
                insert_fmt_string1s = "insert into [%s].sdf_data "
                            "select [row name], ([%s] != 0) %s (? != 0) from [%s].sdf_data";
            } else {
                insert_fmt_string1s = "insert into [%s].sdf_data "
                            "select [row name], (? != 0) %s ([%s] != 0) from [%s].sdf_data";
            }
            insert_fmt_string2s = "insert into [%s].sdf_data values(?, (? != 0) %s (? != 0))";
            /* no need for char version of fmt_string2, since %% only occurs for functype==0 */
            insert_fmt_string1c = insert_fmt_string1s;
            insert_fmt_string2c = insert_fmt_string2s;
        } else {
            if (!reversed) {
                insert_fmt_string1s = "insert into [%s].sdf_data "
                            "select [row name], [%s] %s ? from [%s].sdf_data";
                insert_fmt_string1c = "insert into [%s].sdf_data "
                            "select [row name], %s([%s],?) from [%s].sdf_data";
            } else {
                insert_fmt_string1s = "insert into [%s].sdf_data "
                            "select [row name], ? %s [%s] from [%s].sdf_data";
                insert_fmt_string1c = "insert into [%s].sdf_data "
                            "select [row name], %s(?,[%s]) from [%s].sdf_data";
            }
            /* fmt_string2c format operator from a char, used for %% and %/% */
            insert_fmt_string2c = "insert into [%s].sdf_data values(?, %s(?,?))";
            insert_fmt_string2s = "insert into [%s].sdf_data values(?, ? %s ?)";
        }

        /* for 2* insert statements, the values below specifies the index of bind().
         * if !reversed, then e1 is svec, e2 is anything. otherwise, e2 is the svec */
        if (reversed) { sdf_idx = 3; vec_idx = 2; }
        else { sdf_idx = 2; vec_idx = 3; }

        if (IS_NUMERIC(op2)) {
            op2_len = LENGTH(op2);

            if (op2_len == 1) {
                if (funcname[0] == '%') {
                    sprintf(g_sql_buf[2], insert_fmt_string1c, iname, 
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv", varname_src, iname_src);
                } else {
                    sprintf(g_sql_buf[2], insert_fmt_string1s, iname,
                            varname_src, funcname, iname_src);
                }
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                _sqlite_error(res);
                sqlite3_bind_double(stmt, 1, REAL(op2)[0]);
                sqlite3_step(stmt);
                sqlite3_finalize(stmt);
            } else if (op2_len <= svec_len) {  /* recycle op2 */
                sprintf(g_sql_buf[2], "select [row name], [%s] from [%s].sdf_data",
                        varname_src, iname_src);
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, 0);
                _sqlite_error(res);

                if (funcname[0] == '%') {
                    sprintf(g_sql_buf[2], insert_fmt_string2c, iname, 
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv");  
                    _sqlite_begin;
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);

                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2);

                        sqlite3_reset(stmt);
                        sqlite3_bind_text(stmt, 1, (char *)sqlite3_column_text(stmt2, 0), -1, SQLITE_STATIC);
                        sqlite3_bind_int(stmt, sdf_idx, (int)sqlite3_column_double(stmt2, 1));
                        sqlite3_bind_int(stmt, vec_idx, (int)REAL(op2)[i % op2_len]);
                        sqlite3_step(stmt);
                    }
                    sqlite3_finalize(stmt);
                    sqlite3_finalize(stmt2);
                    _sqlite_commit;
                } else { /* non-integer binary operation */
                    sprintf(g_sql_buf[2], insert_fmt_string2s, iname, funcname);
                    _sqlite_begin;
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);

                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2);

                        sqlite3_reset(stmt);
                        sqlite3_bind_text(stmt, 1, (char *)sqlite3_column_text(stmt2, 0), -1, SQLITE_STATIC);
                        sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 1));
                        sqlite3_bind_double(stmt, vec_idx, REAL(op2)[i % op2_len]);
                        sqlite3_step(stmt);
                    }
                    sqlite3_finalize(stmt);
                    sqlite3_finalize(stmt2);
                    _sqlite_commit;
                }
            } else { /* op2_len > svec_len, recycle svec */
                sprintf(g_sql_buf[1], "select [%s] from [%s].sdf_data",
                        varname_src, iname_src);

                if (funcname[0] == '%') {
                    sprintf(g_sql_buf[2], insert_fmt_string2c, iname, 
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv");
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);
                } else {
                    sprintf(g_sql_buf[2], insert_fmt_string2s, iname, funcname);
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);
                }

                i = 0; _sqlite_begin;
                while (i < op2_len) {
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, 0);
                    _sqlite_error(res);

                    if (funcname[0] == '%') {
                        for ( ; i < op2_len && sqlite3_step(stmt) == SQLITE_ROW ; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            /* simplify my life, just use ints for row names so that we
                             * don't have to worry about duplicates */
                            sqlite3_bind_int(stmt, 1, i); 
                            sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 1));
                            sqlite3_bind_int(stmt, vec_idx, (int)REAL(op2)[i % op2_len]);
                            sqlite3_step(stmt);
                        }
                    } else {
                        for ( ; i < op2_len && sqlite3_step(stmt) == SQLITE_ROW ; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 1));
                            sqlite3_bind_double(stmt, vec_idx, REAL(op2)[i % op2_len]);
                            sqlite3_step(stmt);
                        }
                    }

                    /* recycle on svec if we loop again */
                    sqlite3_finalize(stmt2);
                }

                sqlite3_finalize(stmt);
                _sqlite_commit;
            }

        } else if (IS_INTEGER(op2) || IS_LOGICAL(op2)) { /* I N T E G E R */
            /* we are taking advantage of the fact that logicals are stored as int.
             * if this becomes untrue in the future, then this is a bug */
            op2_len = LENGTH(op2);

            if (op2_len == 1) {
                if (funcname[0] == '%') {
                    sprintf(g_sql_buf[2], insert_fmt_string1c, iname, 
                            varname_src, (funcname[1] == '%') ? "r_mod" : "r_intdiv", iname_src);
                } else {
                    sprintf(g_sql_buf[2], insert_fmt_string1s, iname,
                            varname_src, funcname, iname_src);
                }
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                _sqlite_error(res);
                sqlite3_bind_double(stmt, 1, (double)INTEGER(op2)[0]);
                sqlite3_step(stmt);
            } else if (op2_len <= svec_len) {
                sprintf(g_sql_buf[2], "select [row name], [%s] from [%s].sdf_data",
                        varname_src, iname_src);
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, 0);
                _sqlite_error(res);

                if (funcname[0] == '%') {
                    /* sqlite does int div with "/" if both args are int. mod is % also.
                     * fortunately, in both cases the sqlite op is 2nd char of funcname */
                    sprintf(g_sql_buf[2], insert_fmt_string2c, iname, 
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv");
                    _sqlite_begin;
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);

                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2);

                        sqlite3_reset(stmt);
                        sqlite3_bind_text(stmt, 1, (char *)sqlite3_column_text(stmt2, 0), -1, SQLITE_STATIC);
                        sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 1));
                        sqlite3_bind_int(stmt, vec_idx, INTEGER(op2)[i % op2_len]);
                        sqlite3_step(stmt);
                    }
                    _sqlite_commit;
                } else {
                    sprintf(g_sql_buf[2], insert_fmt_string2s, iname, funcname);
                    _sqlite_begin;
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);

                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2);

                        sqlite3_reset(stmt);
                        sqlite3_bind_text(stmt, 1, (char *)sqlite3_column_text(stmt2, 0), -1, SQLITE_STATIC);
                        sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 1));
                        sqlite3_bind_double(stmt, vec_idx, (double)INTEGER(op2)[i % op2_len]);
                        sqlite3_step(stmt);
                    }
                    _sqlite_commit;
                }
                sqlite3_finalize(stmt2);
            } else {
                sprintf(g_sql_buf[1], "select [%s] from [%s].sdf_data",
                        varname_src, iname_src);

                if (funcname[0] == '%') {
                    sprintf(g_sql_buf[2], insert_fmt_string2c, iname,
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv");
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);
                } else {
                    sprintf(g_sql_buf[2], insert_fmt_string2s, iname, funcname);
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                    _sqlite_error(res);
                }
               
                i = 0; _sqlite_begin; 
                while (i < op2_len) {
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, 0);
                    _sqlite_error(res);

                    if (funcname[0] == '%') {
                        for ( ; i < op2_len && sqlite3_step(stmt) == SQLITE_ROW ; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            /* simplify my life, just use ints for row names so that we
                             * don't have to worry about duplicates */
                            sqlite3_bind_int(stmt, 1, i); 
                            sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 1));
                            sqlite3_bind_int(stmt, vec_idx, INTEGER(op2)[i % op2_len]);
                            sqlite3_step(stmt);
                        }
                    } else {
                        for ( ; i < op2_len && sqlite3_step(stmt) == SQLITE_ROW ; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 1));
                            sqlite3_bind_double(stmt, vec_idx, (double)INTEGER(op2)[i % op2_len]);
                            sqlite3_step(stmt);
                        }
                    }

                    /* recycle on svec if we loop again */
                    sqlite3_finalize(stmt2);
                }
                _sqlite_commit;
            }

            sqlite3_finalize(stmt);
        } else if (inherits(op2, "sqlite.vector")) { 
            /* op2 is surely not a factor, as handled by the R wrapper */
            /* even though it is impossible for reversed to be FALSE, still use
             * sdf_idx and vec_idx so that code would be less confusing */
            char *iname_op2, *varname_op2;
            sqlite3_stmt *stmt3;
            iname_op2 = SDF_INAME(op2);
            varname_op2 = SVEC_VARNAME(op2);

            if (!USE_SDF1(iname_op2, TRUE, TRUE)) {
                /* delete created sqlite.vector */
                Rprintf("Warning: detaching created result SDF %s\n", iname);
                sdf_detach_sdf(mkString(iname));
                return R_NilValue;
            }
            _sqlite_begin;

            op2_len = _get_row_count2(iname_op2, 1);
            
            sprintf(g_sql_buf[2], "select [%s] from [%s].sdf_data", varname_src, iname_src);
            res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, 0);
            _sqlite_error(res);

            sprintf(g_sql_buf[2], "select [%s] from [%s].sdf_data", varname_op2, iname_op2);
            res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt3, 0);
            _sqlite_error(res);

            if (funcname[0] == '%') {
                sprintf(g_sql_buf[2], insert_fmt_string2c, iname,
                            (funcname[1] == '%') ? "r_mod" : "r_intdiv");
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                _sqlite_error(res);

                if (svec_len == op2_len) {
                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2); sqlite3_step(stmt3);

                        sqlite3_reset(stmt);
                        sqlite3_bind_int(stmt, 1, i);
                        sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 0));
                        sqlite3_bind_int(stmt, vec_idx, sqlite3_column_int(stmt3, 0));
                        sqlite3_step(stmt);
                    }
                } else if (svec_len < op2_len) {
                    i = 0;
                    while (i < op2_len) {
                        for ( ; sqlite3_step(stmt2) == SQLITE_ROW && i < op2_len; i++) {
                            sqlite3_step(stmt3);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 0));
                            sqlite3_bind_int(stmt, vec_idx, sqlite3_column_int(stmt3, 0));
                            sqlite3_step(stmt);
                        }
                        sqlite3_reset(stmt2);
                    }
                } else {
                    i = 0;
                    while (i < svec_len) {
                        for ( ; sqlite3_step(stmt3) == SQLITE_ROW && i < svec_len; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_int(stmt, sdf_idx, sqlite3_column_int(stmt2, 0));
                            sqlite3_bind_int(stmt, vec_idx, sqlite3_column_int(stmt3, 0));
                            sqlite3_step(stmt);
                        }
                        sqlite3_reset(stmt3);
                    }
                }
                _sqlite_commit;
            } else { /* not an integer op %% or %/% */
                sprintf(g_sql_buf[2], insert_fmt_string2s, iname, funcname);
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
                _sqlite_error(res);

                if (svec_len == op2_len) {
                    for (i = 0; i < svec_len; i++) {
                        sqlite3_step(stmt2); sqlite3_step(stmt3);

                        sqlite3_reset(stmt);
                        sqlite3_bind_int(stmt, 1, i);
                        sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 0));
                        sqlite3_bind_double(stmt, vec_idx, sqlite3_column_double(stmt3, 0));
                        sqlite3_step(stmt);
                    }
                } else if (svec_len < op2_len) { /* recycle svec */
                    i = 0;
                    while (i < op2_len) {
                        for ( ; sqlite3_step(stmt2) == SQLITE_ROW && i < op2_len; i++) {
                            sqlite3_step(stmt3);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 0));
                            sqlite3_bind_double(stmt, vec_idx, sqlite3_column_double(stmt3, 0));
                            sqlite3_step(stmt);
                        }
                        sqlite3_reset(stmt2);
                    }
                } else { /* svec_len > op2_len, recycle op2 */
                    i = 0;
                    while (i < svec_len) {
                        for ( ; sqlite3_step(stmt3) == SQLITE_ROW && i < svec_len; i++) {
                            sqlite3_step(stmt2);

                            sqlite3_reset(stmt);
                            sqlite3_bind_int(stmt, 1, i);
                            sqlite3_bind_double(stmt, sdf_idx, sqlite3_column_double(stmt2, 0));
                            sqlite3_bind_double(stmt, vec_idx, sqlite3_column_double(stmt3, 0));
                            sqlite3_step(stmt);
                        }
                        sqlite3_reset(stmt3);
                    }
                }
            }

            sqlite3_finalize(stmt);
            sqlite3_finalize(stmt2);
            sqlite3_finalize(stmt3);
            _sqlite_commit;

            UNUSE_SDF2(iname_op2);
        }
    } else if (functype == 1 && funcname[0] != '!') { /* unary not operator */
        iname = _create_svector1(R_NilValue, "bit", NULL, TRUE);
        sprintf(g_sql_buf[2], "insert into [%s].sdf_data "
                    "select [row name], [%s] == 0 from [%s].sdf_data", 
                    iname, varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
        _sqlite_error(res);
        sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    }

    if (iname != NULL) {
        return _create_svector_sexp(iname, "V1", 
                (functype == 0) ? "numeric" : "logical");
        UNUSE_SDF2(iname);
    }

    UNUSE_SDF2(iname_src);

    return R_NilValue;

}

SEXP sdf_do_variable_summary(SEXP func, SEXP vector, SEXP na_rm) {
    char *iname_src, *varname_src, *funcname;
    int res;
    sqlite3_stmt *stmt;
    double _ret = NA_REAL; SEXP ret;

    /* get data from arguments (function name and sqlite.vector stuffs) */
    funcname = CHAR_ELT(func, 0);
    iname_src = SDF_INAME(vector);
    varname_src = SVEC_VARNAME(vector);

    if (!USE_SDF1(iname_src, TRUE, FALSE)) return R_NilValue;

    g_narm = LOGICAL(na_rm)[0];
    if (strcmp(funcname, "range") == 0) {
        /* special handling for range. use min then max */
        g_start = 1;
        g_accumulator = 0.0;
        sprintf(g_sql_buf[0], "select min_df([%s]) from [%s].sdf_data", varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, NULL);
        if (_sqlite_error(res)) return R_NilValue;
        sqlite3_step(stmt); sqlite3_finalize(stmt);
        _ret = g_accumulator;

        if (R_IsNA(_ret) && !g_narm) {
            PROTECT(ret = NEW_NUMERIC(2));
            REAL(ret)[0] = REAL(ret)[1] = R_NaReal;
            goto __sdf_do_variable_summary_out;
        }
        
        g_start = 1;
        g_accumulator = 0.0;
        sprintf(g_sql_buf[0], "select max_df([%s]) from [%s].sdf_data", varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, NULL);
        if (_sqlite_error(res)) return R_NilValue;
        sqlite3_step(stmt); sqlite3_finalize(stmt);

        /* if there is NA, then the if above should have caught it already */
        PROTECT(ret = NEW_NUMERIC(2));
        REAL(ret)[0] = _ret;
        REAL(ret)[1] = g_accumulator;
    } else {
        g_start = 1;  /* we'll use these instead of sqlite3_aggregate_context */
        g_accumulator = 0.0;
        sprintf(g_sql_buf[0], "select %s_df([%s]) from [%s].sdf_data", funcname, varname_src, iname_src);
        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, NULL);
        if (_sqlite_error(res)) return R_NilValue;
        res = sqlite3_step(stmt); sqlite3_finalize(stmt);

        if (strcmp(funcname, "all") == 0 || strcmp(funcname, "any") == 0) {
            PROTECT(ret = NEW_LOGICAL(1));
            if (R_IsNA(g_accumulator)) LOGICAL(ret)[0] = NA_INTEGER;
            else LOGICAL(ret)[0] = !(g_accumulator == 0);
        } else {
            PROTECT(ret = NEW_NUMERIC(1));
            REAL(ret)[0] = g_accumulator;
        }
    }

__sdf_do_variable_summary_out:
    UNPROTECT(1);
    return ret;
}


SEXP sdf_sort_variable(SEXP svec, SEXP decreasing) {
    char *iname, *iname_src, *varname_src, *type;
    sqlite3_stmt *stmt;
    int res;

    iname_src = SDF_INAME(svec);
    varname_src = SVEC_VARNAME(svec);

    if (!USE_SDF1(iname_src, TRUE, TRUE)) return R_NilValue;

    /* determine type of svec */
    sprintf(g_sql_buf[0], "select [%s] from [%s].sdf_data limit 1", varname_src, iname_src);
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
    _sqlite_error(res);
    sqlite3_step(stmt);
    strcpy(g_sql_buf[0], sqlite3_column_decltype(stmt, 0));
    sqlite3_finalize(stmt);

    /* create a new vector of that type */
    iname = _create_svector1(mkString("tmp-sort"), g_sql_buf[0], NULL, TRUE);

    /* insert to new sdf ordered */
    sprintf(g_sql_buf[0], "insert into [%s].sdf_data "
            "select [row name], [%s] from [%s].sdf_data "
            "order by [%s] %s", iname, varname_src, iname_src, varname_src,
            (LOGICAL(decreasing)[0]) ? "desc" : "asc");
    res = _sqlite_exec(g_sql_buf[0]);
    _sqlite_error(res);

    if (TEST_SDFVECTORTYPE(svec, "factor")) { /* copy factor table to iname */
        if (TEST_SDFVECTORTYPE(svec, "ordered")) type = "ordered";
        else type = "factor";
        _copy_factor_levels2(type, iname_src, varname_src, iname, "V1");
    } else type = CHAR_ELT(GET_CLASS(svec), 1);

    UNUSE_SDF2(iname_src);
    UNUSE_SDF2(iname);

    return _create_svector_sexp(iname, "V1", type);
}

/****************************************************************************
 * VECTOR MATH/OPS/GROUP OPERATIONS
 ****************************************************************************/

R_INLINE int __vecmath_checkarg(sqlite3_context *ctx, sqlite3_value *arg, double *value) {
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

#define SQLITE_MATH_FUNC2(name, func) static void __vecmath_ ## name(\
        sqlite3_context *ctx, int argc, sqlite3_value **argv) { \
    double value1, value2; \
    if (__vecmath_checkarg(ctx, argv[0], &value1) && \
        __vecmath_checkarg(ctx, argv[1], &value2)) { \
        sqlite3_result_double(ctx, func((long double)value1, (long double)value2)); \
    }  \
}

#define SQLITE_MATH_FUNC_CUM(name, func) static void __vecmath_ ## name(\
        sqlite3_context *ctx, int argc, sqlite3_value **argv) { \
    double value; \
    if (__vecmath_checkarg(ctx, argv[0], &value)) { \
        if (g_start) { g_start = 0; g_accumulator = value; } \
        else g_accumulator = func(g_accumulator, value); \
        sqlite3_result_double(ctx, g_accumulator); \
    }  \
}

#define LOGBASE(a, b) log(a)/log(b)

/* SQLITE_MATH_FUNC1(abs, abs)   in SQLite */
SQLITE_MATH_FUNC1(sign, sign)   /* in R */
SQLITE_MATH_FUNC1(sqrt, sqrt)
SQLITE_MATH_FUNC1(floor, floor)
SQLITE_MATH_FUNC1(ceiling, ceil)
SQLITE_MATH_FUNC1(trunc, ftrunc) /* in R */
/*SQLITE_MATH_FUNC2(round, fprec)  2 arg, in SQLite, but override with R's version */
SQLITE_MATH_FUNC2(signif, fround) /* 2 arg, in R */
SQLITE_MATH_FUNC1(exp, exp)
SQLITE_MATH_FUNC2(log, LOGBASE) /* 2 arg */
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

#define SUM(a, b)  (a) + (b)
#define PROD(a, b) (a) * (b)
#define MIN(a, b) ((a) <= (b)) ? (a) : (b)
#define MAX(a, b) ((a) >= (b)) ? (a) : (b)
#define ALL(a, b) (((a) == 0) || ((b) == 0)) ? 0 : 1
#define ANY(a, b) (((a) == 0) && ((b) == 0)) ? 0 : 1

SQLITE_MATH_FUNC_CUM(cumsum, SUM)
SQLITE_MATH_FUNC_CUM(cumprod, PROD)
SQLITE_MATH_FUNC_CUM(cummin, MIN)
SQLITE_MATH_FUNC_CUM(cummax, MAX)

#define SQLITE_SUMMARY_FUNC(name, func) static void __vecsummary_ ## name(\
        sqlite3_context *ctx, int argc, sqlite3_value **argv) { \
    double value; \
    if (!g_narm && R_IsNA(g_accumulator)) return; /* NA if na.rm=F & NA found */ \
    if (sqlite3_value_type(argv[0]) != SQLITE_NULL) {  \
        if (sqlite3_value_type(argv[0]) == SQLITE_INTEGER) {  \
            int tmp = sqlite3_value_int(argv[0]); \
            value = (tmp == NA_INTEGER) ? R_NaReal : tmp; \
        } else value = sqlite3_value_double(argv[0]);  \
        if (R_IsNA(value)) { \
            if (!g_narm) g_accumulator = value; return; \
        } else if (g_start) { g_start = 0; g_accumulator = value; } \
        else g_accumulator = func(g_accumulator, value); \
    } \
}

SQLITE_SUMMARY_FUNC(all_df, ALL)
SQLITE_SUMMARY_FUNC(any_df, ANY)
SQLITE_SUMMARY_FUNC(sum_df, SUM)
SQLITE_SUMMARY_FUNC(prod_df, PROD)
SQLITE_SUMMARY_FUNC(min_df, MIN)
SQLITE_SUMMARY_FUNC(max_df, MAX)

static void __vecsummary_finalize(sqlite3_context *ctx) {
    /* g_accumulator already summarizes it. just return that */
    sqlite3_result_double(ctx, g_accumulator);
}

static void __r_modulo(sqlite3_context *ctx, int argc, sqlite3_value **argv) {
    if (sqlite3_value_type(argv[0]) == SQLITE_NULL ||
        sqlite3_value_type(argv[1]) == SQLITE_NULL) { 
        sqlite3_result_null(ctx); 
    } else {
        double v1, v2, q, tmp;
        if (sqlite3_value_type(argv[0]) == SQLITE_INTEGER && 
            sqlite3_value_type(argv[1]) == SQLITE_INTEGER) {
            int i1, i2;
            i1 = sqlite3_value_int(argv[0]);
            i2 = sqlite3_value_int(argv[1]);
            if (i1 > 0 && i2 > 0) { 
                sqlite3_result_int(ctx, i1 % i2); return;
            }
            v1 = i1; v2 = i2;
        } else {
            v1 = sqlite3_value_double(argv[0]);
            v2 = sqlite3_value_double(argv[1]);
        }
        /* copied from myfmod() in src/main/arithmetic.c */
        q = v1 / v2; 
        if (v2 == 0) sqlite3_result_double(ctx, R_NaN);
        tmp = v1 - floor(q) * v2;
        /* checking omitted */
        q = floor(tmp/v2);
        sqlite3_result_double(ctx, tmp - q*v2);

    }
}

static void __r_intdiv(sqlite3_context *ctx, int argc, sqlite3_value **argv) {
    if (sqlite3_value_type(argv[0]) == SQLITE_NULL ||
        sqlite3_value_type(argv[1]) == SQLITE_NULL) { 
        sqlite3_result_null(ctx); 
    } else {
        double v1, v2;
        if (sqlite3_value_type(argv[0]) == SQLITE_INTEGER && 
            sqlite3_value_type(argv[1]) == SQLITE_INTEGER) {
            int i1, i2;
            i1 = sqlite3_value_int(argv[0]);
            i2 = sqlite3_value_int(argv[1]);
            if (i1 == NA_INTEGER || i2 == NA_INTEGER) {
                sqlite3_result_int(ctx, NA_INTEGER); return;
            } else if (i2 == 0) {
                sqlite3_result_int(ctx, 0); return;
            }
            v1 = i1; v2 = i2;
        } else {
            v1 = sqlite3_value_double(argv[0]);
            v2 = sqlite3_value_double(argv[1]);
        }
        /* copied from IDIVOP cases in src/main/arithmetic.c */
        sqlite3_result_double(ctx, floor(v1/v2));
    }
}
#define VMENTRY1(func)  {#func, __vecmath_ ## func}
#define VSENTRY1(func)  {#func, __vecsummary_ ## func}
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
        VMENTRY1(trigamma),
        VMENTRY1(cumsum),
        VMENTRY1(cumprod),
        VMENTRY1(cummin),
        VMENTRY1(cummax)
    }, arr_func2[] =  {
        {"r_mod", __r_modulo},
        {"r_intdiv", __r_intdiv},
        /*VMENTRY1(round),*/
        VMENTRY1(signif),
        VMENTRY1(log)
    }, arr_sum1[] = {
        VSENTRY1(all_df),  /* can't override sum, min, max */
        VSENTRY1(any_df),
        VSENTRY1(sum_df),
        VSENTRY1(prod_df),
        VSENTRY1(min_df),
        VSENTRY1(max_df)
    };

    int len = sizeof(arr_func1) / sizeof(arr_func1[0]);

    for (i = 0; i < len; i++) {
        res = sqlite3_create_function(g_workspace, arr_func1[i].name, 1, 
                SQLITE_ANY, NULL, arr_func1[i].func, NULL, NULL);
        _sqlite_error(res);
    }

    len = sizeof(arr_func2) / sizeof(arr_func2[0]);
    for (i = 0; i < len; i++) {
        res = sqlite3_create_function(g_workspace, arr_func2[i].name, 2, 
                SQLITE_ANY, NULL, arr_func2[i].func, NULL, NULL);
        _sqlite_error(res);
    }

    len = sizeof(arr_sum1) / sizeof(arr_sum1[0]);
    for (i = 0; i < len; i++) {
        res = sqlite3_create_function(g_workspace, arr_sum1[i].name, 1, 
                SQLITE_ANY, NULL, NULL, arr_sum1[i].func, __vecsummary_finalize);
        _sqlite_error(res);
    }
}
