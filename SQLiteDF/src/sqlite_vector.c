#include "sqlite_dataframe.h"

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
            case VAR_INTEGER: class = mkChar("numeric"); break;
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
                SET_STRING_ELT(*ret, 0, mkChar(sqlite3_column_text(stmt, 0)));
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
                    mkChar(sqlite3_column_text(stmt, 0)));
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
    int index, idxlen, i, retlen, res;

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




