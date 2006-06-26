#include <stdio.h>
#include <string.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"
#include "sqlite_dataframe.h"

SEXP sdf_create_sdf(SEXP df, SEXP name, SEXP sql_create) {
    SEXP ret;
    char *rname, *iname, *biname; 
    int file_idx = 0, namelen, res;

    /* check if arg name is supplied */
    if (IS_CHARACTER(name)) {
        rname = CHAR(STRING_ELT(name,0));
        if (!_is_r_sym(rname)) { /* error */ return R_NilValue; }
        namelen = strlen(rname);
        iname = R_alloc(namelen + 9, sizeof(char)); /* <name>10000.db\0 */
        sprintf(iname, "%s.db", rname);
    } else if (name == R_NilValue) {
        rname = "data";
        namelen = 5;
        iname = R_alloc(13, sizeof(char)); /* data10000.db\0 */
        file_idx = 1;
        sprintf(iname, "%s%d.db", rname, file_idx);
    } else {
        /* how to error??? */
        return R_NilValue;
    }

    /* here, iname will contain #{iname}.db */
    do {
        if (!_file_exists(iname)) break;
        namelen = sprintf(iname, "%s%d.db", rname, ++file_idx) - 3;
    } while (file_idx < 10000);
    
    /* found a free number */
    if (file_idx < 10000) {
        /* create a sdf db file by running "attach" statement to non-existent
         * db file
         */
        int sql_len, sql_len2, sql_len3, i, j;
        sqlite3_stmt *stmt;

        iname[namelen] = 0;  /* remove ".db" */
        sql_len = sprintf(g_sql_buf[0], "attach '%s.db' as %s", iname, iname);

        res = _sqlite_exec(g_sql_buf[0]);
        if (_sqlite_error(res)) return R_NilValue; /* duplicate dbname */
            

        /* 
         * create tables for the sdf db
         */
        SEXP names = GET_NAMES(df), variable, levels;
        int ncols = GET_LENGTH(names), type, *types;
        char *col_name, *class, *factor;

        /* TODO: put constraints on table after inserting everything? */

        /* create sdf_attributes table */
        sprintf(g_sql_buf[0], "create table [%s].sdf_attributes(attr text, "
                "value text, primary key (attr))", iname);
        res = _sqlite_exec(g_sql_buf[0]);
        sprintf(g_sql_buf[0], "insert into [%s].sdf_attributes values ('name',"
                "'%s');", iname, iname);
        _sqlite_exec(g_sql_buf[0]);


        /* 
         * create the create table and insert sql scripts by looping through
         * the columns of df
         */
        types = (int *)R_alloc(ncols, sizeof(int));
        sql_len = sprintf(g_sql_buf[0], "create table [%s].sdf_data ([row name] text", iname);
        sql_len2 = sprintf(g_sql_buf[1], "insert into [%s].sdf_data values(?", iname);

        for (i = 0; i < ncols; i++) {
            col_name = CHAR(STRING_ELT(names, i));

            /* add column definition to the create table sql */
            _expand_buf(0, sql_len+strlen(col_name)+10);

            variable = _getListElement(df, col_name);
            class = CHAR(STRING_ELT(GET_CLASS(variable), 0));
            type = TYPEOF(variable);
            types[i] = type;

            sql_len += sprintf(g_sql_buf[0]+sql_len, ", [%s] %s", col_name, 
                    _get_column_type(class, type));

            /* add handler to insert table sql */
            _expand_buf(1, sql_len+5);
            strcpy(g_sql_buf[1]+sql_len2, ",?");
            sql_len2 += 2; 

            /* create separate table for factors decode */
            if (strcmp(class, "factor") == 0 || strcmp(class, "ordered") == 0){
                sprintf(g_sql_buf[2], "create table [%s].[factor %s] ("
                        "level int, label text, primary key(level), "
                        "unique(label));", iname, col_name);
                res = _sqlite_exec(g_sql_buf[2]);
                if (_sqlite_error(res)) return R_NilValue; /* dup tbl name? */

                levels = GET_LEVELS(variable);
                sql_len3 = sprintf(g_sql_buf[2], "insert into [%s].[factor %s] values(?, ?);", iname, col_name);
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], sql_len3, &stmt, NULL);
                if (_sqlite_error(res)) return R_NilValue; /* dup tbl name? */

                for (j = 0; j < GET_LENGTH(levels); j++) {
                    sqlite3_reset(stmt);
                    factor = CHAR(STRING_ELT(levels, j));
                    sqlite3_bind_int(stmt, 1, j+1);
                    sqlite3_bind_text(stmt, 2, factor, strlen(factor), SQLITE_STATIC);
                    sqlite3_step(stmt);
                }
                sqlite3_finalize(stmt);
            }
        }
        
        _expand_buf(0,sql_len+35);
        /* sql_len += sprintf(g_sql_buf[0]+sql_len, ", primary key([row name]));");*/
        sql_len += sprintf(g_sql_buf[0]+sql_len, ");");
        res = _sqlite_exec(g_sql_buf[0]);
        if (_sqlite_error(res)) return R_NilValue; /* why? */

        /*
         * add the data in df to the sdf
         */
        SEXP rownames = getAttrib(df, R_RowNamesSymbol);
        int nrows = GET_LENGTH(rownames);
        char *row_name;

        sql_len2 += sprintf(g_sql_buf[1]+sql_len2, ")");
        res = sqlite3_prepare(g_workspace, g_sql_buf[1], sql_len2, &stmt, NULL);

        for (i = 0; i < nrows; i++) {
            row_name = CHAR(STRING_ELT(rownames, i));
            sqlite3_reset(stmt);

            if (*row_name)
                sqlite3_bind_text(stmt, 1, row_name, strlen(row_name), SQLITE_STATIC);
            else
                sqlite3_bind_int(stmt, 1, i);

            for (j = 0; j < ncols; j++) {
                variable = VECTOR_ELT(df, j);
                switch(types[j]) {
                    case INTSXP : 
                        sqlite3_bind_int(stmt, j+2, INTEGER(variable)[i]);
                        break;
                    case REALSXP:
                        sqlite3_bind_double(stmt, j+2, REAL(variable)[i]);
                        break;
                    case CHARSXP:
                        col_name = CHAR(STRING_ELT(variable,i));
                        sqlite3_bind_text(stmt, j+2, col_name, strlen(col_name), SQLITE_STATIC);
                }
                /* TODO: handle NA's & NULL's */
            }

            res = sqlite3_step(stmt);
            if (res != SQLITE_DONE) { 
                sqlite3_finalize(stmt);
                Rprintf("ERROR: %s\n", sqlite3_errmsg(g_workspace));
                return R_NilValue; /* why? */
            }
        }
                
        sqlite3_finalize(stmt);

        /*
         * add the new sdf to the workspace
         */
        sprintf(g_sql_buf[0], "insert into workspace(filename, internal_name)"
                " values('%s.db','%s')", iname, iname);
        res = _sqlite_exec(g_sql_buf[0]);
        if (_sqlite_error(res)) return R_NilValue; /* why? */

        /*
         * create a new object representing the sdf
         */
        ret = _create_sdf_sexp(iname);

    } else {
        Rprintf("ERROR: cannot find a free internal name for %s", rname);
        ret = R_NilValue;
    }
        
    return ret;
}

SEXP sdf_get_names(SEXP sdf) {
    char *iname = CHAR(STRING_ELT(_getListElement(sdf, "iname"),0));
    int len = sprintf(g_sql_buf[0], "select * from [%s].sdf_data;", iname);

    sqlite3_stmt *stmt;
    int res;
   
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], len, &stmt, NULL);
    if (_sqlite_error(res)) return R_NilValue;

    SEXP ret;

    len = sqlite3_column_count(stmt)-1;
    PROTECT(ret = NEW_CHARACTER(len));

    for (int i = 0; i < len; i++) {
        SET_STRING_ELT(ret, i, mkChar(sqlite3_column_name(stmt, i+1)));
    }

    sqlite3_finalize(stmt);
    UNPROTECT(1);
    return ret;
}

SEXP sdf_get_length(SEXP sdf) {
    char *iname = CHAR(STRING_ELT(_getListElement(sdf, "iname"),0));
    int len = sprintf(g_sql_buf[0], "select * from [%s].sdf_data;", iname);

    sqlite3_stmt *stmt;
    int res;
   
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], len, &stmt, NULL);
    if (_sqlite_error(res)) return R_NilValue;

    SEXP ret;

    len = sqlite3_column_count(stmt)-1;
    PROTECT(ret = NEW_INTEGER(1));
    INTEGER(ret)[0] = len;

    sqlite3_finalize(stmt);
    UNPROTECT(1);
    return ret;
}

SEXP sdf_get_rows(SEXP sdf) {
    char *iname = CHAR(STRING_ELT(_getListElement(sdf, "iname"),0));
    char **out;
    int res, ncol, nrow;
    SEXP ret;
   
    sprintf(g_sql_buf[0], "select count(*) from [%s].sdf_data;", iname);
    res = sqlite3_get_table(g_workspace, g_sql_buf[0], &out, &nrow, &ncol, NULL);
    if (_sqlite_error(res)) return R_NilValue;


    if (nrow != 1) ret = R_NilValue;
    else {
        PROTECT(ret = NEW_INTEGER(1));
        INTEGER(ret)[0] = atoi(out[1]);
        UNPROTECT(1);
    }

    sqlite3_free_table(out);
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
