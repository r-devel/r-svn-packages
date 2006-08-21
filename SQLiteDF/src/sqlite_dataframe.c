#include <stdio.h>
#include <string.h>
#include "sqlite_dataframe.h"

/****************************************************************************
 * UTILITY FUNCTIONS
 ****************************************************************************/

/* if user supplied a name, return that name. otherwise, use the 
 * default "data" */
static int _check_sdf_name(SEXP name, char **rname, char **iname, int *file_idx) {
    int namelen = 0;

    /* check if arg name is supplied */
    if (name == R_NilValue) {
        *rname = "data";
        namelen = 5;
        *iname = (char*)R_alloc(13, sizeof(char)); /* data10000.db\0 */
        *file_idx = 1;
        sprintf(*iname, "%s%d.db", *rname, *file_idx);
    } else if (IS_CHARACTER(name)) {
        *rname = CHAR_ELT(name,0);
        if (!_is_r_sym(*rname)) { 
            Rprintf("Error: supplied name \"%s\"is not a valid R symbol.\n", 
                    *rname); 
            return FALSE; 
        }
        namelen = strlen(*rname);
        *iname = (char*)R_alloc(namelen + 9, sizeof(char)); /* <name>10000.db\0 */
        sprintf(*iname, "%s.db", *rname);
    } else {
        Rprintf("Error: the supplied value for arg name is not a string.\n");
    }
    return namelen;
}

static int _find_free_filename(char *rname, char **iname, int *namelen, int *file_idx) {
    do {
        if (!_file_exists(*iname)) break;
        *namelen = sprintf(*iname, "%s%d.db", rname, ++(*file_idx)) - 3;
    } while (*file_idx < 10000);

    return *file_idx;
}

/* creates an sdf attribute table */
static int _create_sdf_attribute2(const char *iname) {
    sprintf(g_sql_buf[2], "create table [%s].sdf_attributes(attr text, "
            "value text, primary key (attr))", iname);
    int res;
    res = _sqlite_exec(g_sql_buf[2]);
    if (res == SQLITE_OK) {
        sprintf(g_sql_buf[2], "insert into [%s].sdf_attributes values ('name',"
                "'%s');", iname, iname);
        res = _sqlite_exec(g_sql_buf[2]);
    }
    return res;
}

char *_create_sdf_skeleton1(SEXP name, int *onamelen, int protect) {
    char *iname, *rname;
    int namelen, file_idx = 0, res;

    namelen = _check_sdf_name(name, &rname, &iname, &file_idx);

    if (!namelen) return NULL;

    _find_free_filename(rname, &iname, &namelen, &file_idx);

    if (file_idx >= 10000) { 
        Rprintf("Error: cannot find free SDF name.\n");
        return NULL;
    }

    /* add to workspace */
    strcpy(g_sql_buf[3], iname);
    iname[namelen] = 0; /* remove ".db" */
    res = _add_sdf1(g_sql_buf[3], iname);
    if (_sqlite_error(res)) return NULL;

    /* detach SDF's if necessary to attach this one. if file does not
     * exist, then a sqlite db will be created after we make our 1st table */
    if (!USE_SDF1(iname, FALSE, protect)) { /* _delete_sdf2(iname); */ return NULL; }


    /* create attributes table */
    res = _create_sdf_attribute2(iname);
    if (_sqlite_error(res)) {
        sprintf(g_sql_buf[2], "detach '%s'", iname);
        _delete_sdf2(iname);
        return NULL; 
    }

    if (onamelen) *onamelen = namelen;
    return iname;
}
/* checks if a column has a corresponding factor|ordered table */
int _is_factor2(const char *iname, const char *factor_type, const char *colname) {
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[2], "select * from [%s].[%s %s]", iname, factor_type,
            colname);
    int res;
    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
    sqlite3_finalize(stmt);
    return res == SQLITE_OK;
}

/* create a factor|ordered table */
int _create_factor_table2(const char *iname, const char *factor_type, 
        const char *colname) {
    sqlite3_stmt *stmt;
    sprintf(g_sql_buf[2], "create table [%s].[%s %s] (level int, label text, "
            "primary key(level), unique(label));", iname, factor_type, colname);
    int res;
    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
    if (res == SQLITE_OK) sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    return res; /* error on dup name? */
}

/* copy a factor|ordered table from a sdf(db) to another sdf(db) */
int _copy_factor_levels2(const char *factor_type, const char *iname_src,
        const char *colname_src, const char *iname_dst, const char *colname_dst) {
    sqlite3_stmt *stmt;
    int res;
    res = _create_factor_table2(iname_dst, factor_type, colname_dst);
    if (res == SQLITE_OK) {
        sprintf(g_sql_buf[2], "insert into [%s].[%s %s] select * from [%s].[%s %s]",
                iname_src, factor_type, colname_src, iname_dst, factor_type,
                colname_dst);
        res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, 0);
        if (res == SQLITE_OK) sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    }
    return res; /* error on dup name? */
}


/* assumes stmt has col_cnt + 1 columns, col 0 is [row name]. returned SEXP is 
 * not UNPROTECT-ed, user will have to do that. will create the names &
 * attach factor infos */
static SEXP _setup_df_sexp1(sqlite3_stmt *stmt, const char *iname, 
        int col_cnt, int row_cnt, int *dup_indices) {
    SEXP ret, names, value, class = R_NilValue;
    int i, type;
    const char *colname, *coltype;

    PROTECT(ret = NEW_LIST(col_cnt));

    /* set up names. */
    PROTECT(names = NEW_CHARACTER(col_cnt));
    for (i = 0; i < col_cnt; i++) {
        colname = sqlite3_column_name(stmt, i+1);  /* +1 bec [row name is 0] */
        coltype = sqlite3_column_decltype(stmt, i+1);

        if (dup_indices[i] == 0) {
            strcpy(g_sql_buf[1], colname);
        } else {
            sprintf(g_sql_buf[1], "%s.%d", colname, dup_indices[i]);
        }

        SET_STRING_ELT(names, i, mkChar(g_sql_buf[1]));

        if (strcmp(coltype, "text") == 0) {
            PROTECT(value = NEW_CHARACTER(row_cnt));
        } else if (strcmp(coltype, "double") == 0) {
            PROTECT(value = NEW_NUMERIC(row_cnt));
        } else if (strcmp(coltype, "bit") == 0) {
            PROTECT(value = NEW_LOGICAL(row_cnt));
        } else if (strcmp(coltype, "integer") == 0 || 
                   strcmp(coltype, "int") == 0) {
            PROTECT(value = NEW_INTEGER(row_cnt));
            type = _get_factor_levels1(iname, colname, value);
            if (type == VAR_FACTOR) {
                PROTECT(class = mkString("factor"));
                SET_CLASS(value, class);
                UNPROTECT(1);
            } else if (type == VAR_ORDERED) {
                PROTECT(class = NEW_CHARACTER(2));
                SET_STRING_ELT(class, 0, mkChar("ordered"));
                SET_STRING_ELT(class, 1, mkChar("factor"));
                SET_CLASS(value, class);
                UNPROTECT(1);
            }
        } else {
            Rprintf("Error: not supported type %s for %s\n", coltype, colname);
            UNPROTECT(2); /* unprotect ret, names */
            return R_NilValue;
        }

        SET_VECTOR_ELT(ret, i, value);
        UNPROTECT(1); /* unprotect value */
                    
    }
    SET_NAMES(ret, names);
    UNPROTECT(1); /* unprotect names only */
    return ret;
}

/* expected that 1st col of stmt is [row name], so we start with index 1 */
static int _add_row_to_df(SEXP df, sqlite3_stmt *stmt, int row, int ncols) {
    SEXP vec;
    int type = -1, is_null, i;
    for (i = 1; i <= ncols; i++) {
        vec = VECTOR_ELT(df, i-1);
        type = TYPEOF(vec);

        is_null = (sqlite3_column_type(stmt, row) == SQLITE_NULL);

        if (type == CHARSXP) {
            SET_STRING_ELT(vec, row, mkChar((char *)sqlite3_column_text(stmt, i)));
        } else if (type == INTSXP) {
            INTEGER(vec)[row] = sqlite3_column_int(stmt, i);
        } else if (type == REALSXP) {
            REAL(vec)[row] = sqlite3_column_double(stmt, i);
        } else if (type == LGLSXP) {
            LOGICAL(vec)[row] = sqlite3_column_int(stmt, i);
        }
    }
    return type;
}

static void _set_rownames2(SEXP df) {
    SEXP value = VECTOR_ELT(df, 0);
    int len = LENGTH(value), i;
    PROTECT(value = NEW_CHARACTER(len));
    for (i = 0; i < len; i++) {
        sprintf(g_sql_buf[2], "%d", i+1);
        SET_STRING_ELT(value, i, mkChar(g_sql_buf[2]));
    }

    SET_ROWNAMES(df, value);
    UNPROTECT(1);
}


/****************************************************************************
 * SDF FUNCTIONS
 ****************************************************************************/

SEXP sdf_create_sdf(SEXP df, SEXP name) {
    SEXP ret;
    char *iname; 
    int namelen, res, i, j;

    /* find free name, attach sdf, create sdf_attributes */
    iname = _create_sdf_skeleton1(name, &namelen, FALSE);
    
    if (iname != NULL) {
        int sql_len, sql_len2;
        sqlite3_stmt *stmt;

        /* create sdf_data table */
        SEXP names = GET_NAMES(df), variable, levels, var_class;
        int ncols = GET_LENGTH(names), type, *types;
        char *col_name, *class, *factor;

        /* variables for adding rownames */
        SEXP rownames;
        int nrows;
        char *row_name;

        /* TODO: put constraints on table after inserting everything? */


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
            var_class = GET_CLASS(variable);
            class = (var_class == R_NilValue) ? NULL : CHAR(STRING_ELT(var_class, 0));
            type = TYPEOF(variable);
            types[i] = type;

            sql_len += sprintf(g_sql_buf[0]+sql_len, ", [%s] %s", col_name, 
                    _get_column_type(class, type));

            /* add handler to insert table sql */
            _expand_buf(1, sql_len+5);
            strcpy(g_sql_buf[1]+sql_len2, ",?");
            sql_len2 += 2; 

            /* create separate table for factors decode */
            if (class != NULL && (strcmp(class, "factor") == 0 || strcmp(class, "ordered") == 0)){
                if (_create_factor_table2(iname, class, col_name)) 
                    return R_NilValue; /* dup tbl name? */

                _sqlite_exec("begin");
                levels = GET_LEVELS(variable);
                sprintf(g_sql_buf[2], "insert into [%s].[%s %s] values(?, ?);",
                        iname, class, col_name);
                res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, NULL);
                if (_sqlite_error(res)) return R_NilValue; /* dup tbl name? */

                for (j = 0; j < GET_LENGTH(levels); j++) {
                    sqlite3_reset(stmt);
                    factor = CHAR(STRING_ELT(levels, j));
                    sqlite3_bind_int(stmt, 1, j+1);
                    sqlite3_bind_text(stmt, 2, factor, -1, SQLITE_STATIC);
                    sqlite3_step(stmt);
                }
                sqlite3_finalize(stmt);
                _sqlite_exec("commit");
            }
        }
        
        _expand_buf(0,sql_len+35);
        sql_len += sprintf(g_sql_buf[0]+sql_len, ", primary key([row name]));");
        /* sql_len += sprintf(g_sql_buf[0]+sql_len, ");"); */
        res = _sqlite_exec(g_sql_buf[0]);
        if (_sqlite_error(res)) return R_NilValue; /* why? */

        /*
         * add the data in df to the sdf
         */
        rownames = getAttrib(df, R_RowNamesSymbol);
        nrows = GET_LENGTH(rownames);

        _sqlite_error(_sqlite_exec("begin"));
        sprintf(g_sql_buf[1]+sql_len2, ")");
        res = sqlite3_prepare(g_workspace, g_sql_buf[1], -1, &stmt, NULL);

        for (i = 0; i < nrows; i++) {
            sqlite3_reset(stmt);

            /* since this is coming from a real dataframe, we're assured that
             * the rownames are unique */
            if (IS_CHARACTER(rownames)) {
                row_name = CHAR(STRING_ELT(rownames, i));
                if (*row_name) /* if not empty string */
                    sqlite3_bind_text(stmt, 1, row_name, strlen(row_name), SQLITE_STATIC);
                else sqlite3_bind_int(stmt, 1, i);
            } else if (IS_INTEGER(rownames)) {
                sqlite3_bind_int(stmt, 1, INTEGER(rownames)[i]);
            } else sqlite3_bind_int(stmt, 1, i);

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
                _sqlite_exec("ROLLBACK");
                sqlite3_finalize(stmt);
                Rprintf("SQLITE ERROR: %s\n", sqlite3_errmsg(g_workspace));
                return R_NilValue; /* why? */
            }
        }
                
        sqlite3_finalize(stmt);
        _sqlite_error(_sqlite_exec("COMMIT"));

        /* create a new object representing the sdf */
        ret = _create_sdf_sexp(iname);

    } else {
        Rprintf("ERROR: unable to create a sqlite data frame.\n");
        ret = R_NilValue;
    }
        
    return ret;
}

SEXP sdf_get_names(SEXP sdf) {
    char *iname;
    iname = SDF_INAME(sdf);
    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

    int len;
    len = sprintf(g_sql_buf[0], "select * from [%s].sdf_data;", iname);

    sqlite3_stmt *stmt;
    int res, i;
   
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], len, &stmt, NULL);
    if (_sqlite_error(res)) return R_NilValue;

    SEXP ret;

    len = sqlite3_column_count(stmt)-1;
    PROTECT(ret = NEW_CHARACTER(len));

    for (i = 0; i < len; i++) {
        SET_STRING_ELT(ret, i, mkChar(sqlite3_column_name(stmt, i+1)));
    }

    sqlite3_finalize(stmt);
    UNPROTECT(1);
    return ret;
}

SEXP sdf_get_length(SEXP sdf) {
    char *iname;
    int len;
    sqlite3_stmt *stmt;
    int res;

    iname  = CHAR(STRING_ELT(_getListElement(sdf, "iname"),0));
    sprintf(g_sql_buf[0], "select * from [%s].sdf_data;", iname);
    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, NULL);
    if (_sqlite_error(res)) return R_NilValue;

    len = sqlite3_column_count(stmt)-1;
    sqlite3_finalize(stmt);
    return ScalarInteger(len);
}

/* get row count */
SEXP sdf_get_row_count(SEXP sdf) {
    char *iname = CHAR(STRING_ELT(_getListElement(sdf, "iname"),0));
    if (!USE_SDF1(iname, TRUE, FALSE)) return R_NilValue;

    int nrow;
    SEXP ret;
   
    nrow = _get_row_count2(iname, TRUE);

    if (nrow < 1) ret = R_NilValue;
    else ret = ScalarInteger(nrow);

    return ret;
}

    
SEXP sdf_import_table(SEXP _filename, SEXP _name, SEXP _sep, SEXP _quote, 
        SEXP _rownames, SEXP _colnames) {
    char *filename = CHAR_ELT(_filename, 0);
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
        Rprintf("Error: File %s does not exist.", filename);
        return R_NilValue;
    }

    /*char *sep = CHAR_ELT(_sep, 0), *quote = CHAR_ELT(_quote,0);
    char *name = CHAR_ELT(_name, 0);*/

    /* create the table */
    /* insert the stuffs */
    /* register to workspace */

    return R_NilValue;
}

SEXP sdf_get_index(SEXP sdf, SEXP row, SEXP col) {
    SEXP ret = R_NilValue;
    char *iname = SDF_INAME(sdf);
    if (!USE_SDF1(iname, TRUE, TRUE)) return R_NilValue;

    sqlite3_stmt *stmt;
    int buflen = 0, idxlen, col_cnt, row_cnt, index;
    int i, j,res;
    int *col_indices, col_index_len, *dup_indices;
    int row_index_len = 0;

    sprintf(g_sql_buf[0], "select * from [%s].sdf_data ", iname);
    res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
    if (_sqlite_error(res)) return R_NilValue;

    buflen = sprintf(g_sql_buf[0], "select [row name]");
    idxlen = LENGTH(col);
    col_cnt = sqlite3_column_count(stmt) - 1;

    /*
     * PROCESS COLUMN INDEX
     */
    if (col == R_NilValue) {
        for (i = 1; i < col_cnt+1; i++) {
            buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", 
                    sqlite3_column_name(stmt, i)); 
        }
        _expand_buf(0, buflen+20+strlen(iname));
        buflen += sprintf(g_sql_buf[0]+buflen, " from [%s].sdf_data", iname);
        col_index_len = col_cnt;
        col_indices = (int *)R_alloc(idxlen, sizeof(int));
        dup_indices = (int *)R_alloc(idxlen, sizeof(int));
        for (i = 0; i < col_cnt; i++) { col_indices[i] = i; dup_indices[i] = 0; }
    } else if (col == R_NilValue || idxlen < 1) {
        sqlite3_finalize(stmt);
        return R_NilValue;
    } else if (IS_NUMERIC(col)) {
        col_indices = (int *)R_alloc(idxlen, sizeof(int));
        dup_indices = (int *)R_alloc(idxlen, sizeof(int));
        col_index_len = 0;

        for (i = 0; i < idxlen; i++) {
            /* no need to correct for 0 base because col 0 is [row name] */
            index = ((int) REAL(col)[i]); 
            if (index > col_cnt) {
                sqlite3_finalize(stmt);
                Rprintf("Error: undefined columns selected\n");
                return R_NilValue;
            } else if (index > 0) {
                buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", 
                        sqlite3_column_name(stmt,index));
                dup_indices[col_index_len] = 0;
                for (j = 0; j < col_index_len; j++) {
                    if (col_indices[j] == index) dup_indices[col_index_len]++;
                }
                col_indices[col_index_len++] = index;
            } else if (index < 0) {
                sqlite3_finalize(stmt);
                Rprintf("Error: negative indices not supported.\n");
                return R_NilValue;
            }
        }

        if (col_index_len == 0) {
            sqlite3_finalize(stmt);
            Rprintf("Error: no indices detected??\n");
            return R_NilValue;
        } else { 
            _expand_buf(0, buflen+20+strlen(iname));
            buflen += sprintf(g_sql_buf[0]+buflen, " from [%s].sdf_data", iname);
        }
    } else if (IS_INTEGER(col)) {
        /* identical logic with IS_NUMERIC, except that we don't have to cast
         * stuffs from SEXP col. the alternative is doing if idxlen times. */
        col_indices = (int *)R_alloc(idxlen, sizeof(int));
        dup_indices = (int *)R_alloc(idxlen, sizeof(int));
        col_index_len = 0;

        for (i = 0; i < idxlen; i++) {
            /* no need to correct for 0 base because col 0 is [row name] */
            index = INTEGER(col)[i];
            if (index > col_cnt) {
                sqlite3_finalize(stmt);
                Rprintf("Error: undefined columns selected\n");
                return R_NilValue;
            } else if (index > 0) {
                buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", 
                        sqlite3_column_name(stmt,index));
                dup_indices[col_index_len] = 0;
                for (j = 0; j < col_index_len; j++) {
                    if (col_indices[j] == index) dup_indices[col_index_len]++;
                }
                col_indices[col_index_len++] = index;
            } else if (index < 0) {
                sqlite3_finalize(stmt);
                Rprintf("Error: negative indices not supported.\n");
                return R_NilValue;
            }
        }

        if (col_index_len == 0) {
            sqlite3_finalize(stmt);
            Rprintf("Error: no indices detected??\n");
            return R_NilValue;
        } else { 
            _expand_buf(0, buflen+20+strlen(iname));
            buflen += sprintf(g_sql_buf[0]+buflen, " from [%s].sdf_data", iname);
        }
    } else if (IS_LOGICAL(col)) {
        /* recycling stuff, so max column output is the # of cols in the df */
        col_indices = (int *)R_alloc(col_cnt, sizeof(int));
        dup_indices = (int *)R_alloc(col_cnt, sizeof(int));
        col_index_len = 0;
        for (i = 0; i < col_cnt; i++) {
            if (LOGICAL(col)[i%idxlen]) {
                buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", 
                        sqlite3_column_name(stmt,i+1));
                dup_indices[col_index_len] = 0;
                col_indices[col_index_len++] = i;
            }
        }

        if (col_index_len == 0) {
            sqlite3_finalize(stmt);
            Rprintf("Warning: no column selected.\n");
            return R_NilValue;
        } else { 
            _expand_buf(0, buflen+20+strlen(iname));
            buflen += sprintf(g_sql_buf[0]+buflen, " from [%s].sdf_data", iname);
        }
    } else {
        sqlite3_finalize(stmt);
        Rprintf("Error: don't know how to handle column index.\n");
        return R_NilValue;
    }

    /* 
     * PROCESS ROW INDEX
     */
    idxlen = LENGTH(row);
    if (row == R_NilValue) {
        if (col_index_len == 1) {
            ret = sdf_get_variable(sdf, mkString(sqlite3_column_name(stmt,col_indices[0])));
            sqlite3_finalize(stmt);
            return ret;
        } if (col_index_len > 1) {
            /* create a new SDF, logic similar to sdf_create_sdf */

            /* find a new name. data<n> ? */
            char *iname2;
            const char *colname;
            int namelen, sql_len, sql_len2;

            iname2 = _create_sdf_skeleton1(R_NilValue, &namelen, TRUE);

            /* create sdf_data table */
            sql_len = sprintf(g_sql_buf[1], "create table [%s].sdf_data ([row name] text", iname2);
            sql_len2 = 0;

            for (i = 0; i < col_index_len; i++) {
                colname = sqlite3_column_name(stmt, col_indices[i]);
                if (dup_indices[i] == 0) {
                    sql_len += sprintf(g_sql_buf[1]+sql_len, ", [%s] %s", colname,
                            sqlite3_column_decltype(stmt, col_indices[i]));
                } else {
                    /* un-duplicate col names by appending num, just like R */
                    sql_len += sprintf(g_sql_buf[1]+sql_len, ", [%s.%d] %s", colname,
                            dup_indices[i], 
                            sqlite3_column_decltype(stmt, col_indices[i]));
                }

                /* deal with possibly factor columns */
                if (_is_factor2(iname, "factor", colname)) {
                    /* found a factor column */

                    if (dup_indices[i] == 0)  {
                        _copy_factor_levels2("factor", iname, colname, iname2, colname);
                    } else {
                        sprintf(g_sql_buf[3], "%s.%d", colname, dup_indices[i]);
                        _copy_factor_levels2("factor", iname, colname, iname2, 
                                g_sql_buf[3]);
                    }
                }

                /* and deal with ordered factors too... life is hard, then you die */
                if (_is_factor2(iname, "ordered", colname)) {
                    
                    if (dup_indices[i] == 0)  {
                        _copy_factor_levels2("ordered", iname, colname, iname2, 
                                colname);
                    } else {
                        sprintf(g_sql_buf[3], "%s.%d", colname, dup_indices[i]);
                        _copy_factor_levels2("ordered", iname, colname, iname2, 
                                g_sql_buf[3]);
                    }
                }
                                
            }

            /* don't need it anymore */
            sqlite3_finalize(stmt);

            /* no table index created, that's for v2 I hope*/
            sql_len += sprintf(g_sql_buf[1]+sql_len, ");");
            res = _sqlite_exec(g_sql_buf[1]);
            if (_sqlite_error(res)) { Rprintf("here? %s\n", g_sql_buf[1]); return R_NilValue; }

            /* insert data (with row names). buf[0] contains a select */
            sprintf(g_sql_buf[1], "insert into [%s].sdf_data %s", iname2,
                    g_sql_buf[0]);
            res = _sqlite_exec(g_sql_buf[1]);
            if (_sqlite_error(res)) return R_NilValue;

            /* add new sdf to workspace */
            sprintf(g_sql_buf[0], "%s.db", iname2);
            res = _add_sdf1(g_sql_buf[0], iname2);
            if (_sqlite_error(res)) return R_NilValue;

            /* remove protection */
            UNUSE_SDF2(iname2);

            /* create SEXP for the SDF */
            ret = _create_sdf_sexp(iname2);
            return ret;
        } else { sqlite3_finalize(stmt); return R_NilValue; }
    } else if (row == R_NilValue || idxlen < 1) {
        sqlite3_finalize(stmt);
        return R_NilValue;
    } else if (IS_NUMERIC(row)) {
        sqlite3_finalize(stmt);

        sprintf(g_sql_buf[1], "[%s].sdf_data", iname);
        row_cnt = _get_row_count2(g_sql_buf[1], 0);

        /* append " limit ?,1" to the formed select statement */
        _expand_buf(0, buflen+10);
        buflen += sprintf(g_sql_buf[0]+buflen, " limit ?,1");
        /*Rprintf("sql [%d]: %s\n", buflen, g_sql_buf[0]); */

        index = ((int) REAL(row)[0]) - 1;
        if (idxlen < 0 && idxlen == 1) return R_NilValue;

        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
        if (_sqlite_error(res)) return R_NilValue;

        sqlite3_bind_int(stmt, 1, index);
        res = sqlite3_step(stmt);

        /* create data frame */
        ret = _setup_df_sexp1(stmt, iname, col_index_len, idxlen, dup_indices);
        if (ret == R_NilValue) { sqlite3_finalize(stmt); return R_NilValue; }

        /* put data in it */
        if (index >= 0 && index < row_cnt) {
            _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
        }

        for (i = 1; i < idxlen; i++) {
            sqlite3_reset(stmt);
            index = ((int) REAL(row)[i]) - 1;
            if (index >= 0 && index < row_cnt) {
                sqlite3_bind_int(stmt, 1, index);
                res = sqlite3_step(stmt);
                _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
            }
        }
        sqlite3_finalize(stmt);

        /* shrink vectors */
        if (row_index_len < idxlen) {
            for (i = 0; i < col_index_len; i++) {
                SET_VECTOR_ELT(ret, i, 
                        _shrink_vector(VECTOR_ELT(ret, i), row_index_len));
            }
        }

        UNPROTECT(1); /* for the ret from _setup_df_sexp1 */
    } else if (IS_INTEGER(row)) {
        /* same stuff as with IS_INTEGER, except for setting of var index*/
        sqlite3_finalize(stmt);

        sprintf(g_sql_buf[1], "[%s].sdf_data", iname);
        row_cnt = _get_row_count2(g_sql_buf[1], 0);

        /* append " limit ?,1" to the formed select statement */
        _expand_buf(0, buflen+10);
        buflen += sprintf(g_sql_buf[0]+buflen, " limit ?,1");

        index = INTEGER(row)[0] - 1;
        if (idxlen < 0 && idxlen == 1) return R_NilValue;

        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
        if (_sqlite_error(res)) return R_NilValue;

        sqlite3_bind_int(stmt, 1, index);
        res = sqlite3_step(stmt);

        /* create data frame */
        ret = _setup_df_sexp1(stmt, iname, col_index_len, idxlen, dup_indices);
        if (ret == R_NilValue) { sqlite3_finalize(stmt); return R_NilValue; }

        /* put data in it */
        if (index >= 0 && index < row_cnt) {
            _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
        }

        for (i = 1; i < idxlen; i++) {
            sqlite3_reset(stmt);
            index = INTEGER(row)[i] - 1;
            if (index >= 0 && index < row_cnt) {
                sqlite3_bind_int(stmt, 1, index);
                res = sqlite3_step(stmt); 
                _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
            }
        }
        sqlite3_finalize(stmt);

        /* shrink vectors */
        if (row_index_len < idxlen) {
            for (i = 0; i < col_index_len; i++) {
                SET_VECTOR_ELT(ret, i, 
                        _shrink_vector(VECTOR_ELT(ret, i), row_index_len));
            }
        }
        UNPROTECT(1); /* for the ret from _setup_df_sexp1 */
    } else if (IS_LOGICAL(row)) {
        /* */
        sqlite3_finalize(stmt);
        int est_row_cnt = 0;

        sprintf(g_sql_buf[1], "[%s].sdf_data", iname);
        row_cnt = _get_row_count2(g_sql_buf[1], 0);

        /* append " limit ?,1" to the formed select statement */
        _expand_buf(0, buflen+10);
        buflen += sprintf(g_sql_buf[0]+buflen, " limit ?,1");

        /* find if there is any TRUE element in the vector */
        for (i = 0; i < idxlen && i < row_cnt; i++) {
            if (LOGICAL(row)[i]) break;
        }

        if (i < idxlen && i < row_cnt) {
            est_row_cnt = (idxlen-i) * (row_cnt/idxlen);
            if (row_cnt%idxlen > i) est_row_cnt += (row_cnt%idxlen - i);

            res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
            if (_sqlite_error(res)) return R_NilValue;

            sqlite3_bind_int(stmt, 1, i);
            sqlite3_step(stmt);
            ret = _setup_df_sexp1(stmt, iname, col_index_len, est_row_cnt, dup_indices);
            _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
            
            for (i++ ; i < row_cnt; i++) {
                if (LOGICAL(row)[i%idxlen]) {
                    sqlite3_reset(stmt);
                    sqlite3_bind_int(stmt, 1, i);
                    sqlite3_step(stmt);
                    _add_row_to_df(ret, stmt, row_index_len++, col_index_len);
                }
            }
        }

        sqlite3_finalize(stmt);

        /* shrink vectors */
        if (est_row_cnt > 0 && row_index_len < est_row_cnt) {
            for (i = 0; i < col_index_len; i++) {
                SET_VECTOR_ELT(ret, i, 
                        _shrink_vector(VECTOR_ELT(ret, i), row_index_len));
            }
        }

        UNPROTECT(1); /* for the ret from _setup_df_sexp1 */
    }

    UNUSE_SDF2(iname);

    if (ret != R_NilValue && col_index_len > 1) {
        SEXP class = mkString("data.frame");
        SET_CLASS(ret, class);
        _set_rownames2(ret);
    } else if (ret != R_NilValue && col_index_len == 1) {
        ret = VECTOR_ELT(ret, 0);
    }

    return ret;
}

SEXP sdf_rbind(SEXP sdf, SEXP data) {
    char *class = CHAR_ELT(GET_CLASS(data), 0);
    char *iname = SDF_INAME(sdf), *colname, *rowname;
    SEXP ret = R_NilValue, col, names, levels, rownames;
    int ncols, buflen, buflen2, i, j, res, nrows, *types;
    sqlite3_stmt *stmt;
    const char *type;

    if (strcmp(class,"data.frame") == 0) {
        /* check table names, types. variables may not be in same order, as
         * long as names and types are identical. */
        names = GET_NAMES(data);
        ncols = GET_LENGTH(names);

        buflen = sprintf(g_sql_buf[0], "select 1");
        for (i = 0; i < ncols; i++) {
            _expand_buf(0, buflen + 100);
            buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", CHAR_ELT(names, i));
        }
        _expand_buf(0, buflen + 100);
        sprintf(g_sql_buf[0]+buflen, " from [%s].sdf_data", iname);

        res = sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, NULL);
        if (res != SQLITE_OK) { /* column name mismatch */
            Rprintf("Error: column names should exactly match the names of the SDF.\n");
            return ret;
        }

        /* now check the type */
        types = (int *)R_alloc(ncols, sizeof(int));
        for (i = 0; i< ncols; i++) {
            type = sqlite3_column_decltype(stmt, i+1); /* +1 bec 1st col is "1" */

            col = VECTOR_ELT(data, i);
            class = CHAR_ELT(GET_CLASS(col),0);

            if (strcmp(type, "double") == 0 && IS_NUMERIC(col)) types[i] = SQLITE_FLOAT;
            else if (strcmp(type, "text") == 0 && IS_CHARACTER(col)) types[i] = SQLITE_TEXT;
            else if (strcmp(type, "int") == 0) { 
                types[i] = SQLITE_INTEGER;
                if (isUnordered(col)) {  /* test if factor */
                    colname = CHAR_ELT(names, i);
                    sqlite3_stmt *stmt2;
                    if (!_is_factor2(iname, "factor", colname)) break;
                    sprintf(g_sql_buf[1], "[%s].[factor %s]", iname, colname);
                    nrows = _get_row_count2(g_sql_buf[1], 0);

                    /* levels in sdf must be >= levels in df */
                    levels = GET_LEVELS(col);
                    if (nrows < GET_LENGTH(levels)) {
                        Rprintf("Error: The data frame variable %s has more levels than"
                                " its SDF counterpart.\n", colname);
                        break;
                    }

                    sprintf(g_sql_buf[2], "select label from %s order by level", g_sql_buf[1]);
                    res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, NULL);
                    _sqlite_error(res);
                    
                    for (j = 0; j < GET_LENGTH(levels); j++)  {
                        sqlite3_step(stmt2);
                        if (strcmp((char *)sqlite3_column_text(stmt2, 0), CHAR_ELT(levels, j)) != 0) {
                            Rprintf("Error: Level label mismatch in variable %s, [%s] in data frame"
                                    " vs [%s] in SDF.\n", colname, (char *)sqlite3_column_text(stmt2, 0),
                                    CHAR_ELT(levels, j));
                            sqlite3_finalize(stmt2);
                            break;
                        }
                    }
                    sqlite3_finalize(stmt2);
                } else if (isOrdered(col)) { /* same with ordered */
                    colname = CHAR_ELT(names, i);
                    sqlite3_stmt *stmt2;
                    if (!_is_factor2(iname, "ordered", colname)) break;
                    sprintf(g_sql_buf[1], "[%s].[ordered %s]", iname, colname);
                    nrows = _get_row_count2(g_sql_buf[1], 0);

                    /* levels in sdf must be >= levels in df */
                    levels = GET_LEVELS(col);
                    if (nrows < GET_LENGTH(levels)) {
                        Rprintf("Error: The data frame variable %s has more levels than"
                                " its SDF counterpart.\n", colname);
                        break;
                    }

                    sprintf(g_sql_buf[2], "select label from %s order by level", g_sql_buf[1]);
                    sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt2, NULL);
                    
                    for (j = 0; j < GET_LENGTH(levels); j++)  {
                        sqlite3_step(stmt2);
                        if (strcmp((char *)sqlite3_column_text(stmt2, 0), CHAR_ELT(levels, j)) != 0) {
                            Rprintf("Error: Level label mismatch in variable %s, [%s] in data frame"
                                    " vs [%s] in SDF.\n", colname, (char *)sqlite3_column_text(stmt2, 0),
                                    CHAR_ELT(levels, j));
                            sqlite3_finalize(stmt2);
                            break;
                        }
                    }
                    sqlite3_finalize(stmt2);
                }  /* else if class == "ordered" */
            } /* if integer */
        } /* for loop for column type checking */

        sqlite3_finalize(stmt);
        if (i < ncols) {
            Rprintf("Error: type mismatch at variable %s\n", CHAR_ELT(names, i));
            return ret;
        }

        /* create the insert statement */
        buflen = sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name]", iname);
        buflen2 = sprintf(g_sql_buf[1], "values(?");
        for (i = 0; i < ncols; i++) {
            _expand_buf(0, buflen+100);
            buflen += sprintf(g_sql_buf[0]+buflen, ",[%s]", CHAR_ELT(names, i));
            _expand_buf(1, buflen2+5);
            buflen2 += sprintf(g_sql_buf[1]+buflen2, ",?");
        }

        _expand_buf(2, buflen + buflen2 + 10);
        sprintf(g_sql_buf[2], "%s) %s)", g_sql_buf[0], g_sql_buf[1]);
        res = sqlite3_prepare(g_workspace, g_sql_buf[2], -1, &stmt, NULL);
        if (_sqlite_error(res)) return ret;

        /* finally, add rows */
        nrows = LENGTH(GET_ROWNAMES(data));
        rownames = GET_ROWNAMES(data);

        _sqlite_exec("begin");
        for (i = 0; i < nrows; i++) {
            sqlite3_bind_text(stmt, 1, CHAR_ELT(rownames, i), -1, SQLITE_STATIC);
            for (j = 0; j < ncols; j++) {
                col = VECTOR_ELT(data, j);
                switch(types[j]) {
                    case SQLITE_FLOAT:
                        sqlite3_bind_double(stmt, j+2, REAL(col)[i]); break;
                    case SQLITE_TEXT:
                        sqlite3_bind_text(stmt, j+2, CHAR_ELT(col, i), -1, SQLITE_STATIC); break;
                    case SQLITE_INTEGER:
                    default:  /* runtime error if we're wrong */
                        sqlite3_bind_int(stmt, j+2, INTEGER(col)[i]); break;
                }
            }

            res = sqlite3_step(stmt);
            if (res != SQLITE_DONE) { /* row name pk conflict */
                rowname = CHAR_ELT(rownames, i);
                for (j = 1; res != SQLITE_DONE; j++) { /* _normally_, j shouldn't overflow */
                    sqlite3_reset(stmt);
                    sprintf(g_sql_buf[2], "%s-%d", rowname, j);
                    sqlite3_bind_text(stmt, 1, g_sql_buf[2], -1, SQLITE_STATIC);
                    res = sqlite3_step(stmt);
                }
            }
            sqlite3_reset(stmt);
        }

        sqlite3_finalize(stmt);
        _sqlite_exec("begin");
        ret = sdf;
    } else if (strcmp(class,"sqlite.data.frame") == 0) {
    }

    return ret;
}


SEXP sdf_get_iname(SEXP sdf) {
    return _getListElement(sdf, "iname");
}
