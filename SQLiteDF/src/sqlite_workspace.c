#include <stdio.h>
#include <string.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

#define __SQLITE_WORKSPACE__
#include "sqlite_dataframe.h"

sqlite3 *g_workspace = NULL;
char *g_sql_buf[NBUFS];
int g_sql_buf_sz[NBUFS];

sqlite3* _is_sqlitedb(char *filename) {
    sqlite3 *db;
    int res;
    res = sqlite3_open(filename, &db);
    if (res != SQLITE_OK) { goto is_sqlitedb_FAIL; }

    sqlite3_stmt *stmt; char *sql = "select * from sqlite_master";
    res = sqlite3_prepare(db, sql, strlen(sql), &stmt, 0);
    if (stmt != NULL) sqlite3_finalize(stmt);
    /*char **result_set;
    char nrow, ncol;
    res = sqlite3_get_table(db, "select * from sqlite_master limit 0", 
            &result_set, &nrow, &ncol, NULL);
    sqlite3_free_table(result_set);*/
    if (res != SQLITE_OK) goto is_sqlitedb_FAIL;

    return db;

is_sqlitedb_FAIL:
    sqlite3_close(db);
    return NULL;
}

sqlite3* _is_workspace(char *filename) {
    sqlite3* db = _is_sqlitedb(filename); 

    if (db != NULL) {
        sqlite3_stmt *stmt;
        char *sql = "select * from workspace";
        int res = sqlite3_prepare(db, sql, strlen(sql), &stmt, 0), ncols;
        if ((res != SQLITE_OK) || /* no workspace table */
              ((ncols = sqlite3_column_count(stmt)) != 2) ||
              /* below also checks the ordering of the columns */
              (strcmp(sqlite3_column_name(stmt, 0), "filename") != 0) ||
              (strcmp(sqlite3_column_decltype(stmt, 0), "text") != 0) ||
              (strcmp(sqlite3_column_name(stmt, 1), "internal_name") != 0) ||
              (strcmp(sqlite3_column_decltype(stmt, 1), "text") != 0)) {
            sqlite3_finalize(stmt); sqlite3_close(db); db = NULL;
        } else {
            sqlite3_finalize(stmt);
        }
    }

    return db;
}

int _is_sdf(char *filename) {
    sqlite3* db = _is_sqlitedb(filename); 
    int ret = (db != NULL);

    if (ret) {
        sqlite3_stmt *stmt;
        char *sql = "select * from sdf_attributes";
        int res = sqlite3_prepare(db, sql, strlen(sql), &stmt, NULL), ncols;
        ret = ((res == SQLITE_OK) && /* no attribute table */
               ((ncols = sqlite3_column_count(stmt)) == 2) &&
               (strcmp(sqlite3_column_name(stmt, 0), "attr") == 0) &&
               (strcmp(sqlite3_column_decltype(stmt, 0), "text") == 0) &&
               (strcmp(sqlite3_column_name(stmt, 1), "value") == 0) &&
               (strcmp(sqlite3_column_decltype(stmt, 1), "text") == 0)) ;
        
        if (!ret) goto _is_sdf_cleanup;
        sqlite3_finalize(stmt); 
        
        
        sql = "select * from sdf_data";
        res = sqlite3_prepare(db, sql, strlen(sql), &stmt, NULL);
        ret = (res == SQLITE_OK);  /* if not, missing data table */

_is_sdf_cleanup:
        sqlite3_finalize(stmt);
        sqlite3_close(db);
    }

    return ret;
}

void _delete_sdf2(char *iname) {
    sprintf(g_sql_buf[2], "delete from workspace where internal_name='%s';", iname);
    _sqlite_exec(g_sql_buf[2]);
}


SEXP sdf_init_workspace() {
    int file_idx = 0, i;
    char *basename = "workspace", *filename;
    FILE *f;
    SEXP ret;

    /* initialize sql_buf */
    for (i = 0; i < NBUFS; i++) {
        if (g_sql_buf[i] == NULL) {
            g_sql_buf_sz[i] = 1024;
            g_sql_buf[i] = Calloc(g_sql_buf_sz[i], char);
        }
    }

    /*
     * check for workspace.db, workspace1.db, ..., workspace9999.db if they
     * are valid workspace file. if one is found, use that as the workspace.
     */
    filename = R_alloc(18, sizeof(char)); /* workspace10000.db\0 */
    sprintf(filename, "%s.db", basename);
    while(_file_exists(filename) && file_idx < 10000) {
        if ((g_workspace = _is_workspace(filename)) != NULL) break;
        /* warn("%s is not a workspace", filename) */
        sprintf(filename, "%s%d.db", basename, ++file_idx);
    }

    PROTECT(ret = NEW_LOGICAL(1));
    if ((g_workspace == NULL) && (file_idx < 10000)) {
        /* no workspace found but there are still "available" file name */
        /* if (file_idx) warn("workspace will be stored at #{filename}") */
        sqlite3_open(filename, &g_workspace);
        _sqlite_exec("create table workspace(filename text, internal_name text)");
        LOGICAL(ret)[0] = TRUE;
    } else if (g_workspace != NULL) {
        /* a valid workspace has been found, load each of the tables */
        int res, nrows, ncols; 
        char **result_set, *fname, *iname;
        
        res = sqlite3_get_table(g_workspace, "select * from workspace", 
                &result_set, &nrows, &ncols, NULL);
        
        if (res == SQLITE_OK && nrows >= 1 && ncols == 2) {
            for (i = 1; i <= nrows; i++) {
                fname = result_set[i*ncols]; iname = result_set[i*ncols+1];
                
                if (!_file_exists(fname)) {
                    Rprintf("Warning: SDF %s does not exist.\n", iname);
                    _delete_sdf2(iname);
                    continue;
                }

                if (!_is_sdf(fname)) {
                    Rprintf("Warning: %s is not a valid SDF.\n", fname);
                    _delete_sdf2(iname);
                    continue;
                }

                sprintf(g_sql_buf[0], "attach '%s' as %s", fname, iname);
                _sqlite_exec(g_sql_buf[0]);
            }
        }
        sqlite3_free_table(result_set);

        LOGICAL(ret)[0] = TRUE;
    } else { /* can't find nor create workspace */
        LOGICAL(ret)[0] = FALSE;
    }

    UNPROTECT(1);
    return ret;
}
        

    
SEXP sdf_finalize_workspace() {
    SEXP ret;
    PROTECT(ret = NEW_LOGICAL(1)); 
    LOGICAL(ret)[0] = (sqlite3_close(g_workspace) == SQLITE_OK);
    UNPROTECT(1);
    return ret;
} 


SEXP sdf_list_sdfs(SEXP pattern) {
    SEXP ret;
    char **result;
    int nrow, ncol, res, i;

    if (TYPEOF(pattern) != STRSXP) {
        res = sqlite3_get_table(g_workspace, "select internal_name from workspace",
                &result, &nrow, &ncol, NULL);
    } else {
        /* since internal_names must be a valid r symbol, 
           did not check for "'" */
        sprintf(g_sql_buf[0], "select internal_name from workspace where "
                "internal_name like '%s\%'", CHAR(STRING_ELT(pattern, 0)));
        res = sqlite3_get_table(g_workspace, g_sql_buf[0], &result, &nrow,
                &ncol, NULL);
    }

    if (_sqlite_error(res)) return R_NilValue;
    PROTECT(ret = NEW_CHARACTER(nrow));
    
    for (i = 0; i < nrow; i++) SET_STRING_ELT(ret, i, mkChar(result[i+1]));

    sqlite3_free_table(result);
    UNPROTECT(1);
    return ret;
}

SEXP sdf_get_sdf(SEXP name) {    
    if (TYPEOF(name) != STRSXP) {
        Rprintf("Error: Argument must be a string containing the SDF name.\n");
        return R_NilValue;
    }

    char *iname = CHAR(STRING_ELT(name, 0));
    SEXP ret;
    sqlite3_stmt *stmt;
    int res;

    res = sqlite3_prepare(g_workspace, "select * from workspace where internal_name=?",
            45, &stmt, NULL);
    if (_sqlite_error(res)) return R_NilValue;

    sqlite3_bind_text(stmt, 1, iname, strlen(iname), SQLITE_STATIC);
    res = sqlite3_step(stmt);

    if (res == SQLITE_ROW) ret = _create_sdf_sexp(iname);
    else ret = R_NilValue;
    
    sqlite3_finalize(stmt);

    return ret;
}
