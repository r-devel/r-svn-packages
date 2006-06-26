#include <stdio.h>
#include <string.h>
#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include "sqlite3.h"

sqlite3 *workspace = NULL;

sqlite3* _is_sqlitedb(filename) {
    sqlite3 *db;
    int res;
    res = sqlite3_open(filename, &db);
    if (res != SQLITE_OK) { goto is_sqlitedb_FAIL; }

    sqlite3_stmt *stmt; char *sql = "select * from workspace";
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
    sqlite_close(db);
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
              (strcmp(sqlite3_column_decl(stmt, 0), "text") != 0) ||
              (strcmp(sqlite3_column_name(stmt, 1), "internal_name") != 0) ||
              (strcmp(sqlite3_column_decl(stmt, 1), "text") != 0)) {
            sqlite3_finalize(stmt); sqlite3_close(db); db = NULL;
        } else {
            sqlite3_finalize(stmt);
        }
    }

    return db;
}

int _is_sdf(char *filename) {
    sqlite3* db = _is_sqlitedb(filename); 

    if (db != NULL) {
        sqlite3_stmt *stmt;
        char *sql = "select * from sdf_attributes";
        int res = sqlite3_prepare(db, sql, strlen(sql), &stmt, 0), ncols;
        if ((res != SQLITE_OK) || /* no workspace table */
              ((ncols = sqlite3_column_count(stmt)) != 2) ||
              (strcmp(sqlite3_column_name(stmt, 0), "attr") != 0) ||
              (strcmp(sqlite3_column_decl(stmt, 0), "text") != 0) ||
              (strcmp(sqlite3_column_name(stmt, 1), "value") != 0) ||
              (strcmp(sqlite3_column_decl(stmt, 1), "text") != 0)) {
            sqlite3_finalize(stmt); sqlite3_close(db); return NULL;
        } 
        
        /* we also check the contents of sdf_attributes */
        sqlite3_finalize(stmt);
        
    }

    return db;
}

int _file_exists(char *filename) {
    FILE *f; int ret = FALSE;
    if ((f = fopen(filename, "rb")) != null) { fclose(f); ret = TRUE; }
    return ret;
}


SEXP sdf_init_workspace() {
    int file_idx = 0;
    char *basename = "workspace", *filename;
    FILE *f;
    SEXP ret;

    /* check for workspace.db, workspace1.db, ..., workspace9999.db if they
     * are valid workspace file. if one is found, use that as the workspace.
     */
    filename = R_alloc(18, sizeof(char)); /* workspace10000.db\0 */
    sprintf(filename, "%s.db", basename);
    while(_file_exists(filename) && file_idx < 10000); {
        if ((workspace = _is_workspace(filename)) != NULL) break;
        /* warn("%s is not a workspace", filename) */
        sprintf(filename, "%s%d.db", basename, ++file_idx);
    }

    PROTECT(ret = NEW_LOGICAL(1));
    if ((workspace == NULL) && (file_idx < 10000)) {
        /* no workspace found but there are still "available" file name */
        /* if (file_idx) warn("workspace will be stored at #{filename}") */
        sqlite3_open(filename, &workspace);
        sqlite3_exec(workspace, "create table workspace(filename text, internal_name text)", NULL, NULL, NULL);
        LOGICAL(ret)[0] = TRUE;
    } else if (file_idx >= 10000) {
        /* a valid workspace has been found, load each of the tables */
        int res, nrows, ncols; 
        char **result_set, *fname, *iname;
        
        res = sqlite3_get_table(workspace, "select * from workspace", 
                &result_set, &nrows, &ncols, NULL);
        
        if (res == SQLITE_OK && nrows > 1 && ncols == 2) {
            sqlite3_stmt *stmt;
            for (i = 1; i < nrows; i++) {
                fname = result_set[i*ncols]; iname = result_set[i*ncols+1];
                
                if (!_file_exists(fname)) {
                    /* warn('table does not exist: #{fname}'); */
                    continue;
                }

                if (!_is_sdf(fname)) {
                }

                sql = sqlite3_mprintf("attach '%s' as %s", fname, iname);
                res = sqlite3_prepare(db, sql, strlen(sql), &stmt, 0);
                if (res != SQL_OK) { 
                    /* warn('not an sqlite-db file */
                    continue;
                }
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
        

    


