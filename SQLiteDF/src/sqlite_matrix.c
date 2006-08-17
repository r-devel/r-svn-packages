#define __SQLITE_MATRIX__
#include "sqlite_dataframe.h"

/* SEXP SDFMatrixDimSymbol =  */
/****************************************************************************
 * utilities
 ****************************************************************************/
static void _installAttrib(SEXP obj, SEXP name, SEXP value) {
    /* to bypass dim & dimnames property checking */
    SEXP attr, cur_attr;
    PROTECT(attr = NEW_LIST(1));
    SETCAR(attr, value);
    SET_TAG(attr, name);
    if (ATTRIB(obj) == R_NilValue) SET_ATTRIB(obj, attr);
    else {
        cur_attr = nthcdr(ATTRIB(obj), length(ATTRIB(obj)) - 1);
        SETCDR(cur_attr, attr);
    }
    UNPROTECT(1);
}

/****************************************************************************
 * SMAT FUNCTIONS
 ****************************************************************************/
SEXP sdf_as_matrix(SEXP sdf, SEXP name) {
    char *iname, *mat_iname, *type;
    const char *dectype, *colname;
    sqlite3_stmt *stmt, *stmt2;
    int ncols, nrows, i;
    SEXP ret, tmp, names;

    iname = SDF_INAME(sdf);
    mat_iname = CHAR_ELT(name, 0);

    if (!USE_SDF1(iname, TRUE, TRUE)) return R_NilValue;

    /* check column types, and determine matrix mode */
    sprintf(g_sql_buf[0], "select * from [%s].sdf_data", iname);
    sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt, 0);
    sqlite3_step(stmt);

    ncols = sqlite3_column_count(stmt);
    type = NULL;
    for (i = 1; i < ncols; i++) {
        dectype = sqlite3_column_decltype(stmt, i);
        colname = sqlite3_column_name(stmt, i);
        if (strcmp(dectype, "double") == 0) {
            type = "double";
        } else if (strcmp(dectype, "int") == 0) {
            if (_is_factor2(iname, "factor", colname) || _is_factor2(iname, "ordered", colname)) {
                type = "text"; break;
            }
            type = "double";
        } else if (strcmp(dectype, "bit") == 0) {
            if (type == NULL) type = "bit";
        } else if (strcmp(dectype, "text") == 0) {
            type = "text"; break;
        }
    }


    /* create sdf with 1 column */
    mat_iname = _create_svector1(name, type, NULL, TRUE);

    /* create and set sdf_matrix_rownames & sdf_matrix_colnames */
    sprintf(g_sql_buf[0], "create table [%s].sdf_matrix_rownames(name text)", mat_iname);
    _sqlite_error(_sqlite_exec(g_sql_buf[0]));
    sprintf(g_sql_buf[0], "insert into [%s].sdf_matrix_rownames "
            "select [row name] from [%s].sdf_data", mat_iname, iname);
    _sqlite_error(_sqlite_exec(g_sql_buf[0]));
    sprintf(g_sql_buf[0], "create table [%s].sdf_matrix_colnames(name text)", mat_iname);
    _sqlite_error(_sqlite_exec(g_sql_buf[0]));
    sprintf(g_sql_buf[0], "insert into [%s].sdf_matrix_colnames values (?)", mat_iname);
    sqlite3_prepare(g_workspace, g_sql_buf[0], -1, &stmt2, 0);

    /* store names in a SEXP */
    names = NEW_CHARACTER(ncols - 1);
    
    /* insert cast-ed values, and column names */
    nrows = _get_row_count2(iname, FALSE);
    _sqlite_error(_sqlite_exec("begin"));
    for (i = 1; i < ncols; i++) {
        colname = sqlite3_column_name(stmt, i);

        sqlite3_reset(stmt2);
        sqlite3_bind_text(stmt2, 1, colname, -1, SQLITE_STATIC);
        sqlite3_step(stmt2);

        SET_STRING_ELT(names, i-1, mkChar(colname));

        if ((_is_factor2(iname, "ordered", colname) && ((dectype = "ordered"))) || 
            ( _is_factor2(iname, "factor", colname) && ((dectype = "factor")))) {
            sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
                    "select [row name]|| %d, [%s].[%s %s].label from [%s].sdf_data "
                    "join [%s].[%s %s] on [%s].sdf_data.[%s]=[%s].[%s %s].level",
                    mat_iname, i, iname, dectype, colname, iname, iname, dectype, 
                    colname, iname, colname, iname, dectype, colname);
        } else {
            sprintf(g_sql_buf[0], "insert into [%s].sdf_data([row name], V1) "
                    "select [row name] || %d, cast([%s] as %s) from [%s].sdf_data", mat_iname, i,
                    colname, type, iname);
        }
        _sqlite_error(_sqlite_exec(g_sql_buf[0]));
    }
    sqlite3_finalize(stmt2);
    sqlite3_finalize(stmt);
    _sqlite_error(_sqlite_exec("commit"));

    /* return smat sexp */
    PROTECT(ret = NEW_LIST(2)); i = 2; /* 1 for names sexp above */
    SET_VECTOR_ELT(ret, 0, mkString(iname));
    SET_VECTOR_ELT(ret, 1, mkString("V1"));

    /* set smat data name */
    PROTECT(tmp = NEW_CHARACTER(2)); i++;
    SET_STRING_ELT(tmp, 0, mkChar("iname"));
    SET_STRING_ELT(tmp, 1, mkChar("varname"));
    SET_NAMES(ret, tmp);

    /* set class */
    PROTECT(tmp = NEW_CHARACTER(2)); i++;
    SET_STRING_ELT(tmp, 0, mkChar("sqlite.matrix"));
    SET_STRING_ELT(tmp, 1, mkChar("matrix"));
    SET_CLASS(ret, tmp);

    /* set smat dim */
    PROTECT(tmp = NEW_INTEGER(2)); i++;
    INTEGER(tmp)[0] = nrows;
    INTEGER(tmp)[1] = ncols;
    /*_installAttrib(ret, R_DimSymbol, tmp);*/
    setAttrib(ret, R_DimSymbol, tmp);

    /* set smat dimname */
    PROTECT(tmp = NEW_LIST(2)); i++;
    SET_VECTOR_ELT(tmp, 0, _get_rownames2(iname));
    SET_VECTOR_ELT(tmp, 1, names);
    /*_installAttrib(ret, R_DimNamesSymbol, tmp);*/
    SET_DIMNAMES(ret, tmp);

    UNPROTECT(i);

    UNUSE_SDF2(iname);
    UNUSE_SDF2(mat_iname);
    return ret;
}
