#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ctype.h>
#include <R.h>  /* includes Rconfig.h */
#include <Rversion.h>
#include <Rdefines.h> /* Rinternals.h + GET_SLOT etc */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#endif

#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
#ifdef __sun
# include <alloca.h>
#endif
#endif

#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

SEXP triangularMatrix_validate(SEXP obj);
SEXP symmetricMatrix_validate(SEXP obj);
SEXP dense_nonpacked_validate(SEXP obj);

/* enum constants from cblas.h and some short forms */
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
#define RMJ CblasRowMajor
#define CMJ CblasColMajor
#define NTR CblasNoTrans
#define TRN CblasTrans
#define CTR CblasConjTrans
#define UPP CblasUpper
#define LOW CblasLower
#define NUN CblasNonUnit
#define UNT CblasUnit
#define LFT CblasLeft
#define RGT CblasRight

#if !defined(R_VERSION) || R_VERSION < R_Version(2, 7, 0)
char La_norm_type(const char *typstr);
char La_rcond_type(const char *typstr);
#endif

double get_double_by_name(SEXP obj, char *nm);
SEXP set_double_by_name(SEXP obj, double val, char *nm);
SEXP as_det_obj(double val, int log, int sign);
SEXP get_factors(SEXP obj, char *nm);
SEXP set_factors(SEXP obj, SEXP val, char *nm);

#if 0
SEXP dgCMatrix_set_Dim(SEXP x, int nrow);
#endif	/* unused */

/* int csc_unsorted_columns(int ncol, const int p[], const int i[]); */
/* void csc_sort_columns(int ncol, const int p[], int i[], double x[]); */
/* SEXP csc_check_column_sorting(SEXP A); */
SEXP Matrix_make_named(int TYP, char **names);
SEXP check_scalar_string(SEXP sP, char *vals, char *nm);
Rboolean equal_string_vectors(SEXP s1, SEXP s2);

void d_packed_getDiag(double *dest, SEXP x, int n);
void l_packed_getDiag(   int *dest, SEXP x, int n);
void tr_d_packed_getDiag(double *dest, SEXP x);
void tr_l_packed_getDiag(   int *dest, SEXP x);

SEXP Matrix_getElement(SEXP list, char *nm);

#define PACKED_TO_FULL(TYPE)						\
TYPE *packed_to_full_ ## TYPE(TYPE *dest, const TYPE *src,		\
			     int n, enum CBLAS_UPLO uplo)
PACKED_TO_FULL(double);
PACKED_TO_FULL(int);
#undef PACKED_TO_FULL

#define FULL_TO_PACKED(TYPE)						\
TYPE *full_to_packed_ ## TYPE(TYPE *dest, const TYPE *src, int n,	\
			      enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)
FULL_TO_PACKED(double);
FULL_TO_PACKED(int);
#undef FULL_TO_PACKED


extern	 /* stored pointers to symbols initialized in R_init_Matrix */
#include "Syms.h"

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/* number of elements in one triangle of a square matrix of order n */
#define PACKED_LENGTH(n)   ((n) * ((n) + 1))/2

/* duplicate the slot with name given by sym from src to dest */

#define slot_dup(dest, src, sym)  SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

/* is not yet used: */
#define slot_nonNull_dup(dest, src, sym)			\
    if(GET_SLOT(src, sym) != R_NilValue)			\
	SET_SLOT(dest, sym, duplicate(GET_SLOT(src, sym)))

/* TODO: Make this faster for the case where dimnames = list(NULL,NULL)
 *       and hence don't have to be set ! */
#define SET_DimNames(dest, src) slot_dup(dest, src, Matrix_DimNamesSym)


#define uplo_P(_x_) CHAR(STRING_ELT(GET_SLOT(_x_, Matrix_uploSym), 0))
#define diag_P(_x_) CHAR(STRING_ELT(GET_SLOT(_x_, Matrix_diagSym), 0))
#define class_P(_x_) CHAR(asChar(getAttrib(_x_, R_ClassSymbol)))

/* should also work for "matrix" matrices: */
#define Real_KIND(_x_)	(IS_S4_OBJECT(_x_) ? Real_kind(_x_) : \
			 (isReal(_x_) ? 0 : (isLogical(_x_) ? 1 : -1)))
/* This one gives '0' also for integer "matrix" :*/
#define Real_KIND2(_x_)	(IS_S4_OBJECT(_x_) ? Real_kind(_x_) : \
			 (isLogical(_x_) ? 1 : 0))

/* requires 'x' slot: */
#define Real_kind(_x_)	(isReal(GET_SLOT(_x_, Matrix_xSym)) ? 0	:	\
			 (isLogical(GET_SLOT(_x_, Matrix_xSym)) ? 1 : -1))

#define DECLARE_AND_GET_X_SLOT(__C_TYPE, __SEXP)	\
    __C_TYPE *xx = __SEXP(GET_SLOT(x, Matrix_xSym))


/**
 * Check for valid length of a packed triangular array and return the
 * corresponding number of columns
 *
 * @param len length of a packed triangular array
 *
 * @return number of columns
 */
static R_INLINE
int packed_ncol(int len)
{
    int disc = 8 * len + 1;	/* discriminant */
    int sqrtd = (int) sqrt((double) disc);

    if (len < 0 || disc != sqrtd * sqrtd)
	error(_("invalid 'len' = %d in packed_ncol"));
    return (sqrtd - 1)/2;
}

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * NOTE:  GET_SLOT(x, what)        :== R_do_slot       (x, what)
 * ----   SET_SLOT(x, what, value) :== R_do_slot_assign(x, what, value)
 * and the R_do_slot* are in src/main/attrib.c
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

/**
 * Expand compressed pointers in the array mp into a full set of indices
 * in the array mj.
 *
 * @param ncol number of columns (or rows)
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] to hold the result
 *
 * @return mj
 */
static R_INLINE
int* expand_cmprPt(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}

void make_d_matrix_triangular(double *x, SEXP from);
void make_i_matrix_triangular(   int *x, SEXP from);

void make_d_matrix_symmetric(double *to, SEXP from);
void make_i_matrix_symmetric(   int *to, SEXP from);

SEXP Matrix_expand_pointers(SEXP pP);

SEXP dup_mMatrix_as_dgeMatrix(SEXP A);
SEXP dup_mMatrix_as_geMatrix (SEXP A);

SEXP new_dgeMatrix(int nrow, int ncol);

static R_INLINE SEXP
mMatrix_as_dgeMatrix(SEXP A)
{
    return strcmp(class_P(A), "dgeMatrix") ? dup_mMatrix_as_dgeMatrix(A) : A;
}

static R_INLINE SEXP
mMatrix_as_geMatrix(SEXP A)
{
    return strcmp(class_P(A) + 1, "geMatrix") ? dup_mMatrix_as_geMatrix(A) : A;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 *
 * @param class string to match
 * @param valid vector of possible matches terminated by an empty string
 *
 * @return index of match or -1 for no match
 */
static R_INLINE int
Matrix_check_class(const char *class, char **valid)
{
    int ans;
    for (ans = 0; ; ans++) {
	if (!strlen(valid[ans])) return -1;
	if (!strcmp(class, valid[ans])) return ans;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_MUTILS_H_ */
