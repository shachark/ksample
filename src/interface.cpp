#include "Dynslicing.h"

#undef ERROR
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

using namespace std;

extern "C" {

// NOTE:
// In all of the below, R_y is assumed to be sorted in ascending order and
// contain quantiles of F0 at the raw observations.
// R_x is assumed to hold sample identifiers in [0, ..., K-1] ordered according
// to ascending pooled-sample observation order.

SEXP R_C_dynslicing_1eqp(SEXP R_y, SEXP R_lambda) {
	double* y = REAL(R_y);
	int len = length(R_y);
	double lambda = *REAL(R_lambda);

	SEXP R_output;
	PROTECT(R_output = allocVector(REALSXP, 1));
	double* res = REAL(R_output);

	*res = dynslicing_1eqp(y, len, lambda);

	UNPROTECT(1);
	return (R_output);
}

SEXP R_C_dynslicing_one(SEXP R_y, SEXP R_lambda, SEXP R_alpha) {
	double* y = REAL(R_y);
	int len = length(R_y);
	double lambda = *REAL(R_lambda);
	double alpha = *REAL(R_alpha);

	SEXP R_output;
	PROTECT(R_output = allocVector(REALSXP, 1));
	double* res = REAL(R_output);

	*res = dynslicing_one(y, len, lambda, alpha);

	UNPROTECT(1);
	return (R_output);
}

SEXP R_C_dynslicing_keqp(SEXP R_x, SEXP R_dim, SEXP R_lambda, SEXP R_slices_wanted) {
	int* x = INTEGER(R_x);
	int len = length(R_x);
	int dim = *INTEGER(R_dim);
	double lambda = *REAL(R_lambda);
	int slices_wanted = *INTEGER(R_slices_wanted);

	SEXP R_output;
	double stat;

	if (slices_wanted) {
		vector<int> slices;
		stat = dynslicing_keqp(x, len, dim, lambda, slices);
		int nr_slices = slices.size();

		PROTECT(R_output = allocVector(REALSXP, 1 + nr_slices));
		double* res = REAL(R_output);

		res[0] = stat;
		for (int i = 0; i < nr_slices; ++i) {
			res[i + 1] = slices[i];
		}
	} else {
		PROTECT(R_output = allocVector(REALSXP, 1));
		double* res = REAL(R_output);

		*res = dynslicing_keqp(x, len, dim, lambda);
	}

	UNPROTECT(1);
	return (R_output);
}

SEXP R_C_dynslicing_k(SEXP R_x, SEXP R_dim, SEXP R_lambda, SEXP R_slices_wanted) {
	int* x = INTEGER(R_x);
	int len = length(R_x);
	int dim = *INTEGER(R_dim);
	double lambda = *REAL(R_lambda);
	int slices_wanted = *INTEGER(R_slices_wanted);

	SEXP R_output;
	double stat;

	//Rprintf("len = %d, dim = %d, lambda = %g, x[0] = %d, x[1] = %d\n", len, dim, lambda, x[0], x[1]);

	if (slices_wanted) {
		vector<int> slices;
		stat = dynslicing_k(x, len, dim, lambda, slices);
		int nr_slices = slices.size();

		PROTECT(R_output = allocVector(REALSXP, 1 + nr_slices));
		double* res = REAL(R_output);

		res[0] = stat;
		for (int i = 0; i < nr_slices; ++i) {
			res[i + 1] = slices[i];
		}
	} else {
		PROTECT(R_output = allocVector(REALSXP, 1));
		double* res = REAL(R_output);

		*res = dynslicing_k(x, len, dim, lambda);
	}

	UNPROTECT(1);
	return (R_output);
}

}
