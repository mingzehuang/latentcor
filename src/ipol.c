#include "ipol.h"
#define UNUSED(x) (void)(x)


// return index of the smallest element larger or equal to val
// binary search is slower than linear for small n, but we don't optimize that here
// small grids anyway run faster
// we never return 0, if val < arr[0] we return 1
// if val > arr[n-1], we return n-1
// this is suited for our particular application
static inline int binsearch(const int n,double *arr, const double val) {
  int toosmall, toolarge, test;
  toosmall = 0;
  toolarge = n-1;
  while(toosmall+1 < toolarge) {
    test = (toosmall + toolarge) >> 1;
    if(arr[test] >= val) {
      toolarge = test;
    } else if(arr[test] < val) {
      toosmall = test;
    }
  }
  return toolarge;
}

static void C_predmlip(const int rank, double **grid, int *dims, double *values, double *output) {
  int stride[rank+1];
  stride[0] = 1;
  for(int i = 1; i <= rank; ++i) stride[i] = dims[i-1]*stride[i-1];
  int N = stride[rank];
  int dim[rank];
  for(int i = 0; i < rank; ++i) {
    dim[i] = 0;
  }
  /* Loop over entire grid. For each point, predict it as the mean of all its immediate neighbours. */
  for(int j = 0; j < N; j++,dim[0]++) {
    for(int i = 0; i < rank-1; i++) {
      if(dim[i] >= dims[i]) {
	dim[i] -= dims[i];
	dim[i+1]++;
      } else
	break;
    }
    //    printf("j %d: ",j);for(int i = 0; i < rank; i++) printf("%d ",dim[i]); printf("\n");
    double pred = 0;
    int w = 0;
    for(int i = 0; i < rank; i++) {
      if(dim[i] <= 0 || dim[i]+1 >= dims[i]) continue;
      double *gr = grid[i];
      double r = (gr[dim[i]] - gr[dim[i]-1])/(gr[dim[i]+1] - gr[dim[i]-1]);
      w++;
      pred += r*values[j+stride[i]] + (1-r)*values[j-stride[i]];
    }
    if(w == 0)
      output[j] = values[j];
    else
      output[j] = pred/w;
  }
}

static SEXP R_mlippred(SEXP sgrid, SEXP values) {
  int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  double *grid[rank];
  SEXP resvec;
  if(!IS_NUMERIC(values)) error("values must be numeric");

  /* Make some pointers into the grid data */
  for(int i = 0; i < rank; i++) {
    dims[i] = LENGTH(VECTOR_ELT(sgrid,i));
    grid[i] = REAL(VECTOR_ELT(sgrid,i));
    gridsize *= dims[i];
  }

  if(LENGTH(values) != gridsize) error("grid has size %d, you supplied %d values",
				       gridsize,LENGTH(values));

  PROTECT(resvec = NEW_NUMERIC(gridsize));
  (void) C_predmlip(rank,grid,dims,REAL(values),REAL(resvec));
  UNPROTECT(1);
  return resvec;
}

static double C_evalmlip(const int rank, double *x, double **grid, int *dims,
			 double *values, int blend) {

  double weight[rank];
  int valpos = 0;
  double ipval = 0;
  /*
   Find first grid point which is larger or equal, compute
   the weight, store it. As well as the linear position valpos in the grid
   A binary search is faster if there are more grid points
   but up to 10-20 linear search is usually equally fast or faster due to simplicity
   Note that we never check whether x is within the grid. If it's outside we assume
   the function continues linearly outwards.
   Remaining optimizations:  We could have precomputed the stuff we divide by, if we should
   do more than one vector in a row, this would save some time.  Likewise the stride (which is
   integer multiplication, which is slow.)
  */

  int stride = 1;
  for(int i = 0; i < rank; i++) {
    double *gdvec = grid[i];
    // We could use Rs findInterval, but how many ways are there to do this?
    int gp = binsearch(dims[i],gdvec,x[i]);
    weight[i] = (x[i]-gdvec[gp-1])/(gdvec[gp]-gdvec[gp-1]);
    valpos += stride*gp;
    stride *= dims[i];
  }

  // loop over the corners of the box, sum values with weights
  double wsum = 0.0;
  for(int i = 0; i < (1<<rank); i++) {
    // i represents a corner. bit=1 if upper corner, 0 if lower corner.
    // We should find its weight
    // it's a product of the weights from each dimension
    int vpos = valpos;
    int stride = 1;
    double cw = 1;
    for(int g = 0; g < rank; g++) {
      if( (1<<g) & i) {
	cw *= blendfun(weight[g],blend);
      } else {
	cw *= blendfun(1-weight[g],blend);
	vpos -= stride;
      }
      stride *= dims[g];
    }
    wsum += cw;
    ipval += cw*values[vpos];
  }
  return ipval/wsum;
}

/* Then a multilinear approximation */
static SEXP R_evalmlip(SEXP sgrid, SEXP values, SEXP x, SEXP Rthreads, SEXP Sblend) {
  const int rank = LENGTH(sgrid);
  int gridsize = 1;
  int dims[rank];
  int threads = INTEGER(AS_INTEGER(Rthreads))[0];
  double *grid[rank];
  int blend = 0;
  if(!IS_NUMERIC(values)) error("values must be numeric");
  if(!IS_NUMERIC(x)) error("argument x must be numeric");
  if(isMatrix(x) ? (nrows(x) != rank) : (LENGTH(x) != rank))
    error("grid has dimension %d, you supplied a length %d vector",
	  rank, isMatrix(x) ? nrows(x) : LENGTH(x));

  /* Make some pointers into the grid data */
  for(int i = 0; i < rank; i++) {
    dims[i] = LENGTH(VECTOR_ELT(sgrid,i));
    grid[i] = REAL(VECTOR_ELT(sgrid,i));
    gridsize *= dims[i];
  }

  if(LENGTH(values) != gridsize) error("grid has size %d, you supplied %d values",
				       gridsize,LENGTH(values));
  if(!isNull(Sblend)) blend = INTEGER(AS_INTEGER(Sblend))[0];
  const int numvec = isMatrix(x) ? ncols(x) : 1;
  double *xp = REAL(x);
#ifdef RETMAT
  SEXP resvec = PROTECT(allocMatrix(REALSXP, numvec, 1));
#else
  SEXP resvec = PROTECT(NEW_NUMERIC(numvec));
#endif
  double *out = REAL(resvec);
#pragma omp parallel for num_threads(threads) schedule(static) if (numvec > 1 && threads > 1)
  for(int i = 0; i < numvec; i++) {
    out[i] = C_evalmlip(rank,xp+i*rank,grid,dims,REAL(values),blend);
  }
  UNPROTECT(1);
  return resvec;
}

static SEXP R_sqdiffs(SEXP x1, SEXP x2, SEXP Sthreads) {
  // each column in x1 should be subtracted from each column in x2,
  // the squared column sums should be returned.
  const int r1 = nrows(x1), c1 = ncols(x1), r2 = nrows(x2), c2 = ncols(x2);
  int N = c1*c2;
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  SEXP res = PROTECT(NEW_NUMERIC(N));
  double *dres = REAL(res);
  double *np = dres;
#pragma omp parallel for num_threads(threads) schedule(static) if (c1 > 1)
  for(int i = 0; i < c1; i++) {
    double *x1p = REAL(x1) + i*r1;
    for(int j = 0; j < c2; j++) {
      double *x2p = REAL(x2) + j*r2;
      double csum = 0.0;
      for(int k = 0; k < r1; k++) {
	csum += (x1p[k] - x2p[k])*(x1p[k] - x2p[k]);
      }
      np[j + i*c2] = csum;
    }
  }
  SEXP snewdim = PROTECT(NEW_INTEGER(2));
  INTEGER(snewdim)[0] = c1;
  INTEGER(snewdim)[1] = c2;
  setAttrib(res,R_DimSymbol,snewdim);
  UNPROTECT(2);
  return res;
}

// inplace to save memory, do repated multiplication instead of pow(). It's faster.
static SEXP R_phifunc(SEXP Sx, SEXP Sk, SEXP Sthreads) {
  double k = REAL(AS_NUMERIC(Sk))[0];
  double *x = REAL(Sx);
  int threads = INTEGER(AS_INTEGER(Sthreads))[0];
  double *y;
  SEXP res;
  if(XLENGTH(Sx) == 1 && x[0] == 0.0) return ScalarReal((k<0)?1:0);

  if(MAYBE_REFERENCED(Sx)) {
    res = PROTECT(NEW_NUMERIC(XLENGTH(Sx)));
    y = REAL(res);
  } else {
    res = Sx;
    y = x;
  }

  if(k < 0) {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
    for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) y[i] = exp(k*x[i]);
  } else {
    int ki = INTEGER(AS_INTEGER(Sk))[0];
    if(ki % 2 == 1) {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
      for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) {
	// it's the sqrt(x) to ki'th power
	if(x[i] <= 0.0) {y[i]=0.0;continue;}
	y[i] = R_pow_di(sqrt(x[i]), ki);
      }
    } else {
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
      for(R_xlen_t i = 0; i < XLENGTH(Sx); i++) {
	// it's sqrt(x) to ki'th power, multiplied by 0.5 log(x)
	if(x[i] <= 0.0) {y[i]=0.0;continue;}
	y[i] = 0.5*log(x[i]) * R_pow_di(sqrt(x[i]),ki);
      }
    }
  }
  if(y != x) UNPROTECT(1);
  return res;
}

R_CallMethodDef callMethods[] = {
  {"evalmlip", (DL_FUNC) &R_evalmlip, 5},
  {"predmlip", (DL_FUNC) &R_mlippred, 2},
  {"sqdiffs", (DL_FUNC) &R_sqdiffs, 3},
  {"phifunc", (DL_FUNC) &R_phifunc, 3},
  {NULL, NULL, 0}
};


void attribute_visible R_init_latentcor(DllInfo *info) {
  if(info != NULL) {}; // avoid warning about unused parameter
  /* register our routines */
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info, FALSE);
  R_RegisterCCallable("latentcor", "evalmlip", (DL_FUNC) C_evalmlip);
}
