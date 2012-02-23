//
// Wrapper for the LBGFSB fortran (eugh) library
//
#include <stdio.h>
#include <string.h>
#include "mex.h"

const int LBFGSB_TASK_SIZE = 60;

mxArray *rhs[2];

const int    defaultm          = 5;
const int    defaultmaxiter    = 20;
const double defaultfactr      = 1e7;
const double defaultpgtol      = 1e-5;
const int    defaultprintlevel = -1;

enum lbfgsb_res_t {
	CONVERGED, 
	ABNORMAL_TERMINATION,
	ERROR,
	MAX_ITER,
	UNHANDLED
};

// objective & gradient evaluation function
typedef double (*lbfgsb_eval_t)(const int N, const double *x, double *g);
// iteration notification callback
typedef void (*lbfgsb_iter_t)(const int t, const double f, 
		const double *g);

/* optional params */
typedef struct {
	lbfgsb_iter_t iter;
	int maxiter;
	int m;
	double pgtol;
	double factr;
} lbfgsb_options_t;

// This is the L-BFGS-B routine implemented in Fortran 77.
extern "C" void setulb_ (int* n, int* m, double x[], double l[], 
			 double u[], int nbd[], double* f, double g[], 
			 double* factr, double* pgtol, double wa[], 
			 int iwa[], char task[], int* iprint, 
			 char csave[], bool lsave[], int isave[], 
			 double dsave[]);

class LBFGSB {
private:
	lbfgsb_eval_t eval_func;
	lbfgsb_iter_t iter_func;
	int     n;       // The number of variables.
	double* x;       // The current point.
	double* lb;      // The lower bounds.
	double* ub;      // The upper bounds.
	int*    btype;   // The bound types.
	double  f;       // The value of the objective.
	double* g;       // The value of the gradient.
	int     iprint;  // The print level.
	int     maxiter; // The maximum number of iterations.
	double  pgtol;   // Convergence parameter passed to L-BFGS-B.
	double  factr;   // Convergence parameter passed to L-BFGS-B.
	int     m;       // The number of variable corrections to the 
					 // limited-memory approximation to the Hessian.

	// These are structures used by the L-BFGS-B routine.
	double* wa;
	int*    iwa;
	char    task[LBFGSB_TASK_SIZE+1];
	char    csave[60];
	bool    lsave[4];
	int     isave[44];
	double  dsave[29];
public:
	void init(int n, double *x, double *lb, double *ub, 
			lbfgsb_eval_t eval, lbfgsb_options_t *op);
	LBFGSB(int n, double *x, double *lb, double *ub, lbfgsb_eval_t eval);
	LBFGSB(int n, double *x, double *lb, double *ub, lbfgsb_eval_t eval, 
			lbfgsb_options_t *op);
	~LBFGSB(void);
	void call(const char* cmd  = NULL);
	bool isTask(const char* cstr);
	lbfgsb_res_t solve(double& );
};

void LBFGSB::init(int n, double *x, double *lb, double *ub, 
		lbfgsb_eval_t eval, lbfgsb_options_t *op) {
	this->n       = n;
	this->x       = x;
	this->lb      = lb;
	this->ub      = ub;
	this->btype   = new int[n];
	this->m       = op->m;
	this->maxiter = op->maxiter;
	this->factr   = op->factr;
	this->pgtol   = op->pgtol;
	this->iprint  = defaultprintlevel;
	
	for (int i = 0; i < LBFGSB_TASK_SIZE; i++)
		task[i] = ' ';
	
	task[LBFGSB_TASK_SIZE] = '\0';
	/* set evaluation function */
	this->eval_func = eval;
	this->iter_func = op->iter;
	
	/* set btype to 2 = both upper and lower bounds */
	for(int i=0; i < n; i++)
		this->btype[i] = 2;
	
	f   = 0;
	g   = new double[n];
	wa  = new double[(2*m*n + 11*m*m + 5*n + 8*m)];
	iwa = new int[3*n];
}

LBFGSB::LBFGSB(int n, double *x, double *lb, double *ub, 
		lbfgsb_eval_t eval) {
	lbfgsb_options_t defaults = {NULL, defaultmaxiter, defaultm, 
		defaultpgtol, defaultfactr};
	init(n, x, lb, ub, eval, &defaults);
}

LBFGSB::LBFGSB(int n, double *x, double *lb, double *ub, 
		lbfgsb_eval_t eval, lbfgsb_options_t *op) {
	init(n, x, lb, ub, eval, op);
}

LBFGSB::~LBFGSB(void) {
	delete [] this->btype;
	delete [] this->g;
	delete [] this->wa;
	delete [] this->iwa;
}

/* call the fortran library */
void LBFGSB::call(const char* cmd) {
	if(cmd) {
		// Get the length of the source C string.
		int nsource = strlen(cmd);
		// Only perform the copy if the source can fit into the destination.
		if (nsource < 60) {
			for(int i = 0; i < nsource; i++)
				task[i] = cmd[i];

			// Fill in the rest of the string with blanks.
			for (int i = nsource; i < 60; i++)
				task[i] = ' ';
		}
	}
	setulb_(&n,&m,x,lb,ub,btype,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,
		csave,lsave,isave,dsave);
}

bool LBFGSB::isTask(const char* cstr) {
  return !strncmp(task,cstr,strlen(cstr));
}

lbfgsb_res_t LBFGSB::solve(double& fc) {
	lbfgsb_res_t res = UNHANDLED;
	// Initialize the objective function and gradient to zero.
	f = 0;
	for (int i = 0; i < n; i++)
		g[i] = 0;
	
	call("START");
	// Repeat until we've reached the maximum number of iterations.
	int t = 0;
	while (true) {
		// Do something according to the "task" from the previous call to
		// L-BFGS.
		if (isTask("FG")) {
			// Evaluate the objective function and the gradient of the
			// objective at the current point.
			f = eval_func(n, x, g);
		} else if (isTask("NEW_X")) {
			t++;
			if(iter_func != NULL)
				iter_func(t, f, g);
			
			if (t == maxiter) {
				call("STOP");
				res = MAX_ITER;
				break;
			}
		} else if (isTask("CONV")) {
			res = CONVERGED;
			break;
		} else if (isTask("ABNO")) {
			//mexPrintf("Task: %s\n", task);
			res = ABNORMAL_TERMINATION;
			break;
		} else if (isTask("ERROR")) {
			//mexPrintf("Task: %s\n", task);
			res = ERROR;
			break;
		}
		// Call L-BFGS again.
		call();
	}
	fc = f;
	return res;
}

static double evaluate(const int n, const double *x, double *g) {
	double fx = 0.0;

    mxArray *lhs[2];
    memcpy(mxGetPr(rhs[1]),x,n*sizeof(double));

    // Call matlab function and return function value and gradient
    // rhs contains: (function handle,x)
    // lhs contains: (function value,gradient)
    if (mexCallMATLAB(2,lhs,2,rhs,"feval"))
        mexErrMsgTxt("Error while evaluating objective function.");        
    
    // Retrieve function value and gradient
    memcpy(g,mxGetPr(lhs[1]),n*sizeof(double));
    fx = mxGetScalar(lhs[0]); 
    mxDestroyArray(lhs[0]); mxDestroyArray(lhs[1]);
    
    return fx;
}

static void iter(const int t, const double f, const double *g) {
	mexPrintf("Iter: %d F = %2.9f\n", t, f);
}

void mexFunction(int nlhs, mxArray *plhs[], 
		int nrhs, const mxArray *prhs[]) {
    mwSize N;
    double *x, *lb, *ub, fx;
	lbfgsb_options_t op;
    
    if (nrhs<4)
        mexErrMsgTxt("lbfgs_ requires exactly 4 input arguments.");
    if( mxGetClassID(prhs[0]) != mxFUNCTION_CLASS )    
        mexErrMsgTxt("The first argument must be a function handle.");
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("The second argument must be a column vector.");
    if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("The third argument must be a column vector.");
    if (!mxIsDouble(prhs[3]))
        mexErrMsgTxt("The fourth argument must be a column vector.");
	
	if (!mxIsStruct(prhs[4]))
        mexErrMsgTxt("The fith argument must be a struct.");     
    
	N = mxGetM(prhs[1])*mxGetN(prhs[1]);
	
    // Get function handle for use with mexCallMATLAB
    rhs[0] = mxDuplicateArray(prhs[0]);    
    rhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
	
	op.m = (int) mxGetScalar(mxGetField(prhs[4],0,"m"));
	op.maxiter = (int) mxGetScalar(mxGetField(prhs[4],0,"maxiter"));
	op.pgtol = (double) mxGetScalar(mxGetField(prhs[4],0,"pgtol"));
	op.factr = (double) mxGetScalar(mxGetField(prhs[4],0,"factr"));
	op.iter = iter;
    
    // Allocate memory for solution vector and initialize it to x0
	// do the same for the lower and upper bounds	
    x = (double*)mxMalloc(N*sizeof(double));
	lb = (double*)mxMalloc(N*sizeof(double));
	ub = (double*)mxMalloc(N*sizeof(double));
    memcpy(x, mxGetPr(prhs[1]),N*sizeof(double));
	memcpy(lb,mxGetPr(prhs[2]),N*sizeof(double));
	memcpy(ub,mxGetPr(prhs[3]),N*sizeof(double));

	LBFGSB lfbsgb(N, x, lb, ub, evaluate, &op);
    
	lbfgsb_res_t res = lfbsgb.solve(fx);
	mxArray *res_str;
	switch(res) {
		case CONVERGED:
			mexPrintf("L-BFGS-B converged\n");
			res_str = mxCreateString("converged");
			break;
		case ABNORMAL_TERMINATION:
			mexPrintf("L-BFGS-B abnormal termination\n");
			res_str = mxCreateString("abnormal termination");
			break;
		case ERROR:
			mexPrintf("L-BFGS-B error\n");
			res_str = mxCreateString("error");
			break;
		case MAX_ITER:
			mexPrintf("L-BFGS-B max iterations reached\n");
			res_str = mxCreateString("max iterations reached");
			break;
		case UNHANDLED:
		default:
			mexWarnMsgTxt("Unknown result state");
			res_str = mxCreateString("Unknown result state");
			break;
	}
    
    // Allocate outputs
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
	if(nlhs > 1)
		plhs[1] = mxCreateDoubleScalar(fx);
	
    if(nlhs > 2)
		plhs[2] = mxCreateDoubleScalar(res);
    // Copy current iterate to output
    memcpy(mxGetPr(plhs[0]),x,N*sizeof(double));
}