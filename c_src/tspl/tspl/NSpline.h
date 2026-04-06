// NSpline.h: interface for the NSpline class.
//
//////////////////////////////////////////////////////////////////////
/*-------------------------
Brief description (G.Beliakov, 10.9.01)

  This code is based on my papers on least squares constrained splines

  There are 3 main classes

  LSNSpline  - least squares univariate constrained spline
			   handles monotonicity and convexity constraints
  LSNSpline2 - least squares bi-linear constrained spline
               tensor product. handles monotonicity constraints
  LSNSplineB - least squares simple B-spline (no constraints)

  NSpline is the working horse: this class represents N-spline of arbitrary degree
  on any sequence of knots (could be multiple, in this case deficient spline)

  LSNSpline class builds constrained LS spline. 
  Example of usage:
  1) define spline knots t, data to be fitted x,y and data to be interpolated xe,ye
  2) declare class LSNSpline S;
  3) build spline S.Build(kind, order, t,x,y,xe,ye); or
	S.Build(kind, order,t,x,y) if no interpolation conditions.
	kind is interpreted as below in class declaration
  4) if returned value is 0 (success) compute spline value and 
    its derivative at any point p, S.Value(p) or S.ValueDer(p, 2)

	LSNSplineB is a smaller and faster version. Does not use constraints
	Same calling sequence, but kind variable is ignored. 
	Same result as LSNSpline, as if kind=0, but faster

  LSNSpline2 class builds bi-linear constrained spline (can do other orders, 
  but not yet implemented, except when unconstrained). It's a tensor product
  spline, so two LSNSpline classes are used internally.

  Calling sequence is similar, but now S.Build must specify
  2 spline orders n1,n2, two knots sequences t1,t2, and data abscisae are
  matrices (dMat) and not vectors. Example
  S.Build(kind, orderX, orderY, tX, tY, X,y, XE,ye)

  Similarly, call S.Value(pX, pY) calculates spline value and
  S.ValueDer(pX, pY, VAR) calculates 1 partial derivative wtr to VAR (0 or 1)



-------------------------*/
#ifndef LSNSPL_H
#define LSNSPL_H
// for faster operation uncomment this line
#define TNT_NO_BOUNDS_CHECK

//#include "tnt/tnt.h" // this is a NIST set of templates for matrix-vector
//#include "tnt/vec.h"
//#include "tnt/cmat.h"

#include <vector>

#include "math.h"
#include "lsei.h"


#include "matrix2d.cpp"

using namespace std; // We need a couple of classes from here (lists)

typedef vector<double>	dVec;
typedef vector<int>		iVec;
typedef Array2D<double>		dMat;

#define Infinity 1.0E12

// return values
#define INVALID_CODE 10

// used in B-splines to calculate the external knots, if we want them to be multiple
// has litte effect on spline though
#define DELTA_T 0.00001

#include "Tensor.h" // my small service class to handle tensor products


// this class is the working horse to evaluate B-splines
class NSpline  
{
public:
	dVec		knots;			// knot sequence
	dVec		auxknots;		// simple knots (with no multiplicities) 
	int		order;				// order of the spline
	long	N, Naux;			// size of knots array
	int		deficiency;			// calculated. Total deficiency of spline. 1 for normal splines	
	int		invalid;			// cannot be evaluated before initialised
private:
	int		total_deficiency;

public:
	int	Init(int M, int n);		// constructor. Creates uniform mesh of knots from 0 to n
	int Init(int, dVec&);		// constructor. deficiency is calculated from knots

	void Value(double, dVec&, iVec&, int M ); // calculates the value of all nonzero B-splines of order M
	void Value(double, dVec&, iVec&); // calculates the value of all nonzero B-splines
//	at the point t. dVec is output array of size (M+1)*deficiency. Contains all nonzero splines
//  of order M, and then (next M+1 elements) of order M-1, etc. (for deficient splines)
//  iVec contains pairs of integers indicating the indices of the first and last nonzero B splines 
//  of order M, M-1,...M-def, these correspond to the indices of spline coefficients the splines need 
//  to be multiplied. I.e. iVec[0] and iVec[1] contain l and r, such that Bl,Bl+1,...Br are nonzero
//  and the values are given in dVec[0],[1],...
//  iVec[2] and iVec[3] again contain l and r for B-splines of order M-1, etc. till iVec[(def-1)*2],iVec[(def-1)*2+1]
//  to find deficiency one can use "deficiency" member variable or calculate it from size of iVec

private:
	int previousL;		// accelerates search for the interval-the interval previous point was found

	void FindLp(double t, int &l, dVec *Knots);
		// finds the interval for t. checks result on the previous iteration. Faster than general
	void FindL(double t, int& l, dVec* Knots,  int ri, int le=0);
		// finds the interval for t. General subdivision method
public:
	virtual double Mknots(int i);		// use instead of knots - returns knots # ...-2,-1,0,1,...N,N+1,N+2
//private:
	virtual double Mknotsaux(int i);		// use instead of knots - returns knots # ...-2,-1,0,1,...N,N+1,N+2
	virtual int GetSize() {return N;};
};


class LSNSpline: public NSpline {
friend class LSNSpline2; // bi-linear spline uses 2 LSNSplines
friend class LSNSplineB;
friend class LSNSplineT;
public:
	LSNSpline() { der_boundlow.resize(0); der_boundup.resize(0); 
	kind=0; equality=0; der_boundlow_scalar=der_boundup_scalar=0;
	der_boundlow_scalar_specified=der_boundup_scalar_specified=0;};

	double		Value(double t);		// value of the spline
	double		ValueDer(double t,  int Nders); // value of the spline and derivatives
	// value of the spline and N-th derivative at point t

	// LS + equality conditions
	int	Build(int kind, int n, dVec& t, dVec& x,  dVec& y, dVec& xe,  dVec& ye);
// this is the main routine to compute spline coefs.
// kind as described below, n - spline order 
// t- sequence of knots, x - data abscisae, y - data values
// xe, ye - data to be interpolated
// no data is affected in any way

	// normal LS as above, with no equality conds
	int	Build(int kind, int n, dVec& t,  dVec& x, dVec& y);

/*   Interpretation of Kind
	0 = LS spline, no constraints 
	1 = monotone increasing, 2 = monotone decreasing
	3 = monotone increasing, strict, 4= monotone decreasing strict
	5 = bounded derivative
	6 = concave up, 7= concave down
	8 = concave up strict, 9 = concave down strict (bounds are specified)
	10 = second derivative bounded
*/

	// assigns bounds for derivatives. Service methods
	// parameters can be scalars (bounds on the whole interval) or vectors
	// (of size N) - nonuniform bounds
	void		Bounds(dVec& td, dVec& td1);
	void		Bounds(double td, double td1);
	void		Boundlow(dVec& td);
	void		Boundlow(double td);
	void		Boundup(dVec& td);
	void		Boundup(double td);

	virtual int GetSize() {return Nknots;};

private:
	long		size;      // size of the vector of coefs
	dVec		coefs;     // c[-m+1],...,c[N-1]
	dMat		coefs_der;  // coeficients of the derivatives

	dVec		D[2]; // diagonal matrix Dk
	
	int			equality;  // 1 if interpolation conds are specified
	int			kind;		// see above

	long		rowsG, rowsEQ, rowsA;	// size of the matrix

// Our working horse. Populates the matrix of the system and solves it
// by LSEI or WNNLS methods Returns error code (MODE), 0 if OK, 2,3 if system inconsistent
	virtual int	FindCoefficients(int kind, dVec& x, dVec& y, dVec& xe, dVec& ye);

// aux routine for bi-linear spline
	dVec		RowOfValues(double t);
//private:
	int			startindex,Nknots;
	dVec		Bsplines;
	iVec		Index;	
	// auxiliary stuff for LSEI
	dVec		progopts;
	lsei		m_LSEI;

private:
	dVec		der_boundlow, der_boundup; // bounds on derivative
	double		der_boundlow_scalar,der_boundup_scalar;
	int			der_boundlow_scalar_specified,der_boundup_scalar_specified;

	// computation and multiplication by a diagonal matrix D
	void multiplyD(dMat& equations);
	void computeD(int locorder);
	void computeD();
};

class LSNSplineB: public LSNSpline
// this class is derived from LSNSpline. It uses B-splines representation 
// and no constraints
{
public:
	int	FindCoefficients(int kind, dVec& x, dVec& y, dVec& xe, dVec& ye);
	virtual double Mknots(int i);		// use instead of knots - returns knots # ...-2,-1,0,1,...N,N+1,N+2
private:
	virtual double Mknotsaux(int i);		// use instead of knots - returns knots # ...-2,-1,0,1,...N,N+1,N+2
};

class LSNSpline2: public NSpline {
// two-dimensional monotone spline of order 2 

public:
	LSNSpline2() {der_boundup1=der_boundup2=der_bounddown1=der_bounddown2=0.0;}

	LSNSpline	spl1, spl2;
	double		Value(double t1, double t2);		// value of the spline
	double		ValueDer(double t1,  double t2, int var); 
	// value of the spline and 1st derivative with respect to var

	int Build(int kind, int n1, int n2, dVec& t1,  dVec& t2, 
					   dMat& x,  dVec& y, dMat& xe, dVec& ye); 
// this is the main routine to compute spline coefs.
// kind as described below, n1,n2 - spline order (set to 2 if constrained)
// t1,t2 - sequence of knots, x- data abscisae, y data values
// xe, ye - data to be interpolated
// no data is affected in any way
	int Build(int kind, int n1, int n2, dVec& t1, dVec& t2, dMat& x,  dVec& y); 
	// a variant of above with no interpolation conditions

	int			kind;  // kind of spline
/*	0 = LS spline,  no constraints. order may be >2
	otherwise order is forced to be 2

	interpretation of kind: 
	1*   - will be first var code + 10* will be the second var. code
	code means:
	1 = monotone LS, 2 = monotone decreasing
	3 = monotone increasing, strict, 4= monotone decreasing strict
	5 = bounded derivative
	in cases 3,4,5 the bounds on the derivative must be specified (0 default value)
	Examples:
	01 - unconstrained in y, monotone in x
	21 - decreasing in y increasing in x
	33 - strictly increasing in both vars
	55 - bounded derivatives in both vars
	50, 15 35, etc are all reduced to 55, and bounds are set automatically (+/- Infinity)
*/

// service methods to set bounds on derivatives
	void	BoundUpX(double d) {der_boundup1=d;}
	void	BoundDownX(double d) {der_bounddown1=d;}
	void	BoundsX(double ddown, double dup) {der_bounddown1=ddown; der_boundup1=dup;}
	void	BoundUpY(double d) {der_boundup2=d;}
	void	BoundDownY(double d) {der_bounddown2=d;}
	void	BoundsY(double ddown, double dup) {der_bounddown2=ddown; der_boundup2=dup;}

private:
	dMat		coefs; // here matrix of spline coefs is stored
	dMat		coefs_der[2]; // here coefs of derivatives

	long		size, rowsG, rowsEQ, rowsA;	// size of the matrix
	int			invalid; // is 0 after coefs are computed, 1 otherwise

	double		der_boundup1,der_bounddown1;
	double		der_boundup2,der_bounddown2;

// this is our working horse method. Calculates the coefficients setting the
// entries of the matrix and solving LSEI problem. Returns error code (INVALID_CODE)
// or MODE output parameter from LSEI or WNNLS routines (0 if OK, 2,3 if inconsistent)
	int	FindCoefficients(int kind, dMat& x, dVec& y, dMat& xe, dVec& ye);
};

#endif
