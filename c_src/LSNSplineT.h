// LSNSplineT.h: interface for the LSNSplineT class.
//
//////////////////////////////////////////////////////////////////////
#ifndef LSNSPLT_H
#define LSNSPLT_H

#include "nspline.h"
#include "tensor.h"

typedef vector<LSNSpline*> OneSplines;
typedef vector<dVec*> VecArray;
typedef vector<dMat*> MatArray;
typedef vector<Tensor1*> TensorArray; // tensor1


void deleteVecArray(VecArray& _t);

extern void MultiplyByL(dVec& r);
extern dVec TensorProductRow(dVec& r1, dVec& r2);

class LSNSplineT:  public NSpline 
{
public:
	LSNSplineT();
	OneSplines	spl;

	double		Value(dVec t);		// value of the spline
	double		ValueDer(dVec t, int var); 
	// value of the spline and 1st derivative with respect to var

	int Build(int kind, int dim, VecArray& t, VecArray& data, VecArray& dataE);
// t[i] points to a vector of knots wrt i-th variable
// data[i] points to a vector of size dim+1, containing abscisae and the value
// of the i-th data point. Same for dataE (interpolation data)
// this is the main routine to compute spline coefs.
// kind as described below, dim- number of variables
// no data is affected in any way
	int Build(int kind, int dim, VecArray& t, VecArray& data);
	// a variant of above with no interpolation conditions

	int Build(int kind, int dim, dMat& t, dMat& data, dMat& dataE);

	int			kind;  // kind of spline
	int			dim;	// number of variables

/*	0 = LS spline,  no constraints. order may be >2
	otherwise order is forced to be 2

	interpretation of kind: 
	1*   - will be first var code + 10* will be the second var. code + 100* third, etc
	e.g. 1204 4 vars, 1 4th var, 2 third 0 second 4-first variable code
	code means:
	1 = monotone LS, 2 = monotone decreasing
not implemented
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
	void	BoundUpN(double d, int var) {der_boundup[var]=d;}
	void	BoundDownN(double d, int var) {der_bounddown[var]=d;}
	void	BoundsN(double ddown, double dup, int var) {der_bounddown[var]=ddown; der_boundup[var]=dup;}

private:
	Tensor1		coefs; // here tensor of spline coefs is stored
	TensorArray		coefs_der; // here array of coefs of derivatives

	long		size, rowsG, rowsEQ, rowsA;	// size of the matrix
	int			invalid; // is 0 after coefs are computed, 1 otherwise

	dVec		der_boundup, der_bounddown;

// this is our working horse method. Calculates the coefficients setting the
// entries of the matrix and solving LSEI problem. Returns error code (INVALID_CODE)
// or MODE output parameter from LSEI or WNNLS routines (0 if OK, 2,3 if inconsistent)
	int	FindCoefficients(int kind, int _dim, VecArray& t, VecArray& data, VecArray& dataE);

	double NestedSum(int level, iVec& index);
	double NestedSumDer(int level, iVec& index, int var);
};

#endif