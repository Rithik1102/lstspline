// NSpline.cpp: implementation of the NSpline class.
// Also implements LSNSpline, LSNSplineB and LSNSpline2
//////////////////////////////////////////////////////////////////////

#include "NSpline.h"

const double Tolerance0=1.e-5;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
int	NSpline::Init(int M, int n)
{
	N=n;
	knots.resize(N);
	auxknots=knots;
	Naux=N;
//	multiplicity.newsize(N);
	order=M;
	deficiency=1;
	previousL=0;

	return 0;
}

int NSpline::Init(int M,dVec &data)
{
		N=data.size();
		knots=data;	
		auxknots.resize(N);
		int k,i,j=0;
		deficiency=1;

		invalid=0;

		auxknots[0]=knots[0];
		k=1;
		for(i=1;i<knots.size();i++) {
			if(knots[i]>auxknots[j]) {j++; auxknots[j]=knots[i]; k=1; } 
			else {if(knots[i]<auxknots[j]) knots[i]=auxknots[j]; // for safety
				k++; if(deficiency<k) deficiency=k;}
		}
		Naux=j+1;

		order=M;
		previousL=0;
//		multiplicity.newsize(N);
		return 0;
}


void NSpline::Value(double t, dVec &B, iVec &index)
{
	Value(t,B,index,order);
}

void NSpline::Value(double t, dVec &B, iVec &index, int M)
{
	// index contains [start,end]
	// the indices of nonzero elements of Bm, Bm-1, Bm-2,...


	int i,j,k,l;

		// find l
		index.resize(2);
		FindLp(t,l,&knots); // l is 1-based
		k=0;
		// process multiplicity
		while(l>1 && Mknots(l-1) == Mknots(previousL)) {l--; k++;}
		//previousL=l;

		index[0]=l+1-M;
		index[1]=l+k;
		B.resize(M+1+k);

		fill(B.begin(), B.end(),0.0);
		B[M-1]=0;
		B[M+k-1]=1.; 
		//if(k) B(M+1)=1.;

		for(j=2;j<=M;j++) {
			for(i=M-j+1; i<=M+k;i++) {

				if(Mknots(i+l+j-1-M)-Mknots(i+l-M) < Tolerance0) 
				{
					B[i-1]=(t<= Mknots(i+l+j-M) && t>=Mknots(i+l+1-M-1) ? pow( (Mknots(i+l+j-M)-t)/(Mknots(i+l+j-M)-Mknots(i+l+1-M-1)) ,j-1):0.0);
				}
				else if(Mknots(i+l+j-M)-Mknots(i+l+1-M) < Tolerance0) 
				{
					B[i - 1] =(t< Mknots(i+l+j-1-M+1) && t>=Mknots(i+l-M) ?pow((t-Mknots(i+l-M) )/(Mknots(i+l+j-1-M+1)-Mknots(i+l-M)),j-1):0.0);}
				else
					B[i - 1] = B[i - 1] *(t-Mknots(i+l-M) )/(Mknots(i+l+j-1-M)-Mknots(i+l-M))
					+ B[i ] *(Mknots(i+l+j-M)-t)/(Mknots(i+l+j-M)-Mknots(i+l+1-M));
			}
		}
}

void NSpline::FindL(double t, int &l, dVec *Knots, int ri, int le) // l is 1-based
{
	int  u,mid;
	//l=0; u=Knots->size()-1;
	l=le; u=ri-1;
l2:
	if(u-l <=1) { l++; previousL=l; return;}
	mid=(l+u)/2;
	if(t<(*Knots)[mid]) u=mid; else l=mid;
	goto l2;
}

void NSpline::FindLp(double t, int &l, dVec *Knots) // l is 1-based
{
	// using previous l
	if(t<=(*Knots)[0]) {l=1; return;}
	if(t>=(*Knots)[GetSize()-1]) {l=GetSize(); return;}

	if(previousL==0) {FindL(t,l,Knots,GetSize()); return;}
	l=previousL;
	if(t>(*Knots)[l-1]) {
		if(t<(*Knots)[l+0]) { return;}
		if(t<(*Knots)[l+1]) {l++; previousL=l; return;}
		FindL(t,l,Knots,  GetSize() , l+1 ); return;
	}
	if(l>1 && t>(*Knots)[l-2]) {l--; previousL=l; return; }
	FindL(t,l,Knots, l-1, 0); 
}



double NSpline::Mknots(int i)
{	// this metod returns knots(i), unless i is <0 or > N
	// in this case the neighbouring knots are mirrowed
	if(i>=1 && i<=N) return knots[i-1];
	if(i<1) return(2*auxknots[0]-Mknotsaux(2-i));
	return Mknotsaux(i+Naux-N);//(2*auxknots(Naux)-Mknotsaux(2*Naux-i));
}

double NSpline::Mknotsaux(int i)
{	// this metod returns knots(i), unless i is <0 or > N
	// in this case the neighbouring knots are mirrowed
	if(i>=1 && i<=Naux) return auxknots[i-1];
	if(i<1) return(2*auxknots[1-1]-auxknots[2-i-1]);
	return (2*auxknots[Naux-1]-Mknotsaux(2*Naux-i));
}



// LSNSpline
double		LSNSpline::Value(double t)		// value of the spline
{
	if(invalid) return 0;

	NSpline::Value(t,Bsplines,Index);
	// now calculate
	double f=0;
	int j,i2;

	j=0;

	for(i2=Index[0];i2<=Index[1];i2++) 
		f += coefs[i2 - startindex]*Bsplines[j++];

	return f;
}


double		LSNSpline::ValueDer(double t,  int Nders) // value of the spline and derivatives
{
	if(invalid) return 0;
	if(Nders==0) return Value(t);

	if(deficiency > order - Nders) return 0;
	NSpline::Value(t,Bsplines,Index,order-Nders);
	// now calculate
	double f=0;
	int j,i2;

	j=0;

	for(i2=Index[0];i2<=Index[1];i2++) 
		f += coefs_der[i2 - startindex /*+ Nders*/ ][ Nders-1]*Bsplines[j++];

	for(i2=1;i2<=Nders;i2++)
		f *= (order-i2);

	return f;
}

	// assigns bounds
void		LSNSpline::Bounds(dVec& td, dVec& td1)
{
	der_boundlow=td;
	der_boundup=td1;
	der_boundlow_scalar_specified=der_boundup_scalar_specified=0;
}
void		LSNSpline::Bounds(double td, double td1)
{
	der_boundlow_scalar=td;
	der_boundup_scalar=td1;
	der_boundlow_scalar_specified=der_boundup_scalar_specified=1;
}
void		LSNSpline::Boundlow(dVec& td) 
	{der_boundlow=td; der_boundlow_scalar_specified=0;}
void		LSNSpline::Boundup(dVec& td) 
	{der_boundup=td; der_boundup_scalar_specified=0;}
void		LSNSpline::Boundlow(double td) 
	{der_boundlow_scalar=td; der_boundlow_scalar_specified=1;}
void		LSNSpline::Boundup(double td) 
	{der_boundup_scalar=td; der_boundup_scalar_specified=1;}

int	LSNSpline::Build(int kind, int n, dVec& t,  dVec& x, dVec& y, dVec& xe, dVec& ye) 
// kind of constructor
{
	NSpline::Init(n,t);
	if(deficiency > order) invalid=1;
	if(invalid) return INVALID_CODE;
	equality=1;
	return FindCoefficients(kind,x,y,xe,ye);
}

int	LSNSpline::Build(int kind, int n, dVec& t,  dVec& x, dVec& y) // kind of constructor
{
	NSpline::Init(n,t);
	if(deficiency > order) invalid=1;
	if(invalid) return INVALID_CODE;
	equality=0;
	return FindCoefficients(kind,x,y, x,y); // dummy 4,5 params
}

dVec LSNSpline::RowOfValues(double t)
{
	dVec	R(size,0);
	int i2,j1=0;
	NSpline::Value(t, Bsplines, Index);

	//R=0;
	int I=Index[1];
	if(I-startindex>=size) I--;
//	for(i2=startindex;i2<Index[0];i2++) R[i2 - startindex]=0;
	for(i2=Index[0];i2<=I;i2++) R[i2 - startindex]=Bsplines[j1++];
//	for(i2=Index[1]+1;i2 < size+startindex; i2++) R[i2 - startindex]=0;
	return R;
}

int	LSNSpline::FindCoefficients(int kind, dVec& x, dVec& y, dVec& xe, dVec& ye)
{
	// Build the matrix
	Nknots=N; // or N-1
	size=order + Nknots - 1; // B N-th not counted, as b=tN
	startindex=(2-order)+0;			 // where c starts

	LSNSpline::kind=kind;

	rowsA =		x.size();
	rowsEQ=0;	if(equality) rowsEQ=xe.size();
	rowsG=0;	
	if(kind>0 && kind<6) rowsG=size-1;//N-1;//Naux-1; // only normal knots are checked, except 1
	else if(kind>=6 && kind <11) rowsG=size-1;

	int i,j;

	if(kind==5 || kind==10) rowsG *=2;

	// 0 = LS spline, 1 = monotone LS, 2 = monotone dec
	// 3 = monotone, strict, 4= monotone dec strict
	// 5 = bounded derivative
	// 6 = concave up, 7= concave down
	// 8 = concave up strict, 9 = concave down strict
	// 10 = second derivative bounded
	// 11 = LS spline, but in T-bar representation (for testing only)
	long rows = rowsEQ+rowsA+rowsG;

	dMat equations(size+1,rows);
	// matrix w (M,N+1), M=ME+MA+MG, input to to LSEI algorithms
//(E  F) 
//(A  B) 
//(G  H) 

	// fill the matrix
	int i1,j1,i2;
	// matrix E
	for(i=0;i<rowsEQ;i++)
	{
			//j=coefs_division(1);   //startindex+N+1;
// i-th row
			j1=0;
			NSpline::Value(xe[i], Bsplines, Index);
			for(i2=startindex;i2<Index[0];i2++) equations[i2 - startindex][i]=0;
			for(i2=Index[0];i2<=Index[1];i2++) 
				equations[i2 - startindex][i]=Bsplines[j1++];
			for(i2=Index[1]+1;i2 < size+startindex; i2++) equations[i2 - startindex][i]=0;

			equations[(int)size][i]=ye[i]; // RHS
	}

	// same for matrix A
	for(i=0;i<rowsA;i++)
	{
			i1=i+rowsEQ;
			//j=coefs_division(1);  //startindex+N+1;
// i-th row
			j1=0;
			NSpline::Value(x[i], Bsplines, Index);
			for(i2=startindex;i2<Index[0];i2++) equations[i2 - startindex][i1]=0;
			for(i2=Index[0];i2<=Index[1];i2++) 
				equations[i2 - startindex][i1]=Bsplines[j1++];
			for(i2=Index[1]+1;i2 < size + startindex; i2++) equations[i2 - startindex][i1]=0;

			equations[(int)size][i1]=y[i]; // RHS
	}

	// and now matrix of inequalities
	i1=rowsEQ+rowsA;

	if(der_boundup_scalar_specified) 
		{der_boundup.resize(size,der_boundup_scalar);}
	if(der_boundlow_scalar_specified) 
		{der_boundlow.resize(size, der_boundlow_scalar);}

	int bound_up_specified =(der_boundup.size()>0);
	int bound_low_specified =(der_boundlow.size()>0);

	//??
//	if(deficiency>1 && kind==5) kind=1;
	if(deficiency>1 || kind==1 || kind==2 || kind==6 || kind==7) 
							bound_low_specified=bound_up_specified=0;

	int sign=1;
	if(kind==4 || kind==2 || kind==7 || kind==9) sign=-1; 

	switch(kind) {
	case 0: break;
	case 1: 
	case 3:
		for(i=1;i<rowsG;i++) { // first spline coef not important
			for(j=0;j<=size;j++) equations[j][i1]=0;
			equations[i][i1]=sign;
			equations[(int)size][i1]=sign*((bound_low_specified && der_boundlow.size()>=i) ? der_boundlow[i-1]:0);
			i1++;	
		}
		break;
	case 2: // from up
	case 4:		// rowsG must be size + startindex
		for(i=1;i<rowsG;i++) { // first spline coef not important
			for(j=0;j<=size;j++) equations[j][i1]=0;
			equations[i][i1]=sign;
			equations[(int)size][i1]=sign*((bound_up_specified && der_boundup.size()>=i) ? der_boundup[i-1]:0);
			i1++;	
		}
		break;
	case 5: // monotone bounded, no deficiency
		for(i=1;i<rowsG/2;i++) {
			i2=(i1+rowsG/2);
			for(j=0;j<=size;j++) equations[j][i1]=equations[j][i2]=0;
			equations[i][i1]=1;
			equations[i][i2]=-1.;

			equations[(int)size][i1]=((bound_low_specified && der_boundlow.size()>=i) ? der_boundlow[i-1]:0);
			equations[(int)size][i2]= - ((bound_up_specified && der_boundup.size()>=i) ? der_boundup[i-1]:0);
			i1++;
		}
		break;
	case 6: // low
	case 8:
		for(i=2;i<rowsG;i++) { // first spline coef not important
			for(j=0;j<=size;j++) equations[j][i1]=0;
			equations[i][i1]=sign;
			equations[(int)size][i1]=sign*((bound_low_specified && der_boundlow.size()>=i) ? der_boundlow[i-2]:0);
			i1++;	
		}
		break;
	case 7: // up
	case 9:
		for(i=2;i<rowsG;i++) { // first spline coef not important
			for(j=0;j<=size;j++) equations[j][i1]=0;
			equations[i][i1]=sign;
			equations[(int)size][i1]=sign*((bound_up_specified && der_boundup.size()>=i) ? der_boundup[i-2]:0);
			i1++;	
		}
		break;
	case 10: // convex bounded, no deficiency
		for(i=2;i<rowsG/2;i++) {
			i2=(i1+rowsG/2);
			for(j=0;j<=size;j++) equations[j][i1]=equations[j][i2]=0;
			equations[i][i1]=1;
			equations[i][i2]=-1.;

			equations[(int)size][i1]=((bound_low_specified && der_boundlow.size()>=i) ? der_boundlow[i-2]:0);
			equations[(int)size][i2]= - ((bound_up_specified && der_boundup.size()>=i) ? der_boundup[i-2]:0);
			i1++;
		}
		break;

	}

// Now transform it to T-splines by multiplying by L

	for(i=0;i<rowsEQ+rowsA;i++) {
		j1=size-1;
		for(j=j1; j> 0; j--) equations[j-1][i] += equations[j][i];
	}

	switch(kind) {
	case 3:
	case 4:
	case 5:
		// calculate matrix D and multiply the right hand side
		computeD();
		multiplyD(equations);
		break;
// convex splines
	case 6:
	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
		// Calculate matrix D L- and multiply T-splines  (T-bar representation)
		computeD();
		for(i=0;i<rowsEQ+rowsA;i++) {
			j1=size-1;
			equations[j1][i] *= D[0][size-1];
			for(j=j1; j> 0; j--) {
				equations[j-1][i] = equations[j][i] + equations[j-1][i] * D[0][j-1];
			} 
		}
		// calculate matrix D and multiply the right hand side
		computeD(order-1);
		multiplyD(equations);
		break;
	}

	double rnorme, rnorml;
	long mode;

	progopts.resize(4);
	progopts[0]=4;
	progopts[1]=2; // scale the matrix entries
	progopts[2]=1;
	progopts[3]=1;

	coefs.resize(size,0.0);	//coefs=0;

// save a copy just in case
	dMat equationsC(equations);

	m_LSEI.lsei_(equations.data(), &rows, &rowsEQ,  &rowsA, &rowsG, &size, 
		progopts.data(), coefs.data(), &rnorme, &rnorml, &mode);

	if(mode != 0) { // in case LSEI did not find the feasible solution
		// try a more expensive method

//		printf("Attempting wnnls\n");
		//equations.resize(0,0);

		j1=size+rowsG;  // new size of the syetem
		i1=rowsEQ+rowsG; //new size of equality conditions
		long j1S=j1; // as above
		long i1S=i1;
		int j2;
		dMat equations1(j1+1,rows);
		coefs.resize(j1); // will contain slack variables

		//progopts[0]=1; // no options
		progopts[1]=6; // scale

		for(j=0;j<size+1;j++) {
			for(i=0; i<rowsEQ; i++) equations1[j][i]=equationsC[j][i];
			for(i=0,i2=rowsEQ,j2=rowsEQ+rowsA; i<rowsG; i++,i2++,j2++) 
				equations1[j][i2]=equationsC[j][j2];
			for(i=0,i2=rowsEQ+rowsG,j2=rowsEQ; i<rowsA; i++,i2++,j2++) 
				equations1[j][i2]=equationsC[j][j2];
		}
		for(i=0;i<rows;i++) equations1[j1][i]=equations1[(int)size][i]; // move RHS
		for(j=size;j<j1;j++) for(i=0;i<rows;i++) equations1[j][i]=0;  // clear
			// slack variables
		for(j=size,i=rowsEQ;j<j1;j++,i++) equations1[j][i]=-1; 

	
		m_LSEI.wnnls_(equations1.data(), &rows, &i1S,  &rowsA, &j1S, &size, // coefs startring with size+1 are >0
			progopts.data(), coefs.data(),&rnorml, &mode);

		// coefs is now larger: contains slack vars. So what?
		// clean up
		equations1.newsize(0,0);
		equationsC.newsize(0,0);
	}

	// interpret the results
	if (kind>=6) {
// translate coefs to T-spline repsn. using b=DLc
		for(i=1;i < size;i++) coefs[i] += coefs[i-1];
		for(i=1;i <= size;i++) coefs[i-1] *=  D[0][i-1] ;
	}

// translate coefs to B-spline repsn. using a=Lb
	for(i=1;i < size;i++)  coefs[i] += coefs[i-1];
	
// now find coefs of the derivatives
	coefs_der.newsize(size,order);
	for(i=1;i<size;i++)
		coefs_der[i][0]=(coefs[i]-coefs[i-1])/(Mknots(i+order-1+startindex)-Mknots(i+startindex));
	coefs_der[0][0]=0;
	for(j=1;j<order;j++) {
		for(i=j+1;i<size;i++)
		 coefs_der[i][j]=(coefs_der[i][j-1]-coefs_der[i-1][j-1])/
			(Mknots(i+order-j+startindex-1)-Mknots(i+startindex));
	}

	equations.newsize(0,0);
	return mode;
}

void LSNSpline::computeD() { computeD(order); }

void LSNSpline::computeD(int locorder)
{
	int i,k;
	k=order-locorder;
	D[k].resize(size);
	// this is only done for defficiency less than order-1
	if(deficiency < locorder)
		for(i=1;i<= size; i++) D[k][i-1]=(Mknots(i)-Mknots(i-locorder+1))/(locorder-1);
/*	else
	  for(k=1;k<=deficiency;k++)
		  for(i=1;i< coefs_division(k);i++) {
			 D(j)=(Mauxknots(i,k)-Mauxknots(i-(order-k),k))/(order-k);
			 j++;
		  }
		  */
}

void LSNSpline::multiplyD(dMat& equations)
{
	int i;
	int sh=rowsEQ+rowsA;
	// multiplication of RHS by Dk, only for deficiency 1
//	if(deficiency < order) return;

	switch(kind){
	case 5:
		for(i=2;i<=N;i++) 
			equations(size+1,i+sh+rowsG/2-1) *= D[0][i-1];
	case 3:// this is done for one-side bounds case as well
	case 4:
		for(i=2;i<=N;i++)
			  equations(size+1,i+sh-1) *= D[0][i-1];
		break;
	case 10:
		for(i=2;i<=N;i++) 
			equations(size+1,i+sh+rowsG/2-1) *= D[1][i-1];
	case 9:
	case 8:
		for(i=2;i<=N;i++)
			  equations(size+1,i+sh-1) *= D[1][i-1];
		break;
	}

}



// LSNSpline
double		LSNSpline2::Value(double t1,double t2)		// value of the spline
{
	if(invalid) return 0;

	spl1.NSpline::Value(t1,spl1.Bsplines,spl1.Index);
	spl2.NSpline::Value(t2,spl2.Bsplines,spl2.Index);

	// now calculate
	double f=0;
	double g;
	int j,i2, j1, i21 ;

	j=j1=0;

	for(i21=spl2.Index[0];i21<=spl2.Index[1];i21++) // second
	{
		g=0;	j=0;
		for(i2=spl1.Index[0];i2<=spl1.Index[1];i2++) // first
			g += coefs(i2 - spl1.startindex+1, i21-spl2.startindex+1) * (spl1.Bsplines[j++]);
		f += g*spl2.Bsplines[j1++];
	}
	return f;
}


double		LSNSpline2::ValueDer(double t1,  double t2, int var) 
// value of the spline and 1st derivative with respect to var (0 or 1)
{
	if(invalid) return 0;

	switch(var) {
	case 0:
		if(spl1.deficiency > spl1.order-1) return 0;
		break;
	case 1:
		if(spl2.deficiency > spl2.order-1) return 0;
		break;
	}

	double g,f=0;
	int j,i2, j1,i21;


	switch (var) {
	case 0:
		spl1.NSpline::Value(t1,spl1.Bsplines,spl1.Index,spl1.order-1);
		spl2.NSpline::Value(t2,spl2.Bsplines,spl2.Index);
		break;
	case 1:
		spl1.NSpline::Value(t1,spl1.Bsplines,spl1.Index,spl1.order);
		spl2.NSpline::Value(t2,spl2.Bsplines,spl2.Index,spl2.order-1);
		break;
	}

	// now calculate
	j=j1=0;
	for(i21=spl2.Index[0];i21<=spl2.Index[1];i21++) // second
	{
		g=0;	j=0;
		for(i2=spl1.Index[0];i2<=spl1.Index[1];i2++) // first
			g += (coefs_der[var])(i2 - spl1.startindex+1, i21-spl2.startindex+1)
				* spl1.Bsplines[j++];
		f += g*spl2.Bsplines[j1++];
	}

	switch (var) {
	case 0:
		f *= (spl1.order-1); break;
	case 1:
		f *= (spl2.order-1); break;
	}

	return f;
}

int LSNSpline2::Build(int kind, int n1, int n2, dVec& t1, dVec& t2, dMat& x,  dVec& y) 
{
	dMat xe(0,0); // dummy interpolation conditions
	dVec ye(0);
	return Build(kind,n1,n2,t1,t2,x,y,xe,ye);
}

int LSNSpline2::Build(int kind, int n1, int n2, dVec& t1,  dVec& t2, 
					   dMat& x,  dVec& y, dMat& xe, dVec& ye) 
// kind of constructor
{
	// we only process bi-spline of order 2 (bilinear)
	if(kind>0) n1=n2=2;
// It is possible to handle bicuadratic splines
// but the construction of G needs to be modified
// instead of (0|I) x L we have to use I x NL
// where N is the matrix of values of cuadratic spline at the knots
// i.e. the values of T-spline (and similar for the second variable)
// Almost L but the upper codiagonal is not 0
// this is postponed for the future.

	invalid=1;
	spl1.NSpline::Init(n1,t1);
	if(spl1.invalid) return INVALID_CODE;
	spl2.NSpline::Init(n2,t2);
	if(spl2.invalid) return INVALID_CODE;

	if(spl1.deficiency > spl1.order || spl2.deficiency > spl2.order) return INVALID_CODE;

	invalid=0;
	return FindCoefficients(kind,x,y,xe,ye);
}

void MultiplyByL(dVec& r)
{
	// multiplies a row vector r by low triangular L and returns r
	int j;
	for(j=r.size()-1; j> 0; j--) 	r[j-1] += r[j];
}

dVec TensorProductRow(dVec& r1, dVec& r2)
{
	// returns a row which is the tensor product of 2 rows
	dVec r(r1.size()*r2.size());
	int i,j,k;
	k=0;
	for(i=0;i<r1.size(); i++) {
		for(j=0;j<r2.size(); j++) {
			r[k]=r1[i]*r2[j];
			k++;
		}
	}
	return r;
}

void TensorProduct(dMat& m1, dMat& m2, dMat& m)
{
	// returns a tensor product of two matrices in m
	int i1,i2,j1,j2,i,j,I,J;
	I=m2.dim(1); J=m2.dim(2);
	for(i1=0;i1<m1.dim(1); i1++) 
	for(j1=0;j1<m1.dim(2); j1++) 

		for(i2=0;i2<I; i2++) 
		for(j2=0;j2<J; j2++) {
			i=i1*I+i2; j=j1*J+j2;	
			m[i][j] = m1[i1][j1]*m2[i2][j2];
		}
}

int LSNSpline2::FindCoefficients(int kind, dMat& x, dVec& y, dMat& xe, dVec& ye)
{
	// Build the matrix
	spl1.Nknots=spl1.NSpline::N-1; // or N-1
	spl2.Nknots=spl2.NSpline::N-1; // or N-1

	spl1.startindex=(2-spl1.order)+0;			 // where c starts
	spl2.startindex=(2-spl2.order)+0;			 // where c starts

	spl1.size=spl1.order + spl1.Nknots - 1; // B N-th not counted, as b=tN
	spl2.size=spl2.order + spl2.Nknots - 1; // B N-th not counted, as b=tN

	size = spl1.size * spl2.size;

	LSNSpline2::kind=kind;

	// 0 = LS spline, 
	// 1*   - will be first var + 10* will be the second
	// 1 = monotone LS, 2 = monotone dec
	// 3 = monotone, strict, 4= monotone dec strict
	// 5 = bounded derivative
	// no convexity - too much
	// so that
	// 11  means both monotone, 21 means y dec and x inc, etc

	div_t d=div(kind,10);

// handle erroneous input (like 54, 25, etc)
	if(d.rem==5 && d.quot!=5 ) {
		switch(d.quot) {
		case 0: // unconstrained
			der_bounddown2=-Infinity;
			der_boundup2=Infinity;
			break;
		case 1: // increasing
			der_bounddown2=0;
		case 3:
			der_boundup2=Infinity;
			break;
		case 2:
			der_boundup2=0;
		case 4:
			der_bounddown2=-Infinity;
			break;
		}
		kind=55; d=div(kind,10); // force kind 55
	}
	if(d.rem!=5 && d.quot==5 ) {
		switch(d.rem) {
		case 0: // unconstrained
			der_bounddown1=-Infinity;
			der_boundup1=Infinity;
			break;
		case 1: // increasing
			der_bounddown1=0;
		case 3:
			der_boundup1=Infinity;
			break;
		case 2:
			der_boundup1=0;
		case 4:
			der_bounddown1=-Infinity;
			break;
		}
		kind=55; d=div(kind,10); // force kind 55
	}


	rowsA =		x.dim(1);
	rowsEQ=0;	rowsEQ=xe.dim(1); // must be set
	rowsG=0;	
	if(kind>0) rowsG=(spl1.size-1)*(spl2.size)+(spl1.size)*(spl2.size-1);//N-1;//Naux-1; // only normal knots are checked, except 1
	if(kind==55 ) rowsG*=2;

	int i,j,k;
	int i1,j1,i2,j2;

	long rows = rowsEQ+rowsA+rowsG;

	// matrix w (M,N+1), M=ME+MA+MG, input to to LSEI algorithms
	dMat equations(size+1,rows);
//(E  F) 
//(A  B) 
//(G  H) 
FILE *ff;
//		ff=fopen("dump.txt","w");

	// fill the matrix
	dVec R1(spl1.size),R2(spl2.size), R3;
	// matrix E
	for(i=0;i<rowsEQ;i++)
	{
// i-th row
			R1=spl1.RowOfValues(xe[i][0]); 
			MultiplyByL(R1);
			R2=spl2.RowOfValues(xe[i][1]); MultiplyByL(R2);
			R3=TensorProductRow(R1,R2);
			for(j=0;j<size;j++) equations[j][i]=R3[j];
			equations[(int)size][i]=ye[i];
	}
//fclose(ff);
	// same for matrix A
	for(i=0;i<rowsA;i++)
	{
			i1=i+rowsEQ;
			R1=spl1.RowOfValues(x[i][0]); MultiplyByL(R1);
			R2=spl2.RowOfValues(x[i][1]); MultiplyByL(R2);
			R3=TensorProductRow(R1,R2);
			for(j=0;j<size;j++) equations[j][i1]=R3[j];
			equations[(int)size][i1]=y[i];
	}

	// and now matrix of inequalities
	i1=rowsEQ+rowsA;

	int bound_up_specified1 =(der_boundup1 !=0) ;
	int bound_low_specified1 =(der_bounddown1 !=0) ;
	int bound_up_specified2 =(der_boundup2 !=0) ;
	int bound_low_specified2 =(der_bounddown2 !=0) ;

	int sign1=1;
	int sign2=1;
	if(d.rem==4 || d.rem==2)  sign1=-1;
	if(d.quot==4 || d.quot==2) sign2=-1; 

	// current bounds
	double bound1=der_bounddown1;
	if(d.rem==4) bound1=der_boundup1; 

	double bound2=der_bounddown2;
	if(d.quot==4) bound2=der_boundup2; 

	// handles unconstrained case in one variable. Not very nice though: redundant eqs.
	if(d.rem==0)  bound1=-Infinity;
	if(d.quot==0) bound2=-Infinity; 

	Tensor T;
	dMat I1(spl1.size-1,spl1.size),I2(spl2.size,spl2.size);

	i1++; // using fortran indices
	switch(kind) {
	case 0: break; // do nothing here

	case 55: // special twice as many
		// first process lower bounds
		sign1=sign2=1;
		bound1=der_bounddown1;
		bound2=der_bounddown2;

		// derivative wrt to x
		//I1=0; I2=0;
		fill(I1.begin(), I1.end(), 0);
		fill(I2.begin(), I2.end(), 0);

		for(i=0;i<spl1.size-1;i++) I1[i][i+1]=1; // matrices I and L
		for(i=0;i<spl2.size;i++) for(j=0;j<=i;j++)  I2[i][j]=1;
		// build tensor
		T.AddTail(I1); T.AddTail(I2);
		// populate matrix of eqns
		for(i=1;i<=(spl1.size-1)*(spl2.size);i++) { 
			for(j=1;j<=(spl1.size)*(spl2.size);j++) 
				equations(j,i1)=sign1*T(i,j);
			equations((int)(size+1),i1)=sign1*bound1; // RHS
			i1++;	// next row
		}
// cleanup
		T.RemoveTail(); T.RemoveTail(); 

		I1.newsize(spl1.size,spl1.size); I2.newsize(spl2.size-1,spl2.size);
		// derivative wrt to y
		//I1=0; I2=0;
		fill(I1.begin(), I1.end(), 0);
		fill(I2.begin(), I2.end(), 0);

		for(i=0;i<spl1.size;i++) for(j=0;j<=i;j++) 	I1[i][j]=1;
		for(i=0;i<spl2.size-1;i++)  I2[i][i+1]=1;
		// build tensor
		T.AddTail(I1); T.AddTail(I2);
		// populate matrix of eqns
		for(i=1;i<=(spl1.size)*(spl2.size-1);i++) { 
			for(j=1;j<=(spl1.size)*(spl2.size);j++) 
				equations(j,i1)=sign2*T(i,j);
			equations((int)(size+1),i1)=sign2*bound2; // RHS
			i1++;	
		}
		T.RemoveTail(); T.RemoveTail(); 
		// continue with the upper bounds
		bound1=der_boundup1; 
		bound2=der_boundup2; 
		sign1=sign2= -1;
		I1.newsize(spl1.size-1,spl1.size); I2.newsize(spl2.size,spl2.size);
		//break;  no break here

	default:
		// rowsG must be size + startindex unless continuing kind 55
		// derivative wrt to x
		//I1=0; I2=0;
		fill(I1.begin(), I1.end(), 0);
		fill(I2.begin(), I2.end(), 0);
		for(i=0;i<spl1.size-1;i++) I1[i][i+1]=1;
		for(i=0;i<spl2.size;i++) for(j=0;j<=i;j++)  I2[i][j]=1;
		// build tensor
		T.AddTail(I1); T.AddTail(I2);
		// populate matrix of eqns
		for(i=1;i<=(spl1.size-1)*(spl2.size);i++) { 
			for(j=1;j<=(spl1.size)*(spl2.size);j++) 
				equations(j,i1)=sign1*T(i,j);
			equations((int)(size+1),i1)=sign1 * bound1;  // RHS
			i1++;	
		}
// cleanup
		T.RemoveTail(); T.RemoveTail(); 
		I1.newsize(spl1.size,spl1.size); I2.newsize(spl2.size-1,spl2.size);

		// derivative wrt to y
		//I1=0; I2=0;
		fill(I1.begin(), I1.end(), 0);
		fill(I2.begin(), I2.end(), 0);
		for(i=0;i<spl1.size;i++) for(j=0;j<=i;j++) 	I1[i][j]=1;
		for(i=0;i<spl2.size-1;i++)  I2[i][i+1]=1;
		// build tensor
		T.AddTail(I1); T.AddTail(I2);
		// populate matrix of eqns
		for(i=1;i<=(spl1.size)*(spl2.size-1);i++) { 
			for(j=1;j<=(spl1.size)*(spl2.size);j++) 
				equations(j,i1)=sign2*T(i,j);
			equations((int)(size+1),i1)=sign2* bound2; // RHS
			i1++;	
		}
		// cleanup
		T.RemoveTail(); T.RemoveTail(); 
		I1.newsize(0,0); I2.newsize(0,0);
		break;
	}

	// process derivative wrt to variable x 
	// (only if the derivative value is given)
	switch(d.rem) {
	case 3:
	case 4:
	case 5:
		spl1.computeD();	// multiply by D
		i1=rowsEQ+rowsA; j=2;
		for(i=2;i<=(spl1.size-1)*(spl2.size)+1;i++)
		{	if(i>(spl2.size+1)*(j-1)) j++;
			equations((int)(size+1),i+i1-1) *= spl1.D[0][j-1];  
		}
	}

	// process der wrt to second variable y
	switch(d.quot) {
	case 3:
	case 4:
	case 5:
		spl2.computeD();
		i1=rowsEQ+rowsA+(spl1.size-1)*(spl2.size);	j=1;
		for(i=2;i<=(spl1.size)*(spl2.size-1)+1;i++)
		{	j++; if(j>spl2.size) j=2;
			equations((int)(size+1),i+i1-1) *= spl2.D[0][j-1];
		}
	}

	if(kind == 55) // repeat that (bounds from both sides)
	{
		i1=rowsEQ+rowsA+(rowsG /2);  j=2;
		for(i=2;i<=(spl1.size-1)*(spl2.size)+1;i++)
		{	if(i>(spl2.size+1)*(j-1)) j++;
			equations((int)(size+1),i+i1-1) *= spl1.D[0][j-1];
		}
		i1=rowsEQ+rowsA +(rowsG /2) + (spl1.size-1)*(spl2.size); j=1;
		for(i=2;i<=(spl1.size)*(spl2.size-1)+1;i++)
		{	j++; if(j>spl2.size) j=2;
			equations((int)(size+1),i+i1-1) *= spl2.D[0][j-1];
		}
	}
	
	
	double rnorme, rnorml;
	long mode;

	spl1.progopts.resize(4);
	spl1.progopts[0]=4;
	spl1.progopts[1]=2; // scale all entries
	spl1.progopts[2]=1;
	spl1.progopts[3]=1;

	spl1.coefs.resize(size);	fill(spl1.coefs.begin() , spl1.coefs.end(), 0);

	// make a copy of equation just in case for wmmls (it is destroyed in lsei)
	dMat equationsC(equations);

	// attempt to solve the system quickly
	spl1.m_LSEI.lsei_(equations.data(), &rows, &rowsEQ,  &rowsA, &rowsG, &size, 
		spl1.progopts.data(), spl1.coefs.data(), &rnorme, &rnorml, &mode);
//mode=2;
	if(mode == 0) { // in case LSEI did not find the feasible solution
		// try a more expensive method

//		printf("Attempting wnnls\n");
		equations.newsize(0,0); // cleanup

		j1=size+rowsG;  // new size ( adding slack variables)
		i1=rowsEQ+rowsG; // new size of equalities system
		long j1S=j1;  // same: for passing parameters
		long i1S=i1;

		dMat equations1(j1+1,rows); // new system of eqs
		spl2.coefs.resize(j1);     // here goes the solution
		spl1.progopts[0]=1; // 1 no options, 6 - scale
		spl1.progopts[1]=6; // 1 no options, 6 - scale

		// copy and restructure matrix columnwise (faster, will use system cache)
		for(j=0;j<size+1;j++) {
			for(i=0; i<rowsEQ; i++) equations1[j][i]=equationsC[j][i];
			for(i=0,i2=rowsEQ,j2=rowsEQ+rowsA; i<rowsG; i++,i2++,j2++) 
				equations1[j][i2]=equationsC[j][j2];
			for(i=0,i2=rowsEQ+rowsG,j2=rowsEQ; i<rowsA; i++,i2++,j2++) 
				equations1[j][i2]=equationsC[j][j2];
		}
		for(i=0;i<rows;i++) equations1[j1][i]=equations1[(int)size][i]; // move RHS
		for(j=size;j<j1;j++) for(i=0;i<rows;i++) equations1[j][i]=0; // clear other entries
			// slack variables on diagonal of G
		for(j=size,i=rowsEQ; j<j1; j++,i++) equations1[j][i]=-1; 
/*
	ff=fopen("dump.txt","a");
	for(i=0;i<rows;i++) {
		for(j=0;j<j1+1;j++)
			fprintf(ff,"%3.1f ",equations1[j][i]);
		fprintf(ff,"\n");
	}
	fclose(ff);
*/		
		spl1.m_LSEI.wnnls_(equations1.data(), &rows, &i1S,  &rowsA, &j1S, &size, // coefs startring with size+1 are >0
			spl1.progopts.data(), spl2.coefs.data(),&rnorml, &mode);

		// just a copy
		for(j=0;j<size;j++) spl1.coefs[j]=spl2.coefs[j];

		// clean up
		equations1.newsize(0,0);
		equationsC.newsize(0,0);
	}
	equations.newsize(0,0); // cleanup

//	printf("mode %d\n",mode);

// translate coefs to B-spline repsn. using a=Lb
// this is slow explicit method. We can do better than this 
/*	dMat L1(spl1.size,spl1.size);
	dMat L2(spl2.size,spl2.size);
	dMat L(size,size); L1=0; L2=0;
	for(i=0;i<spl1.size;i++) for(j=0;j<=i;j++) L1[i][j]=1;
	for(i=0;i<spl2.size;i++) for(j=0;j<=i;j++) L2[i][j]=1;
	TensorProduct(L1,L2,L);

    spl1.coefs=L * spl1.coefs;
	L1.newsize(0,0); L2.newsize(0,0); L.newsize(0,0);
*/
// translate to matrix of coefs
	coefs.newsize(spl1.size,spl2.size); k=0;
	for(i=0;i<spl1.size;i++)
		for(j=0;j<spl2.size;j++) coefs[i][j]=spl1.coefs[k++];

//	FILE *ff;
/*	ff=fopen("dump.txt","a");
	for(i=0;i<spl1.size;i++) {
		fprintf(ff,"\n");
	for(j=0;j<spl2.size;j++)
		fprintf(ff,"%3.1f ",coefs[i][j]);
	}
	fclose(ff);
*/

// translate coefs to B-spline repsn. using a=Lb
// a faster method compared with direct tensor product
	spl2.coefs.resize(spl1.size); // just a temporary column
	for(i=0;i < spl2.size;i++)  {
		spl2.coefs[0]=coefs[0][i];
		for(j=1; j<spl1.size;j++) spl2.coefs[j]= coefs[j][i] + spl2.coefs[j-1];
		for(j=0;j<spl1.size;j++) coefs[j][i]=(i>0 ? coefs[j][i-1]:0)+spl2.coefs[j];	
	}

/*	ff=fopen("dump.txt","a");
	for(i=0;i<spl1.size;i++) {
		fprintf(ff,"\n");
	for(j=0;j<spl2.size;j++)
		fprintf(ff,"%3.1f ",coefs[i][j]);
	}
	fclose(ff);

*/	
// now find coefs of the derivatives
	coefs_der[0].newsize(spl1.size,spl2.size);
	coefs_der[1].newsize(spl1.size,spl2.size);
	for(i=1;i<spl1.size;i++)
	   for(j=0;j<spl2.size;j++) {
		coefs_der[0](1,j+1)=0;
		coefs_der[0](i+1,j+1)=(coefs[i][j]-coefs[i-1][j])/(spl1.Mknots(i+spl1.order-1+spl1.startindex)-spl1.Mknots(i+spl1.startindex));
	   }
	for(j=1;j<spl2.size;j++)
	   for(i=0;i<spl1.size;i++) {
		coefs_der[1](i+1,1)=0;
		coefs_der[1](i+1,j+1)=(coefs[i][j]-coefs[i][j-1])/(spl2.Mknots(j+spl2.order-1+spl2.startindex)-spl2.Mknots(j+spl2.startindex));
	   }

	   // finish here
	   // cleanup has been done earlier

	   return mode;
}


// B-spline implementation --------------------------------

double LSNSplineB::Mknots(int i)
{	// this metod returns knots(i), unless i is <0 or > N
	// in this case the neighbouring knots are mirrowed
	if(i>=1 && i<=N) return knots[i-1];
	if(i<1) return(auxknots[1-1]-DELTA_T*(Mknotsaux(2-i)-auxknots[1-1]));
	return Mknotsaux(i+Naux-N);//(2*auxknots(Naux)-Mknotsaux(2*Naux-i));
}

double LSNSplineB::Mknotsaux(int i)
{	// this metod returns knots(i), unless i is <0 or > N
	// in this case the neighbouring knots are mirrowed
	if(i>=1 && i<=Naux) return auxknots[i-1];
	if(i<1) return(auxknots[0]-DELTA_T*(auxknots[2-i-1]-auxknots[0]));
	return (auxknots[Naux-1]-DELTA_T*(Mknotsaux(2*Naux-i)-auxknots[Naux-1]));
}

int	LSNSplineB::FindCoefficients(int kind, dVec& x, dVec& y, dVec& xe, dVec& ye)
{
	// Build the matrix
	Nknots=N-1; // or N-1
	size=order + Nknots - 1; // B N-th not counted, as b=tN
	startindex=(2-order)+0;			 // where c starts

	LSNSpline::kind=kind;
	// ignores kind

	rowsA =		x.size();
	rowsEQ=0;	if(equality) rowsEQ=xe.size();
	rowsG=0;	

	int i,j;

	long rows = rowsEQ+rowsA+rowsG;
	dMat equations(size+1,rows);

	int i1,j1,i2;
	// matrix E
	for(i=0;i<rowsEQ;i++)
	{
// i-th row
			j1=0;
			NSpline::Value(xe[i], Bsplines, Index);
			for(i2=startindex;i2<Index[0];i2++) equations[i2 - startindex][i]=0;
			for(i2=Index[0];i2<=Index[1];i2++) 
				equations[i2 - startindex][i]=Bsplines[j1++];
			for(i2=Index[1]+1;i2 < size+startindex; i2++) equations[i2 - startindex][i]=0;

			equations[(int)size][i]=ye[i];
	}

	// same for matrix A
	for(i=0;i<rowsA;i++)
	{
			i1=i+rowsEQ;
			//j=coefs_division(1);  //startindex+N+1;
// i-th row
			j1=0;
			NSpline::Value(x[i], Bsplines, Index);
			for(i2=startindex;i2<Index[0];i2++) equations[i2 - startindex][i1]=0;
			for(i2=Index[0];i2<=Index[1];i2++) 
				equations[i2 - startindex][i1]=Bsplines[j1++];
			for(i2=Index[1]+1;i2 < size + startindex; i2++) equations[i2 - startindex][i1]=0;

			equations[(int)size][i1]=y[i];
	}

	// no matrix G - no constraints

	double rnorme, rnorml;
	long mode;

	progopts.resize(4);
	progopts[0]=1;
	progopts[1]=2;
	progopts[2]=1;
	progopts[3]=1;


	coefs.resize(size,0); //coefs=0;
	m_LSEI.lsei_(equations.data(), &rows, &rowsEQ,  &rowsA, &rowsG, &size, 
		progopts.data(), coefs.data(),&rnorme, &rnorml, &mode);

//	m_LSEI.wnnls_(equations[0], &rows, &rowsEQ,  &rowsA, &rowsG, &size, progopts.begin(), coefs.begin(),
//	   &rnorml, &mode);
	 	
// now find coefs of the derivatives
	coefs_der.newsize(size,order);
	for(i=1;i<size;i++)
		coefs_der[i][0]=(coefs[i]-coefs[i-1])/(Mknots(i+order-1+startindex)-Mknots(i+startindex));
	coefs_der[0][0]=0;
	for(j=1;j<order;j++) {
		for(i=j+1;i<size;i++)
		 coefs_der[i][j]=(coefs_der[i][j-1]-coefs_der[i-1][j-1])/
			(Mknots(i+order-j+startindex-1)-Mknots(i+startindex));
	}

	// cleanup
	equations.newsize(0,0);
	return mode;
}
