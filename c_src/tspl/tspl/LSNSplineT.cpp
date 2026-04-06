// LSNSplineT.cpp: implementation of the LSNSplineT class.
//
//////////////////////////////////////////////////////////////////////

#include "LSNSplineT.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void deleteVecArray(VecArray& _t)
{
	for (auto i = 0; i < _t.size(); i++)
		//delete (span<double>*)(_t[i]);
		delete (_t[i]);
	_t.clear();
}

LSNSplineT::LSNSplineT()
{

}

double	LSNSplineT::NestedSum(int level, iVec& index)		// value of the spline
{
	double g=0;
	int k;
	int j=0;
	if(level==-1) {
		// only coefficients
		return coefs(index);
	}

	LSNSpline* splt=spl[level];
	for(k=splt->Index[0]; k<= splt->Index[1]; k++) {
		index[level]=k-splt->startindex+1;      // 1-based indices, for similarity with LSNSpline2 method
		g+=NestedSum(level-1, index) * splt->Bsplines[j++];
	}
	return g;
}
double	LSNSplineT::NestedSumDer(int level, iVec& index, int var)		// value of the spline
{
	double g=0;
	int k;
	int j=0;
	if(level==-1) {
		// only coefficients
		return (*coefs_der[var])(index);
	}

	LSNSpline* splt=spl[level];
	for(k=splt->Index[0]; k<= splt->Index[1]; k++) {
		index[level]=k-splt->startindex+1;
		g+=NestedSumDer(level-1, index,var) * splt->Bsplines[j++];
	}
	return g;
}

// LSNSpline
double		LSNSplineT::Value(dVec t)		// value of the spline
{
	if(invalid) return 0;

	int k;
	LSNSpline* splt;
	for(k=0; k<dim; k++) {
		splt=spl[k];
		splt->NSpline::Value(t[k],splt->Bsplines,splt->Index);
	}

	iVec index(dim);
	double f=0;
	// now calculate
	f=NestedSum(dim-1, index);

	return f;
}


double		LSNSplineT::ValueDer(dVec t, int var) 
// value of the spline and 1st derivative with respect to var (0 or 1)
{
	if(invalid) return 0;

	if(spl[var]->deficiency > spl[var]->order-1) return 0;

	int k;
	LSNSpline* splt;
	for(k=0; k<dim; k++) {
		splt=spl[k];
		if(k==var)
		   splt->NSpline::Value(t[k],splt->Bsplines,splt->Index,splt->order-1);
		else
		   splt->NSpline::Value(t[k],splt->Bsplines,splt->Index);
	}

	// now calculate
	iVec index(dim);
	double f=NestedSumDer(dim-1, index, var);

	splt=spl[var];
	f *= (splt->order-1);

	return f;
}

//GB for now working on that part
int LSNSplineT::Build(int kind, int dim, dMat& t, dMat& data, dMat& dataE)
{
	VecArray _t, _data, _dataE;
	t.MakeVecArray(_t);
	data.MakeVecArray(_data);
	dataE.MakeVecArray(_dataE);
	int r= Build(kind, dim, _t, _data, _dataE);
// delete VecArrays
	deleteVecArray(_t);
	deleteVecArray(_data);
	deleteVecArray(_dataE);

	return r;
}


int LSNSplineT::Build(int kind, int dim, VecArray& t, VecArray& data, 
					  VecArray& dataE) 
// kind of constructor
{
	// we only process T-spline of order 2 (T-linear)

	invalid=1;
	spl.reserve(dim);
	int k;

	dVec* ttemp; 
	LSNSpline* splt;
	for(k=0;k<dim;k++) {
		splt=new LSNSpline;
		if(splt==NULL) return INVALID_CODE;
		ttemp=t[k];
		splt->NSpline::Init(2, *ttemp);
		if(splt->invalid) return INVALID_CODE; // perhaps cleanup here
		if(splt->deficiency > splt->order) return INVALID_CODE;
		spl.push_back(splt);
		//spl[k]=splt;
	}

	invalid=0;
	return FindCoefficients(kind,dim,t,data,dataE);
}

// GB easier to code as dMat, rather than vector of pointers, or using span

int LSNSplineT::FindCoefficients(int kind, int _dim, VecArray& t, VecArray& data, 
					  VecArray& dataE)
{
	// Build the matrix
	LSNSpline* splt;
	size=1;
	int k,k1;
	dim=_dim;
		FILE *ff;


	iVec sizes(dim);

	for(k=0;k<dim;k++) {
		splt=spl[k];
		splt->Nknots=splt->NSpline::N-1;           // to get correct size
		splt->startindex=(2-splt->order)+0;	
		splt->size=splt->order + splt->Nknots - 1; // B N-th not counted, as b=tN
		size *= splt->size;
		sizes[k]=splt->size;
	}

	LSNSplineT::kind=kind;

	// 0 = LS spline, 
	// 1*   - will be first var + 10* will be the second
	// 1 = monotone LS, 2 = monotone dec
	
	iVec  d(dim); // contains the kinds of splines wrt to each var
	div_t  dt;
	for(k=0;k<dim;k++) {
		dt=div(kind,10);
		d[k]=dt.rem;
		kind=dt.quot;
	}	
	kind=LSNSplineT::kind;

	rowsA =		data.size();
	rowsEQ=		dataE.size(); // must be set

	rowsG=0;	
	if(kind>0) {
		for(k=0;k<dim;k++) {
			splt=spl[k];
			rowsG += (splt->size-1) * size / splt->size;
		}
	}

	int i,j;
	int i1,j1,i2,j2;

	long rows = rowsEQ+rowsA+rowsG;

	// matrix w (M,N+1), M=ME+MA+MG, input to to LSEI algorithms
	dMat equations(size+1,rows);
//(E  F) 
//(A  B) 
//(G  H) 

	// fill the matrix
	dVec R1, R2, R3;
	// matrix E
	for(i=0;i<rowsEQ;i++)
	{
// i-th row
			k=0;
			splt=spl[k];
			R1.resize(splt->size);
			R1=splt->RowOfValues( (*(dataE[i]))[k] ); MultiplyByL(R1);
			for(k=1;k<dim;k++) {
				splt=spl[k];
				R2.resize(splt->size);
				R2=splt->RowOfValues( (*(dataE[i]))[k] ); MultiplyByL(R2);
				R1=TensorProductRow(R1,R2);
			}
			for(j=0;j<size;j++) equations[j][i]=R1[j];
			equations[(int)size][i]= (*(dataE[i]))[dim];
	}

	// same for matrix A
	for(i=0;i<rowsA;i++)
	{
			i1=i+rowsEQ;
			k=0;
			splt=spl[k];
			R1.resize(splt->size);
			R1=splt->RowOfValues( (*(data[i]))[k] ); MultiplyByL(R1);
			for(k=1;k<dim;k++) {
				splt=spl[k];
				R2.resize(splt->size);
				R2=splt->RowOfValues( (*(data[i]))[k] ); MultiplyByL(R2);
				R1=TensorProductRow(R1,R2);
			}
			for(j=0;j<size;j++) equations[j][i1]=R1[j];
			equations[(int)size][i1]= (*(data[i]))[dim];
	}
//??????????
	///////
	// and now matrix of inequalities
	i1=rowsEQ+rowsA;

	int sign1=1;
	int sign2=1;

	Tensor T;
	dMat *I1;
//	dMat I1(spl1.size-1,spl1.size),I2(spl2.size,spl2.size);

	i1++; // using fortran indices
	switch(kind) {
	case 0: break; // do nothing here
	default:
	  for(k=0;k<dim;k++) { // all derivatives
		// derivative wrt to x[k]
		for(k1=0;k1<dim;k1++) {
			splt=spl[k1];
			if(k==k1) {
				I1=new dMat(splt->size-1,splt->size,0);
				//(*I1)=0;
				for(i=1;i<splt->size;i++) (*I1)[i-1][i]=1;  // no first row
			}
			else {
				I1=new dMat(splt->size,splt->size,0);
				//(*I1)=0;
				for(i=0;i<splt->size;i++) for(j=0;j<=i;j++)  (*I1)[i][j]=1;
			}
		// build tensor		
			T.AddTail(*I1);
		} // tensor created
			switch(d[k]) {
			case 0: sign1=0; break;
			case 1: sign1=1; break;
			case 2: sign1=-1; break;
			}
		// populate matrix of eqns
		for(i=1;i<=T.dim(1);i++) { 
			for(j=1;j<=size;j++) 
				equations(j,i1) = sign1 * T.GetAt(i,j);  //??????
			equations((int)(size+1),i1)=0;  // RHS
			i1++;	
		}
// cleanup
		for(k1=0;k1<dim;k1++) {
			I1=T.Tail();
			T.RemoveTail(); 
			I1->newsize(0,0);
			delete I1;
		}
		} // k loop
		break;
	}

	
	double rnorme, rnorml;
	long mode;

	splt=spl[0];
	splt->progopts.resize(4);
	splt->progopts[0]=4;
	splt->progopts[1]=2; // scale all entries
	splt->progopts[2]=1;
	splt->progopts[3]=1;

	splt->coefs.resize(size,0);	//splt->coefs=0;

	// make a copy of equation just in case for wmmls (it is destroyed in lsei)
	dMat equationsC(equations);

	// attempt to solve the system quickly
	splt->m_LSEI.lsei_(equations.data(), &rows, &rowsEQ,  &rowsA, &rowsG, &size, 
		splt->progopts.data(), splt->coefs.data(), &rnorme, &rnorml, &mode);
//mode=2;
	if(mode != 0) { // in case LSEI did not find the feasible solution
//	if(mode == 2) { // in case LSEI did not find the feasible solution
		// try a more expensive method

		printf("Attempting wnnls\n");
		equations.newsize(0,0); // cleanup

		j1=size+rowsG;  // new size ( adding slack variables)
		i1=rowsEQ+rowsG; // new size of equalities system
		long j1S=j1;  // same: for passing parameters
		long i1S=i1;

		dMat equations1(j1+1,rows); // new system of eqs
		//splt->
		splt->coefs.resize(j1,0);     // here goes the solution
		//splt->coefs=0;

		splt->progopts[0]=1; // 1 no options, 6 - scale
		splt->progopts[1]=6; // 1 no options, 6 - scale

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
	ff=fopen("dump.txt","w");
	for(i=0;i<rows;i++) {
		for(j=0;j<j1+1;j++)
			fprintf(ff,"%3.1f ",equations1[j][i]);
		fprintf(ff,"\n");
	}
	fclose(ff);
*/
		splt->m_LSEI.wnnls_(equations1.data(), &rows, &i1S,  &rowsA, &j1S, &size, // coefs startring with size+1 are >0
			splt->progopts.data(), splt->coefs.data(),&rnorml, &mode);

		// clean up
		equations1.newsize(0,0);
		equationsC.newsize(0,0);
	}
	equations.newsize(0,0); // cleanup

	k=0;
//	printf("mode %d\n",mode);
/*	ff=fopen("dump.txt","a");
	for(i=0;i<spl[0]->size;i++) {
		fprintf(ff,"\n");
	for(j=0;j<spl[1]->size;j++)
		fprintf(ff,"%3.1f ",spl[0]->coefs[k++]);
	}
	fclose(ff);
*/

// translate coefs to B-spline repsn. using a=Lb
// this is slow explicit method. We can do better than this 
	k=0;

	T.clear();
	for(k=0;k<dim;k++) { // all dimensions
			splt=spl[k];
			I1=new dMat(splt->size,splt->size,0);
			//(*I1)=0;
			for(i=0;i<splt->size;i++) for(j=0;j<=i;j++) (*I1)[i][j]=1;
		// build tensor		
			T.AddTail(*I1);
	} // tensor created

	coefs.newsize(sizes); 
	splt=spl[0];
	coefs = T * splt->coefs;

// cleanup
	for(k=0;k<dim;k++) {
		I1=T.Tail();
		T.RemoveTail(); 
		I1->newsize(0,0);
		delete I1;
	}


// now find coefs of the derivatives

	// I gave up. I need a nesed lop here for all vars...
/*
	Tensor1 *cder;
	for(k=0;k<dim;k++)
	{
		splt=spl[k];
		cder=new Tensor1;
		cder->newsize(sizes);
	for(i=1;i<spl1.size;i++)
	   for(j=0;j<spl2.size;j++) {
		coefs_der[0](1,j+1)=0;
		coefs_der[0](i+1,j+1)=(coefs[i][j]-coefs[i-1][j])/
			(splt->Mknots(i+splt->order-1+splt->startindex)-splt->Mknots(i+splt->startindex));
	   }


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
	   }
*/
	   // finish here
	   // cleanup has been done earlier

	   return mode;
}

