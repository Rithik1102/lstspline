#ifndef TENSOR_H
#define TENSOR_H


#include <stdlib.h>
#include <vector>

using namespace std;


//typedef TNT::Matrix<double>		dMat;
//typedef TNT::Vector<int>		iVec;
//typedef TNT::Vector<double>		dVec;

typedef vector<dMat*> ArrayOfPointers;
typedef vector<dVec*> ArrayOfPointersV;

typedef int Subscript;
 
class Tensor {
// n-dimensional tensor (actually 2n-dimensional)
// components are pointers to matrices

private:
	ArrayOfPointers	Comp;

public:
	Tensor() {};
	~Tensor() {};

	void AddTail(dMat& M) {Comp.push_back( &M );} 
	void RemoveTail() {Comp.pop_back(  );} 
	dMat* Tail() {return Comp.back(); }
	void	clear() {Comp.clear();}

	int dim(int i)
	{ // returns the size of the matrix i=1 in rows, i=2 in columns
		int s=1,j;
		for (j=0;j<Comp.size();j++) 
			s *= Comp[j]->dim(i);
		return s;
	}

	void SetMatrixAt(int i, dMat& M) {Comp[i] = (&M);} 

    double operator()(iVec& index)
	{
		int i,j;
		double f=1;
		for(i=0,j=0; i<index.size() && j<Comp.size(); i+=2,j++)
			f *=(*Comp[j])(index[i]+1,index[i+1]+1);

		return f;
	}

    double operator()(Subscript i, Subscript j, Subscript k, Subscript l)
	{
		return (*Comp[0])(i,j) * (*Comp[1])(k,l);
	}
    double operator()(Subscript i, Subscript j, Subscript k, Subscript l) const
	{
		return (*Comp[0])(i,j) * (*Comp[1])(k,l);
	}

    double operator()(Subscript i, Subscript j)
    { 
		div_t d1,d2;
		d1=div(i-1,(*Comp[1]).dim(1));
		d2=div(j-1,(*Comp[1]).dim(2));
        return  (*this)(d1.quot+1,d2.quot+1,d1.rem+1,d2.rem+1); 
    }


    
    double operator() (Subscript i, Subscript j) const
    {
		div_t d1,d2;
		d1=div(i,(*Comp[1]).dim(1));
		d2=div(j,(*Comp[1]).dim(2));
        return  (*this)(d1.quot,d2.quot,d1.rem,d2.rem); 
    }

	double GetAt(Subscript i, Subscript j)
	{
		div_t d1,d2;
		i--; j--;
		iVec index(Comp.size()*2);
		int j1,i1,I,J;
		I=1; J=1;
		for(i1=1;i1<Comp.size();i1++) {I *= (*Comp[i1]).dim(1);
									   J *= (*Comp[i1]).dim(2);}
		j1=0;
		for(i1=1;i1<Comp.size();i1++) {
			d1=div(i, I);
			d2=div(j, J);
			index[j1]=d1.quot; j1++;
			index[j1]=d2.quot; j1++;
			i=d1.rem; j=d2.rem;
			I /=(*Comp[i1]).dim(1);
			J /=(*Comp[i1]).dim(2);
		}
		// last ones
		index[j1]=i; j1++;
		index[j1]=j; 
		return (*this)(index);
	}

	dVec operator*( const dVec& r)
	{
		dVec tmp(dim(1));
		double tmp1;
		int i,j;
		for(i=1;i<=dim(1);i++){
			tmp1=0;
			for(j=1;j<=dim(2);j++) tmp1 += r[j-1]*GetAt(i,j);
			tmp[i-1]=tmp1;
		}
		return tmp;
	}

};

// just a n-d tensor product of vectors, a tool for calculating the indices
class Tensor1:public dVec {
public:

	iVec sizes;
	void newsize(iVec& index) 
	{
		int total=1;
		int i;
		for(i=0;i<index.size();i++) total *= index[i];
		dVec::resize(total);
		sizes=index;
	}
//	dVec data;
//	Tensor1() {};
//	Tensor1(int n) {data.newsize(n);}
//	~Tensor1() {data.newsize(0);}

 //   double operator()(int i) {return data(i);}
 //   double operator()(int i) const {return data(i);}

    double operator()(iVec& index)
	{
		int i,p;
		p=index[0]-1; // they were 1-based
		for(i=2;i<=index.size();i++)
			p = p * sizes[i-1] + index[i-1]-1;

		return  at(p);
	}
	Tensor1& operator=(const dVec &A)
	{
		assign(A.begin(), A.end());
		//copy(A.begin(), A.end(), this);
/*		(*this)(A);

		if (this->size() == A.size())         // no need to re-alloc
            copy(A.begin());
       else
        {
            destroy();
            initialize(A.dim());
            copy(A.begin());
        }*/
        return *this;
	}

};

#endif