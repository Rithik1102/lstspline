
#include "cubicsp.h"
#include <math.h>

cubic_spline::cubic_spline()
{
	m_type = 1;
    tauI = 0;
    NN = 0;
    T = 0;
    Taut = 0;
    err = 0;
}

                                     
void	cubic_spline::sortX(void)
{   
	int i,j; 
	Float_ a,b, ay,by;
	i=1;
  while(i<NN) {
	a= u[i];	b= u[i-1];
	ay= v[i]; 	by= v[i-1]; 
	
	if( a< b) {    
		j=i-1; 
		//swap i,j
		while((j>=0) && (a<b)) { 
			u[j+1]=b;  //(*x).put(b,j+1,1);
			v[j+1]=by; //(*y).put(by,j+1,1);
			j--;
			if(j>=0) {b=u[j];by=v[j];}//{ b=(*x).get(j,1); by=(*y).get(j,1); }
		}   
		if(a==b) if(j==0) a=a-0.00001; else a = a + 0.00001;         
		u[j+1]=a;  //(*x).put(a,j+1,1);
		v[j+1]=ay; //(*y).put(ay,j+1,1);
		
	} else if(a == b) u[i]=a+0.00001; 	//(*x).put(a+0.00001,i,1); 
	i++;
  }
}                                     
int		cubic_spline::BuildSpline(Float_ *x, Float_ *y, int N, int type, Float_ tau)
{
	int i;
	if (N <= 0) return -1;

    
    int work_size = 2 * N + 2;

    b.assign(work_size, 0.0);
    c.assign(work_size, 0.0);
    h.assign(work_size, 0.0);
    u.assign(work_size, 0.0);
    v.assign(work_size, 0.0);
    alpha.assign(work_size, 0.0);
    beta.assign(work_size, 0.0);
    gamma.assign(work_size, 0.0);

    m_type=type;

    Taut=(tau>0.01);
    tauI=tau;       
    NN=N;  
    err = 0;
    T = 0;  
	for(i=0;i<N;i++) {  
		u[i]=x[i]; v[i]=y[i];       
	}	
	                    
//	sort 
	sortX();
	if(Taut && N>4) m_type=2;
	if(N<4) m_type=0;
	switch(m_type) {
	case 2: TautSpline(u.data(),v.data(),NN); break;
	case 0:
	case 1: CubSpline(u.data(),v.data(),NN);  break;   
	case 3: BuildMonotonSpline(u.data(),v.data(),NN);  break;   
	}
	return err;		
}   

Float_  cubic_spline::SplineValue( Float_ t)
{   Float_ val;
	switch(m_type) {
	case 2:  val=TautSplineval(u.data(),v.data(),t); break;
	case 0:  val=LinearSpline(u.data(),v.data(),t);  break; 
	case 1:  val=Splineval(u.data(),v.data(),t); break;
	case 3:  val=MonotoneSplineValue(t); break;
	}
	return val;	
}

Float_  cubic_spline::LinearSpline(Float_ *x, Float_ *y, Float_ t)  
{	
	int i;      
	if (NN<1) return 0; else if (NN<2) return y[0]; else
    Bisection(t,x, 1, NN, &i);
    return(c[i-1]*t+h[i-1]);
}
                                                                   
// Procedure to find the number of the interval where the point u is situated.
// Parameters: u - the point of interest, x the vector of knots, ileft,iright - the boundaries, l - the 
// solution.
void cubic_spline::Bisection(Float_ u, Float_ *x, int ileft,int iright,int *l)
{
 int itemp,mid,temp;
 temp=ileft; itemp=iright;
 while ((itemp-temp)>1 ) {
     mid = (int) floor((temp+itemp) /2);
     if (u<x[mid-1]) itemp=mid; else temp=mid;
 }
 *l=temp;
}


// Procedure to compute the coefficients of cubic interpolating spline.
// Parameters: x - vector of knots, y - vector of function values in the knots, N -number of knots, 
// h - the computed distance between the knots at the output, c - vector of spline coefficients.
void cubic_spline::CubSpline(Float_ *x, Float_ *y, int N)
{
 	int i; 
 	if (N<2) return;
	if (N<4) {
		// linear spline,2 or 3 points  c[i]*x+h[i]
		c[0]=(y[0]-y[1])/(x[0]-x[1]);
		h[0]=y[0]-c[0]*x[0];
		if(N==3) {          
			c[1]=(y[1]-y[2])/(x[1]-x[2]);
			h[1]=y[1]-c[1]*x[1];
		}
	return;}
	
	for( i=0;i < N-1;i++) {
      h[i]=(x[i+1]-x[i]);
      b[i]=(6*(y[i+1]-y[i])/h[i]);  
     } 
 	alpha[2 -1]=2*(h[1 -1]+h[2 -1]); 
 	beta[2 -1]=b[2 -1]-b[1 -1];
 	for( i=2;i <  N-1;i++) {
       alpha[i]=(2*(h[i]+h[i-1]) - h[i-1]*h[i-1]/alpha[i-1]);
       beta[i]=(b[i]-b[i-1]-h[i-1]*beta[i-1]/alpha[i-1]);
    }                        
    
 	c[N -1]=0;
 	for( i=N-1; i>=2; i--){ // from N-1 to 2 by -1 do
     	c[i -1]=((beta[i -1]-h[i -1]*c[i+1 -1])/alpha[i -1]);
     }	
 	c[0]=0;
} 	
 
// Procedure to compute the value of spline in the given point.
// Parameters: x -vector of knots, y -vector of funcion values in the knots, N- number of knots, h 
// - vector of distances between the knots, c - vector of spline coefficients, t - the point of 
// interest. As the outpunt the value of spline is returned.
Float_ cubic_spline::Splineval(Float_ *x, Float_ *y, Float_ t)
{  
	Float_ g, A,B,C,u;
	int i;
    Bisection(t,x, 1, NN, &i);
    A=((c[i+1 -1]-c[i -1])/(6*h[i -1])); 
    B=c[i -1]/2.0;   
    C=(-h[i -1]*c[i+1 -1]/6-h[i -1]*c[i -1]/3+(y[i+1 -1]-y[i -1])/h[i -1]);
    u=t-x[i -1];
    g=y[i -1]+u*(C+u*(B+u*A)); 
    return g;
}
 


// Procedure to compute the coefficient of exponential spline in tension.
// Parameters: x - vector of knots, y - vector of function values in the knots, N -number of knots, 
// tau - tension parameter bigger than 0 (actually bigger than 0.01), h - the computed distance 
// between the knots at the output, c - vector of spline coefficients.
void cubic_spline::TautSpline(Float_ *x, Float_ *y, int N)
{	
	int	i;
 	for (i=1;i<= N-1; i++) {
      h[i -1]=(x[i+1 -1]-x[i -1]);
      gamma[i -1]=( tauI*tauI * (y[i+1 -1]-y[i -1]) / h[i -1]);
      alpha[i -1]=(1/h[i -1]-  tauI / sinh(tauI*h[i -1]));
      beta[i -1]=(tauI* cosh(tauI*h[i -1]) / sinh(tauI*h[i -1]) - 1/h[i -1]);
 	}
 
 	u[2 -1]=(beta[1 -1]+beta[2 -1]); 
 	v[2 -1]=gamma[2 -1]-gamma[1 -1];
 	for (i=3;i<= N-1;i++) {
       u[i -1]=((beta[i -1]+beta[i-1 -1])- alpha[i-1 -1]*alpha[i-1 -1]/u[i-1 -1]);
       v[i -1]=(gamma[i -1]-gamma[i-1 -1]  - alpha[i-1 -1] * v[i-1 -1]/u[i-1 -1]);
 	}
 	c[N -1]=0;
 	for( i= N-1; i>=2; i--)
     c[i -1]=((v[i -1]-alpha[i -1]*c[i+1 -1])/u[i -1]);
 
 	c[1 -1]=0;
}
 
// 
// Procedure to compute the value of spline in tension in the given point.
// Parameters: x -vector of knots, y -vector of funcion values in the knots, N- number of knots, 
// tau- the tension parameter, h - vector of distances between the knots, c - vector of spline 
// coefficients, t - the point of interest. As the outpunt the value of spline is returned.
// 
Float_ cubic_spline::TautSplineval(Float_ *x, Float_ *y, Float_ t)
{
 	int i;
 	Float_ g,z;

    Bisection(t,x,1,NN,&i);
    z=(x[i+1 -1]-t)/h[i -1]; g=h[i -1]*tauI;
    
    g=(c[i -1] * sinh(g*z) + c[i+1 -1]*sinh(g*(1-z))  )/(tauI*tauI*sinh(g))  + 
         (y[i -1] - c[i -1]/tauI/tauI)*(z) + (y[i+1 -1] -c[i+1 -1]/tauI/tauI)*(1-z);   
    return g;
}


// monotone spline

void cubic_spline::CalculateBeta()
{
	int i=0;
	while(i<=(NN-2))
	{
		b[i]=(v[i+1] - v[i])/(u[i+1] - u[i]);
		i++;
	};

}

void cubic_spline::CalculateDeltaGeneral()
{
	int i=1;
	while(i<=(NN-2))
	{
		if((b[i-1]*b[i])<=0)
			c[i]=0;
		else
			c[i]=(2*b[i-1]*b[i])/(b[i-1]+b[i]);
		i++;
	}
}


void cubic_spline::CalculateDeltaFirst()
{
	if((b[0]*(2*b[0]-c[1]))<=0)
		c[0]=0;
	else
		c[0]=2*b[0] - c[1];
}
void cubic_spline::CalculateDeltaLast()
{
	if((b[NN-2]*(2*b[NN-2]-c[NN-2]))<=0)
		c[NN-1]=0;
	else
		c[NN-1]=2*b[NN-2] - c[NN-2];
}


void cubic_spline::CalculateCoefficientsKnots()
{
	int i=0,s=0,ch=0;
	while(i<=(NN-2))
	{
		//decide whether an additional knot is necessary
		if(fabs(c[i]+c[i+1]- 2*b[i])<=0.1e-5)
		{
			//no additional knot 
			h[s]=u[i];
			alpha[s]=v[i];
			beta[s]=c[i];
			gamma[s]=(c[i+1]-c[i])/(2*(u[i+1]-u[i]));
			s++;
		}
		else
		{
			//adding a knot
			h[s]=u[i];
			//determine the position of the knot
			if((c[i+1] - b[i])*(c[i] - b[i])<0)
				h[s+1]=u[i+1] + (c[i] - b[i])*(u[i+1]-u[i])/(c[i+1] - c[i]);
			else
				h[s+1]=0.5*(u[i+1] + u[i]);

			//calculate coefficients
			double dtemp = (2*b[i] - c[i+1]) + ((c[i+1] - c[i])*(h[s+1] - u[i]))/
												(u[i+1] - u[i]);
			alpha[s]=v[i];
			beta[s]=c[i];
			gamma[s]=(dtemp - c[i])/(2*(h[s+1] - u[i]));
			alpha[s+1]=v[i] + c[i]*(h[s+1] - u[i]) + (dtemp - c[i])*(h[s+1] - u[i])/2;
			gamma[s+1]=(c[i+1] - dtemp)/(2*(u[i+1] - h[s+1]));
			beta[s+1]=dtemp;

			s++;
			s++;
		}
		i++;
	}
	h[s]=u[NN-1];
	T=s+1;
}
				
int cubic_spline::BuildMonotonSpline(Float_ *xx, Float_ *yy, int NN)
{
	CalculateBeta();
	CalculateDeltaGeneral();
	CalculateDeltaFirst();
	CalculateDeltaLast();
	CalculateCoefficientsKnots();
	return(1);
}


Float_ cubic_spline::MonotoneSplineValue(Float_ t)
{
	double value;
	int i;
	Bisection(t,h.data(),1,T,&i); i--;
	if(i<0)i=0;
	value=(t-h[i]);
	value=alpha[i] + value*(beta[i] + gamma[i]*value);
	return value;
};
