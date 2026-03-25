#define Float_ double
#include <vector>
class cubic_spline {
	public:  
		int		m_type; // 1 = cubic, 2 = taught, 3 - monotone
	 	
		std::vector<Float_> b, c, h, u, v, alpha, beta, gamma;

	 	Float_ tauI; 
	 	int NN,T;          
	 	int	Taut, err;
		void Bisection(Float_ u, Float_ *x, int ileft,int iright,int *l);
		void CubSpline(Float_ *x, Float_ *y, int N);
		Float_ Splineval(Float_ *x, Float_ *y, Float_ t);
		void TautSpline(Float_ *x, Float_ *y, int N);
		Float_ TautSplineval(Float_ *x, Float_ *y, Float_ t);   
		Float_ LinearSpline(Float_ *x, Float_ *y, Float_ t); 
		void	sortX(void);
		// monotone spline
		void CalculateBeta();
		void CalculateDeltaGeneral();
		void CalculateDeltaFirst();
		void CalculateDeltaLast();
		void CalculateCoefficientsKnots();
		int  BuildMonotonSpline(Float_ *xx, Float_ *yy, int NN);
		Float_  MonotoneSplineValue( Float_ t);
	public:    
		cubic_spline();
		int		BuildSpline(Float_ *x, Float_ *y, int N, int type, Float_ tau);
		Float_  SplineValue( Float_ t);
}; 
/*
private:

	void Bisection(Float_ u, int *i);

public:
	monoton_spline();
	virtual ~monoton_spline();
*/
