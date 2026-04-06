
//#include "tnt/tnt.h"
//#include "bspline.h"
#include<vector>

#include "nspline.h"
#include "lsnsplinet.h"

#include "math.h"

//#include "tnt/vec.h"
//#include "tnt/cmat.h"

using namespace std;

#define sqr(a) ((a)*(a))

//#define DATATYPE double

template<class DATATYPE>
DATATYPE TheRand(DATATYPE left=0, DATATYPE right=1)
{	    
//	return y + (x*rand())/(RAND_MAX+y);
	return left +  (right-left)* rand()/double(RAND_MAX);
}

int splineXdtest(char* inp)
{
	FILE	*fp;

	int dim,type;
	dim =3; type = 111;

	int i,j,k,k1;
	VecArray knots,data,exactdata;
	dVec *v1;
	double r;

	for(k=0;k<dim;k++) {
		k1=4;
		v1=new dVec;
		v1->resize(k1);
		for(i=0;i<k1;i++) (*v1)[i]=i*2.0/k1;
		knots.push_back(v1);
	}

	k1=60;
	srand(11);
	double fff;

	for(i=0;i<k1;i++)
	{
		v1=new dVec;
		v1->resize(dim+1);
		fff=1;
		for(j=0;j<dim;j++) {
			r=TheRand(0.0,2.0);
			(*v1)[j]=r;
			fff *= r;
		}

		(*v1)[dim]=fff*(1+TheRand(-1.0,1.0)*0.0001);

		data.push_back(v1);

	}

	k1=0;
// exact data
	v1=new dVec;
	v1->resize(dim+1);
		fff=1;
		for(j=0;j<dim;j++) {
			r=0;
		  (*v1)[j]=r;
			fff *= r;
		}
		(*v1)[dim]=fff;

		exactdata.push_back(v1);

	v1=new dVec;
	v1->resize(dim+1);
		fff=1;
		for(j=0;j<dim;j++) {
			r=2;
		  (*v1)[j]=r;
			fff *= r;
		}
		(*v1)[dim]=fff;

		exactdata.push_back(v1);

	for(i=0;i<k1;i++)
	{
		v1=new dVec;
		v1->resize(dim+1);

		fff=1;
		for(j=0;j<dim;j++) {
			r=TheRand(0,2);
		  (*v1)[j]=r;
			fff *= r;
		}
		(*v1)[dim]=fff;

		exactdata.push_back(v1);
	}

//	fp=fopen(inp,"w");
//	fclose(fp);

	LSNSplineT s;

	int mode=s.Build(type,dim,knots,data,exactdata);
	printf("Exit mode: %d\n",mode);

	double ee,err=0;

	for(i=0;i<10;i++)
	{
		v1=data[i];
		for(j=0;j<=dim;j++) printf("%f ",(*v1)[j]); printf("    ");
		ee=s.Value(*v1);   printf("%f ",ee);
		ee=fabs(ee-(*v1)[dim]);
		printf("\n");
		err +=ee;

	}
	printf("error: %f\n",err);

	for(i=0;i<2;i++)
	{
		v1=exactdata[i];
		for(j=0;j<=dim;j++) printf("%f ",(*v1)[j]); printf("    ");
		ee=s.Value(*v1);   printf("%f ",ee);
		ee=fabs(ee-(*v1)[dim]);
		printf("\n");
		err +=ee;

	}
	printf("error: %f\n",err);

	return 0;
}


int splineXd(char* inp)
{
	FILE	*fp;
	fp=fopen(inp,"r");
	int dim,type;

	fscanf(fp,"dim %d type %d\n",&dim,&type);
	int i,j,k,k1;
	VecArray knots,data,exactdata;
	dVec *v1;
	double r;

	for(k=0;k<dim;k++) {
		fscanf(fp,"knots %d\n",&k1);
		v1=new dVec;
		v1->resize(k1);
		for(i=0;i<k1;i++) (*v1)[i]=1.0*i/k1;
		knots.push_back(v1);
	}

	fscanf(fp,"data %d\n",&k1);
	for(i=0;i<k1;i++)
	{
		v1=new dVec;
		v1->resize(dim+1);

		for(j=0;j<dim;j++) {
		  fscanf(fp,"%lf ",&r);
		  (*v1)[j]=r;
		}
		fscanf(fp,"%lf\n",&r);
		(*v1)[dim]=r;

		data.push_back(v1);

	}

	fscanf(fp,"exact data %d\n",&k1);
	for(i=0;i<k1;i++)
	{
		v1=new dVec;
		v1->resize(dim+1);

		for(j=0;j<dim;j++) {
		  fscanf(fp,"%lf ",&r);
		  (*v1)[j]=r;
		}
		fscanf(fp,"%lf\n",&r);
		(*v1)[dim]=r;

		exactdata.push_back(v1);
	}

	fclose(fp);

	LSNSplineT s;

	int mode=s.Build(type,dim,knots,data,exactdata);
	printf("Exit mode: %d\n",mode);


	dVec t(dim, 0.5);
	//t=0.5;

	s.Value(t);


	return 0;
}


int spline2Xd(char* inp)
// Build 2d spline
{
	FILE	*fp;

	int i,j,k1, k2,M1, M2,type;
// read the data
	fp=fopen(inp,"r");
	fscanf(fp,"knots %d %d order %d %d type%d\n",&k1,&k2, &M1, &M2, &type);

	dVec t1(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf ",&(t1[i]));
	fscanf(fp,"\n");

	dVec t2(k2);
	for(i=0;i<k2;i++) fscanf(fp,"%lf ",&(t2[i]));
	fscanf(fp,"\n");

	fscanf(fp,"data %d\n",&k1);
	dMat x(k1,2);
	dVec y(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf %lf %lf\n",&(x[i][0]),&(x[i][1]),&(y[i]));

	fscanf(fp,"exact data %d\n",&k1);
	dMat xe(k1,2);
	dVec ye(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf %lf %lf\n",&(xe[i][0]),&(xe[i][1]),&(ye[i]));

	VecArray knots,data,exactdata;
	dVec *v1;

	int dim=2;
	knots.push_back(&t1);
	knots.push_back(&t2);

	for(i=0;i<y.size();i++)
	{
		v1=new dVec;
		v1->resize(dim+1);

		for(j=0;j<dim;j++) {
		  (*v1)[j]=x[i][j];
		}
		(*v1)[dim]=y[i];

		data.push_back(v1);
	}

	for(i=0;i<ye.size();i++)
	{
		v1=new dVec;
		v1->resize(dim+1);

		for(j=0;j<dim;j++) {
		  (*v1)[j]=xe[i][j];
		}
		(*v1)[dim]=ye[i];

		exactdata.push_back(v1);
	}

	LSNSplineT s;

	fclose(fp);


	int Ndata=x.dim(1);
	int Nknots1=t1.size();
	int Nknots2=t2.size();
	int NEdata=xe.dim(1);

	int r=s.Build(type,dim,knots,data,exactdata);
	printf("Exit mode: %d\n",r);


//	fp=fopen("c:\\maplev4\\bin.win\\output2.txt","w");
	fp=fopen("output2.txt","w");
	fprintf(fp,"%d\n",Ndata);
	for(i=0;i<Ndata;i++) { 
		fprintf(fp,"%f %f %f\n",x[i][0],x[i][1],y[i]);
	}

	int PTS=20; // in each dim
	double h1=(t1[t1.size()-1]-t1[0])/(PTS-1.);
	double h2=(t2[t2.size()-1]-t2[0])/(PTS-1.);
	double tt1, tt2;
	double f;
	fprintf(fp,"%d\n",PTS);
	tt1=t1[0];
	dVec tt(dim);

	for(i=0;i<PTS;i++) {
	tt2=t2[0];

	for(j=0;j<PTS;j++) {
		tt[0]=tt1; tt[1]=tt2;
		f=s.Value(tt);
		fprintf(fp,"%f %f %f\n",tt1,tt2,f);
		tt2+=h2;
	}
	tt1+=h1;
	}
// derivatives

	fclose(fp);

	return 0;
}


int spline2d(char* inp)
// Build 2d spline
{
	FILE	*fp;

	int i,j,k1, k2,M1, M2,type;
// read the data
	fp=fopen(inp,"r");
	fscanf(fp,"knots %d %d order %d %d type%d\n",&k1,&k2, &M1, &M2, &type);

	dVec t1(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf ",&(t1[i]));
	fscanf(fp,"\n");

	dVec t2(k2);
	for(i=0;i<k2;i++) fscanf(fp,"%lf ",&(t2[i]));
	fscanf(fp,"\n");

	fscanf(fp,"data %d\n",&k1);
	dMat x(k1,2);
	dVec y(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf %lf %lf\n",&(x[i][0]),&(x[i][1]),&(y[i]));

	fscanf(fp,"exact data %d\n",&k1);
	dMat xe(k1,2);
	dVec ye(k1);
	for(i=0;i<k1;i++) fscanf(fp,"%lf %lf %lf\n",&(xe[i][0]),&(xe[i][1]),&(ye[i]));

	fscanf(fp,"bound # %d %d\n",&j,&k1);
	double bou;

	LSNSpline2 s;

	if(j==1) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundDownX(bou);
	}
	if(j==2) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundUpX(bou);
	}
	if(j==3) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundDownX(bou);
		fscanf(fp,"%lf\n",&bou);
		s.BoundUpX(bou);
	}
	if(k1==1) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundDownY(bou);
	}
	if(k1==2) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundUpY(bou);
	}
	if(k1==3) {
		fscanf(fp,"%lf\n",&bou);
		s.BoundDownY(bou);
		fscanf(fp,"%lf\n",&bou);
		s.BoundUpY(bou);
	}


	fclose(fp);


	int Ndata=x.dim(1);
	int Nknots1=t1.size();
	int Nknots2=t2.size();
	int NEdata=xe.dim(1);

	int r=s.Build(type,M1,M2,t1,t2,x,y,xe,ye);
	printf("Exit mode: %d\n",r);


	fp=fopen("output2.txt","w");
	fprintf(fp,"%d\n",Ndata);
	for(i=0;i<Ndata;i++) { 
		fprintf(fp,"%f %f %f\n",x[i][0],x[i][1],y[i]);
	}

	int PTS=20; // in each dim
	double h1=(t1[t1.size()-1]-t1[0])/(PTS-1.);
	double h2=(t2[t2.size()-1]-t2[0])/(PTS-1.);
	double tt1, tt2;
	double f;
	fprintf(fp,"%d\n",PTS);
	tt1=t1[0];
	for(i=0;i<PTS;i++) {
	tt2=t2[0];
	for(j=0;j<PTS;j++) {
		f=s.Value(tt1,tt2);
		fprintf(fp,"%f %f %f\n",tt1,tt2,f);
		tt2+=h2;
	}
	tt1+=h1;
	}
// derivatives
	fprintf(fp,"%d\n",PTS);
	tt1=t1[0];
	for(i=0;i<PTS;i++) {
	tt2=t2[0];
	for(j=0;j<PTS;j++) {
		f=s.ValueDer(tt1,tt2,0);
		fprintf(fp,"%f %f %f\n",tt1,tt2,f);
		tt2+=h2;
	}
	tt1+=h1;
	}
	//y
	fprintf(fp,"%d\n",PTS);
	tt1=t1[0];
	for(i=0;i<PTS;i++) {
	tt2=t2[0];
	for(j=0;j<PTS;j++) {
		f=s.ValueDer(tt1,tt2,1);
		fprintf(fp,"%f %f %f\n",tt1,tt2,f);
		tt2+=h2;
	}
	tt1+=h1;
	}

	fclose(fp);

	return 0;
}

int spline1d(char* inp)
{
	FILE	*fp;

	int i,j,k,M,type;
// read the data
	fp=fopen(inp,"r");
	fscanf(fp,"knots %d order %d type %d\n",&k,&M,&type);

	cout << k << M << type << endl;
	dVec t(k);
	for(i=0;i<k;i++) fscanf(fp,"%lf ",&(t[i]));
	fscanf(fp,"\n");

	fscanf(fp,"data %d\n",&k);
	dVec x(k),y(k);
	for(i=0;i<k;i++) fscanf(fp,"%lf %lf\n",&(x[i]),&(y[i]));

	fscanf(fp,"exact data %d\n",&k);
	dVec xe(k),ye(k);
	for(i=0;i<k;i++) fscanf(fp,"%lf %lf\n",&(xe[i]),&(ye[i]));

	fscanf(fp,"bound # %d %d\n",&j,&k);

	dVec bound1(k),bound2(k);
	switch(j) {
	case 1:
		for(i=0;i<k;i++) fscanf(fp,"%lf\n",&(bound1[i])); break;
	case 2:
		for(i=0;i<k;i++) fscanf(fp,"%lf\n",&(bound2[i])); break;
	case 3:
		for(i=0;i<k;i++) fscanf(fp,"%lf %lf\n",&(bound1[i]),&(bound2[i])); break;
	case 0:;
	}

	fclose(fp);

	int Ndata=x.size();
	int Nknots=t.size();
	int NEdata=xe.size();

//

	LSNSpline s;


//
	// only 1 scalar value
	if(k==1) {
	switch(j) {
	case 1:  s.Boundlow(bound1[0]); break;
	case 2:  s.Boundup(bound2[0]); break;
	case 3:  s.Bounds(bound1[0],bound2[0]); break;
	}
	} else {
	switch(j) {
	case 1:  s.Boundlow(bound1); break;
	case 2:  s.Boundup(bound2); break;
	case 3:  s.Bounds(bound1,bound2); break;
	}
	} // if

	cout << type << endl;

	int r=s.Build(type,M,t,x,y,xe,ye);
	printf("Exit mode: %d\n",r);


	fp=fopen("output1.txt","w");
	fprintf(fp,"%d\n",Ndata);
	for(i=0;i<Ndata;i++) { 
		fprintf(fp,"%f %f\n",x[i],y[i]);
	}

	int PTS=300;
	double h=(x[Ndata-1]-x[0])/(PTS-1.);
	double tt=x[0];
	double f;
	fprintf(fp,"%d\n",PTS);
	for(i=0;i<PTS;i++) {
//		f=s.Value(tt);
		f=s.ValueDer(tt,0);
		fprintf(fp,"%f %f\n",tt,f);
		tt+=h;
	}
	tt=x[0];
	fprintf(fp,"%d\n",PTS);
	for(i=0;i<PTS;i++) {
		f=s.ValueDer(tt,1);
		fprintf(fp,"%f %f\n",tt,f);
		tt+=h;
	}

	fclose(fp);
	

        return 0;
}


#include <span>
int main()
{
	char s[80];

	//test
	VecArray t;
	dMat A(4, 5, 11);

	A[1][0] = 12;



	A.MakeVecArray(t);
	cout << t.size() << endl;

//	cout <<  (*(span<double>*)((t[0]))) [1] << endl;



	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++)
			//cout << (*(span<double>*)((t[i])))[j] << " ";
			cout << (*(t[i]))[j] << " ";
		cout << endl;
	}

	for (auto it = t[1]->begin(); it != t[1]->end(); ++it)
	{
		cout << *it;
	}
	cout << endl;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++)
			cout << *((t[i])->data() +j) << " ";
		cout << endl;
	}

	cout << (*(span<double>*)((t[1]))).size()<<endl;
	cout << (*reinterpret_cast<span<double>*>((t[1]))).size() << endl;

	// no write access
//	(*(span<double>*)((t[0])))[4] = 13;
	A[1][0] = 14;	
	deleteVecArray(t);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}





L1:
	printf("Enter 1 for 1-d spline and 2 for 2-d spline, 3 for 2dtensor 4 for n-d tensor 0 to exit:");
	int i;
	scanf("%d",&i);
	if(i==0) return 0;

	printf("\nInput file name:");
	scanf("%s",s);

	switch(i) {
	case 1:
		spline1d(s); break;
	case 2:
		spline2d(s); break;
	case 3:
		spline2Xd(s); break;
	case 4:
		splineXd(s); break;
	case 5:
		splineXdtest(s); break;
	}
	//"spldata2.txt"
	goto L1;
}

/*
    Vector<double> x(Nx), y(Nx);
    Matrix<double> A(N,N);
    Vector<double> xe(Nxe), ye(Nxe);
    Vector<double> bound1(N), bound2(N);

	dVec	r;
	iVec	inx;
	
	dVec	t(N);

//	t[0]=0; t[1]=1.3; t[2]=1.5; t[3]=2.; t[4]=3.7; t[5]=3.7; t[6]=3.9; t[7]=4.5;
//	t[8]= 4.9; t[9]=5.2; t[10]=5.8; t[11]=6.5;
	t[0]=0; t[1]=2; t[2]=3; t[3]=5; t[4]=6; t[5]=8; t[6]=9; t[7]=11;
	t[8]= 11; t[9]=12; t[10]=14; t[11]=15;

	t[0]=0; t[1]=2; t[2]=3; t[3]=5; t[4]=6; t[5]=8; t[6]=9; t[7]=10;
	t[8]= 11; t[9]=11.5; t[10]=12; t[11]=14; t[12]=15;
	
//	t[0]=0; t[1]=0.5; t[2]=1; t[3]=1.5; t[4]=2; t[5]=2.5; t[6]=3.; t[7]=3.5;
//	t[8]= 4.0; t[9]=4.5; t[10]=5.0; t[11]=5.5;
	for(i=0;i<Nx;i++) { 
		x[i]=t[i];
//		x[i]=i*t(N)/(Nx-1.); 
//		y[i]=fun(x[i]);
		 }
//	t[7]=t[8]=t[9]=11.1; t[10]=11.1;

	y[0]=10; y[1]=10; y[2]=10; y[3]=10; y[4]=10; y[5]=10; y[6]=10.5; y[7]=10.5;
	y[8]= 12; y[9]=12; y[10]=50; y[11]=53; y[12]=85;

//	t.newsize(6);
//	t[0]=0; t[1]=3; t[2]=10; t[3]=11; t[4]=12; t[5]=14; 

	for(i=0;i<Nxe;i++) { xe[i]=(Nxe<2? 0 : i*t(N)/(Nxe-1.)); 
		ye[i]=fun(xe[i]);
	}
	for(i=0;i<N;i++) { bound1[i]=-20.0; 
		bound2[i]=40;
	}
	bound1[5]=-20;

	m[3]=1;
	m[4]=2;
//	s.Init(M,t);


*/

/*
#include "titan.h"

void main()
{
	FILE	*fp;
    const int N=7;
    const int Nx=49;
    Vector<double> x(Nx), y(Nx,gtitan);
    Vector<double> xe(0), ye(0);
	int i,j;

//	for(i=0;i<Nx;i++) x[i]=585.+10*i;
	for(i=0;i<Nx;i++) x[i]=0.+75./(49.-1)*i;

	dVec	r;
	iVec	inx;
	
	dVec	t(N);

	double h;
	h=(x[Nx-1]-x[0])/(N-1.0+0.0);
	for(i=0;i<N;i++) t[i]=x[0]+h*i;

	t[0]=0; t[1]=37.62; t[2]=43.96; t[3]=47.43; t[4]=50.14; t[5]=59.34; t[6]=75;
//	t[0]=0; t[1]=37.65; t[2]=43.97; t[3]=47.37; t[4]=50.21; t[5]=59.2; t[6]=75;
//	t[0]=0; t[1]=14.69; t[2]=37.91; t[3]=44.2; t[4]=47.28; t[5]=48.93; t[6]=75;

	LSNSplineB s;

	int M=4;

	s.Build(0,M,t,x,y,xe,ye);

//	fp=fopen("c:\\maplev4\\bin.win\\output1.txt","w");
	fp=fopen("d:\\maplev4\\bin.win\\output1.txt","w");
	fprintf(fp,"%d\n",Nx);
	for(i=0;i<Nx;i++) { 
		fprintf(fp,"%f %f\n",x[i],y[i]);
	}

	double f;
//	f=s.ValueDer(t[N-1],0);
//	f=s.ValueDer(t[N-1]+1,0);

	int PTS=300;
	h=(t[N-1]-t[0])/(PTS -1);
	double tt=t[0];
	fprintf(fp,"%d\n",PTS);
	for(i=0;i<PTS;i++) {
//		f=s.Value(tt);
		f=s.ValueDer(tt,0);
		fprintf(fp,"%f %f\n",tt,f);
		tt+=h;
	}

	fclose(fp);

	double delta=0;

	for(i=1;i<Nx-1;i++) { 
		f=s.Value(x[i]);
		delta += sqr((f-y[i]));
	}
	f=s.Value(x[0]);
	delta += sqr((f-y[i]))/4;
	f=s.Value(x[Nx-1]);
	delta += sqr((f-y[i]))/4;

	for(i=-3;i<N+3;i++)
		printf("%f ",s.Mknots(i));
	printf("\n");
	//delta=sqrt(delta);
	printf("%f\n",delta);
	

}
*/
