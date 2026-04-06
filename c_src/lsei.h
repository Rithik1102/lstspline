/* 587.f -- translated by f2c (version 20000704).
   You must link the resulting object file with the libraries:
    -lf2c -lm   (in that order)
*/

#ifndef LSEI_H
#define LSEI_H

#include <stdlib.h>
#include <math.h>

// lsei_error.h: interface for the err_message class.
//
//////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <vector>

typedef long int integer ;
typedef long int logical;
typedef double doublereal ;
#define FALSE_ 0
#define TRUE_ 1

#ifndef _WINDEF_
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

#define dmax(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) ((a) <= (b) ? (a) : (b))
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)fabs(x)

using namespace std;
//const int N = 10;

class err_message  
{
public:	
	err_message();
	err_message(const char *messge, integer *nerror, integer *errlevel);
	virtual ~err_message();

	const char* Getmessage();
	int GetErrLevel();
	int GetErrNo();

private:
	int level;
	int nerr;
	string messg;
};


class err_list  
{
public:
	err_list();
	virtual ~err_list();

	void printErrors();
	bool insert(err_message& n);
	
private:
	vector<err_message> v;
};



class lsei {
private:
    err_list errors;
    double M_precision;
public:
    /*constructor*/
    int lout;
    lsei(){ M_precision=0; }
    /*destructor*/
    ~lsei(){;}

/*     SUBROUTINE LSEI(W,MDW,ME,MA,MG,N,PRGOPT,X,RNORME,RNORML,MODE,WS, IP) */
/*     ALGORITHM 587, COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL. 8, NO. 3, */
/*     SEP., 1982, P.323. */
/*     REVISED OCT. 1, 1981. */
private:
    /* Subroutine */ int lsei_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *mg, integer *n, double *prgopt, double *x, double *rnorme, double *
    rnorml, integer *mode, double *ws, integer *ip);

public:
    /* Subroutine */ int lsei_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *mg, integer *n, double *prgopt, double *x, double *rnorme, double *
    rnorml, integer *mode);
/*     DIMENSION W(MDW,N+1),PRGOPT(*),X(N), */
/*     WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2) */
/*     ABOVE, K=MAX(MA+MG,N). */
/*     ABSTRACT */
/*     THIS SUBPROGRAM SOLVES A LINEARLY CONSTRAINED LEAST SQUARES */
/*     PROBLEM WITH BOTH EQUALITY AND INEQUALITY CONSTRAINTS, AND, IF THE */
/*     USER REQUESTS, OBTAINS A COVARIANCE MATRIX OF THE SOLUTION */
/*     PARAMETERS. */

/*     SUPPOSE THERE ARE GIVEN MATRICES E, A AND G OF RESPECTIVE */
/*     DIMENSIONS ME BY N, MA BY N AND MG BY N, AND VECTORS F, B AND H OF */
/*     RESPECTIVE LENGTHS ME, MA AND MG.  THIS SUBROUTINE SOLVES THE */
/*     LINEARLY CONSTRAINED LEAST SQUARES PROBLEM */

/*                   EX = F, (E ME BY N) (EQUATIONS TO BE EXACTLY */
/*                                       SATISFIED) */
/*                   AX = B, (A MA BY N) (EQUATIONS TO BE */
/*                                       APPROXIMATELY SATISFIED, */
/*                                       LEAST SQUARES SENSE) */
/*                   GX.GE.H,(G MG BY N) (INEQUALITY CONSTRAINTS) */

/*     THE INEQUALITIES GX.GE.H MEAN THAT EVERY COMPONENT OF THE PRODUCT */
/*     GX MUST BE .GE. THE CORRESPONDING COMPONENT OF H. */

/*     IN CASE THE EQUALITY CONSTRAINTS CANNOT BE SATISFIED, A */
/*     GENERALIZED INVERSE SOLUTION RESIDUAL VECTOR LENGTH IS OBTAINED */
/*     FOR F-EX. THIS IS THE MINIMAL LENGTH POSSIBLE FOR F-EX. */
/*     ANY VALUES ME.GE.0, MA.GE.0, OR MG.GE.0 ARE PERMITTED.  THE */
/*     RANK OF THE MATRIX E IS ESTIMATED DURING THE COMPUTATION. WE CALL */
/*     THIS VALUE KRANKE. IT IS AN OUTPUT PARAMETER IN IP(1) DEFINED */
/*     BELOW. USING A GENERALIZED INVERSE SOLUTION OF EX=F, A REDUCED */
/*     LEAST SQUARES PROBLEM WITH INEQUALITY CONSTRAINTS IS OBTAINED. */
/*     THE TOLERANCES USED IN THESE TESTS FOR DETERMINING THE RANK */
/*     OF E AND THE RANK OF THE REDUCED LEAST SQUARES PROBLEM ARE */
/*     GIVEN IN SANDIA TECH. REPT. SAND 78-1290. THEY CAN BE */
/*     MODIFIED BY THE USER IF NEW VALUES ARE PROVIDED IN */
/*     THE OPTION LIST OF THE ARRAY PRGOPT(*). */
/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START EDITING AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SDOT/DDOT/, */
/*     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/ ABS/ DABS/,/AMAX1/DMAX1/, */
/*     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/SSWAP/DSWAP/,/E0/D0/, */
/*     /, DUMMY/,SNGL(DUMMY)/,/SRELPR/DRELPR/ */
/*     WRITTEN BY R. J. HANSON AND K. H. HASKELL.  FOR FURTHER MATH. */
/*     AND ALGORITHMIC DETAILS SEE SANDIA LABORATORIES TECH. REPTS. */
/*     SAND 77-0552, (1978), SAND 78-1290, (1979), AND */
/*     MATH. PROGRAMMING, VOL. 21, (1981), P.98-118. */

/*     THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST.. */
/*     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2) */
/*     WHERE K=MAX(MA+MG,N).  THIS ALLOWS FOR A SOLUTION OF A RANGE OF */
/*     PROBLEMS IN THE GIVEN WORKING SPACE.  THE DIMENSION OF WS(*) */
/*     GIVEN IS A NECESSARY OVERESTIMATE.  ONCE A PARTICULAR PROBLEM */
/*     HAS BEEN RUN, THE OUTPUT PARAMETER IP(3) GIVES THE ACTUAL */
/*     DIMENSION REQUIRED FOR THAT PROBLEM. */
/*     THE PARAMETERS FOR LSEI( ) ARE */
/*     INPUT.. */
/*     W(*,*),MDW,   THE ARRAY W(*,*) IS DOUBLY SUBSCRIPTED WITH */
/*     ME,MA,MG,N    FIRST DIMENSIONING PARAMETER EQUAL TO MDW. */
/*                   FOR THIS DISCUSSION LET US CALL M = ME+MA+MG.  THEN */
/*                   MDW MUST SATISFY MDW.GE.M.  THE CONDITION */
/*                   MDW.LT.M IS AN ERROR. */
/*                   THE ARRAY W(*,*) CONTAINS THE MATRICES AND VECTORS */
/*                                  (E  F) */
/*                                  (A  B) */
/*                                  (G  H) */
/*                   IN ROWS AND COLUMNS 1,...,M AND 1,...,N+1 */
/*                   RESPECTIVELY. */
/*                   THE INTEGERS ME, MA, AND MG ARE THE */
/*                   RESPECTIVE MATRIX ROW DIMENSIONS */
/*                   OF E, A AND G. EACH MATRIX HAS N COLUMNS. */
/*     PRGOPT(*)    THIS ARRAY IS THE OPTION VECTOR. */
/*                  IF THE USER IS SATISFIED WITH THE NOMINAL */
/*                  SUBPROGRAM FEATURES SET */
/*                  PRGOPT(1)=1 (OR PRGOPT(1)=1.0) */
/*                  OTHERWISE PRGOPT(*) IS A LINKED LIST CONSISTING OF */
/*                  GROUPS OF DATA OF THE FOLLOWING FORM */
/*                  LINK */
/*                  KEY */
/*                  DATA SET */

/*                  THE PARAMETERS LINK AND KEY ARE EACH ONE WORD. */
/*                  THE DATA SET CAN BE COMPRISED OF SEVERAL WORDS. */
/*                  THE NUMBER OF ITEMS DEPENDS ON THE VALUE OF KEY. */
/*                  THE VALUE OF LINK POINTS TO THE FIRST */
/*                  ENTRY OF THE NEXT GROUP OF DATA WITHIN */
/*                  PRGOPT(*).  THE EXCEPTION IS WHEN THERE ARE */
/*                  NO MORE OPTIONS TO CHANGE.  IN THAT */
/*                  CASE LINK=1 AND THE VALUES KEY AND DATA SET */
/*                  ARE NOT REFERENCED. THE GENERAL LAYOUT OF */
/*                  PRGOPT(*) IS AS FOLLOWS. */
/*               ...PRGOPT(1)=LINK1 (LINK TO FIRST ENTRY OF NEXT GROUP) */
/*               .  PRGOPT(2)=KEY1 (KEY TO THE OPTION CHANGE) */
/*               .  PRGOPT(3)=DATA VALUE (DATA VALUE FOR THIS CHANGE) */
/*               .       . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK1)=LINK2 (LINK TO THE FIRST ENTRY OF */
/*               .                       NEXT GROUP) */
/*               .  PRGOPT(LINK1+1)=KEY2 (KEY TO THE OPTION CHANGE) */
/*               .  PRGOPT(LINK1+2)=DATA VALUE */
/*               ...     . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK)=1 (NO MORE OPTIONS TO CHANGE) */
/*                  VALUES OF LINK THAT ARE NONPOSITIVE ARE ERRORS. */
/*                  A VALUE OF LINK.GT.NLINK=100000 IS ALSO AN ERROR. */
/*                  THIS HELPS PREVENT USING INVALID BUT POSITIVE */
/*                  VALUES OF LINK THAT WILL PROBABLY EXTEND */
/*                  BEYOND THE PROGRAM LIMITS OF PRGOPT(*). */
/*                  UNRECOGNIZED VALUES OF KEY ARE IGNORED.  THE */
/*                  ORDER OF THE OPTIONS IS ARBITRARY AND ANY NUMBER */
/*                  OF OPTIONS CAN BE CHANGED WITH THE FOLLOWING */
/*                  RESTRICTION.  TO PREVENT CYCLING IN THE */
/*                  PROCESSING OF THE OPTION ARRAY A COUNT OF THE */
/*                  NUMBER OF OPTIONS CHANGED IS MAINTAINED. */
/*                  WHENEVER THIS COUNT EXCEEDS NOPT=1000 AN ERROR */
/*                  MESSAGE IS PRINTED AND THE SUBPROGRAM RETURNS. */
/*                  OPTIONS.. */
/*                  KEY=1 */
/*                         COMPUTE IN W(*,*) THE N BY N */
/*                  COVARIANCE MATRIX OF THE SOLUTION VARIABLES */
/*                  AS AN OUTPUT PARAMETER.  NOMINALLY THE */
/*                  COVARIANCE MATRIX WILL NOT BE COMPUTED. */
/*                  (THIS REQUIRES NO USER INPUT.) */
/*                  THE DATA SET FOR THIS OPTION IS A SINGLE VALUE. */
/*                  IT MUST BE NONZERO WHEN THE COVARIANCE MATRIX */
/*                  IS DESIRED.  IF IT IS ZERO, THE COVARIANCE */
/*                  MATRIX IS NOT COMPUTED.  WHEN THE COVARIANCE MATRIX */
/*                  IS COMPUTED, THE FIRST DIMENSIONING PARAMETER */
/*                  OF THE ARRAY W(*,*) MUST SATISFY MDW.GE.MAX0(M,N). */

/*                  KEY=2 */
/*                         SCALE THE NONZERO COLUMNS OF THE */
/*                         ENTIRE DATA MATRIX. */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  TO HAVE LENGTH ONE.  THE DATA SET FOR THIS */
/*                  OPTION IS A SINGLE VALUE.  IT MUST BE */
/*                  NONZERO IF UNIT LENGTH COLUMN SCALING */
/*                  IS DESIRED. */

/*                  KEY=3 */
/*                         SCALE COLUMNS OF THE ENTIRE DATA MATRIX */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  WITH A USER-PROVIDED DIAGONAL MATRIX. */
/*                  THE DATA SET FOR THIS OPTION CONSISTS */
/*                  OF THE N DIAGONAL SCALING FACTORS, ONE FOR */
/*                  EACH MATRIX COLUMN. */

/*                  KEY=4 */
/*                         CHANGE THE RANK DETERMINATION TOLERANCE FOR */
/*                  THE EQUALITY CONSTRAINT EQUATIONS FROM */
/*                  THE NOMINAL VALUE OF SQRT(SRELPR).  THIS QUANTITY CAN */
/*                  BE NO SMALLER THAN SRELPR, THE ARITHMETIC- */
/*                  STORAGE PRECISION.  THE QUANTITY SRELPR IS THE */
/*                  LARGEST POSITIVE NUMBER SUCH THAT T=1.+SRELPR */
/*                  SATISFIES T.EQ.1.  THE QUANTITY USED */
/*                  HERE IS INTERNALLY RESTRICTED TO BE AT */
/*                  LEAST SRELPR.  THE DATA SET FOR THIS OPTION */
/*                  IS THE NEW TOLERANCE. */

/*                  KEY=5 */
/*                         CHANGE THE RANK DETERMINATION TOLERANCE FOR */
/*                  THE REDUCED LEAST SQUARES EQUATIONS FROM */
/*                  THE NOMINAL VALUE OF SQRT(SRELPR).  THIS QUANTITY CAN */
/*                  BE NO SMALLER THAN SRELPR, THE ARITHMETIC- */
/*                  STORAGE PRECISION.  THE QUANTITY USED */
/*                  HERE IS INTERNALLY RESTRICTED TO BE AT */
/*                  LEAST SRELPR.  THE DATA SET FOR THIS OPTION */
/*                  IS THE NEW TOLERANCE. */

/*                  FOR EXAMPLE, SUPPOSE WE WANT TO CHANGE */
/*                  THE TOLERANCE FOR THE REDUCED LEAST SQUARES */
/*                  PROBLEM, COMPUTE THE COVARIANCE MATRIX OF */
/*                  THE SOLUTION PARAMETERS, AND PROVIDE */
/*                  COLUMN SCALING FOR THE DATA MATRIX.  FOR */
/*                  THESE OPTIONS THE DIMENSION OF PRGOPT(*) */
/*                  MUST BE AT LEAST N+9.  THE FORTRAN STATEMENTS */
/*                  DEFINING THESE OPTIONS WOULD BE AS FOLLOWS. */

/*                  PRGOPT(1)=4 (LINK TO ENTRY 4 IN PRGOPT(*)) */
/*                  PRGOPT(2)=1 (COVARIANCE MATRIX KEY) */
/*                  PRGOPT(3)=1 (COVARIANCE MATRIX WANTED) */

/*                  PRGOPT(4)=7 (LINK TO ENTRY 7 IN PRGOPT(*)) */
/*                  PRGOPT(5)=5 (LEAST SQUARES EQUAS. TOLERANCE KEY) */
/*                  PRGOPT(6)=... (NEW VALUE OF THE TOLERANCE) */

/*                  PRGOPT(7)=N+9 (LINK TO ENTRY N+9 IN PRGOPT(*)) */
/*                  PRGOPT(8)=3 (USER-PROVIDED COLUMN SCALING KEY) */

/*                  CALL SCOPY(N,D,1,PRGOPT(9),1)  (COPY THE N */
/*                    SCALING FACTORS FROM THE USER ARRAY D(*) */
/*                    TO PRGOPT(9)-PRGOPT(N+8)) */

/*                  PRGOPT(N+9)=1 (NO MORE OPTIONS TO CHANGE) */

/*                  THE CONTENTS OF PRGOPT(*) ARE NOT MODIFIED */
/*                  BY THE SUBPROGRAM. */
/*                  THE OPTIONS FOR WNNLS( ) CAN ALSO BE INCLUDED */
/*                  IN THIS ARRAY.  THE VALUES OF KEY RECOGNIZED */
/*                  BY WNNLS( ) ARE 6, 7 AND 8.  THEIR FUNCTIONS */
/*                  ARE DOCUMENTED IN THE USAGE INSTRUCTIONS FOR */
/*                  SUBROUTINE WNNLS( ).  NORMALLY THESE OPTIONS */
/*                  DO NOT NEED TO BE MODIFIED WHEN USING LSEI( ). */

/*     IP(1),       THE AMOUNTS OF WORKING STORAGE ACTUALLY */
/*     IP(2)        ALLOCATED FOR THE WORKING ARRAYS WS(*) AND */
/*                  IP(*), RESPECTIVELY.  THESE QUANTITIES ARE */
/*                  COMPARED WITH THE ACTUAL AMOUNTS OF STORAGE */
/*                  NEEDED BY LSEI( ).  INSUFFICIENT STORAGE */
/*                  ALLOCATED FOR EITHER WS(*) OR IP(*) IS AN */
/*                  ERROR.  THIS FEATURE WAS INCLUDED IN LSEI( ) */
/*                  BECAUSE MISCALCULATING THE STORAGE FORMULAS */
/*                  FOR WS(*) AND IP(*) MIGHT VERY WELL LEAD TO */
/*                  SUBTLE AND HARD-TO-FIND EXECUTION ERRORS. */

/*                  THE LENGTH OF WS(*) MUST BE AT LEAST */

/*                  LW = 2*(ME+N)+K+(MG+2)*(N+7) */


/*                  WHERE K = MAX(MA+MG,N) */
/*                  THIS TEST WILL NOT BE MADE IF IP(1).LE.0. */

/*                  THE LENGTH OF IP(*) MUST BE AT LEAST */

/*                  LIP = MG+2*N+2 */
/*                  THIS TEST WILL NOT BE MADE IF IP(2).LE.0. */

/*     OUTPUT.. */

/*     X(*),RNORME,  THE ARRAY X(*) CONTAINS THE SOLUTION PARAMETERS */
/*     RNORML        IF THE INTEGER OUTPUT FLAG MODE = 0 OR 1. */
/*                   THE DEFINITION OF MODE IS GIVEN DIRECTLY BELOW. */
/*                   WHEN MODE = 0 OR 1, RNORME AND RNORML */
/*                   RESPECTIVELY CONTAIN THE RESIDUAL VECTOR */
/*                   EUCLIDEAN LENGTHS OF F - EX AND B - AX.  WHEN */
/*                   MODE=1 THE EQUALITY CONSTRAINT EQUATIONS EX=F */
/*                   ARE CONTRADICTORY, SO RNORME.NE.0. THE RESIDUAL */
/*                   VECTOR F-EX HAS MINIMAL EUCLIDEAN LENGTH. FOR */
/*                   MODE.GE.2, NONE OF THESE PARAMETERS ARE */
/*                   DEFINED. */

/*     MODE          INTEGER FLAG THAT INDICATES THE SUBPROGRAM */
/*                   STATUS AFTER COMPLETION.  IF MODE.GE.2, NO */
/*                   SOLUTION HAS BEEN COMPUTED. */

/*                   MODE = */

/*                   0  BOTH EQUALITY AND INEQUALITY CONSTRAINTS */
/*                      ARE COMPATIBLE AND HAVE BEEN SATISFIED. */

/*                   1  EQUALITY CONSTRAINTS ARE CONTRADICTORY. */
/*                      A GENERALIZED INVERSE SOLUTION OF EX=F WAS USED */
/*                      TO MINIMIZE THE RESIDUAL VECTOR LENGTH F-EX. */
/*                      IN THIS SENSE, THE SOLUTION IS STILL MEANINGFUL. */

/*                   2  INEQUALITY CONSTRAINTS ARE CONTRADICTORY. */

/*                   3  BOTH EQUALITY AND INEQUALITY CONSTRAINTS */
/*                      ARE CONTRADICTORY. */

/*                   THE FOLLOWING INTERPRETATION OF */
/*                   MODE=1,2 OR 3 MUST BE MADE.  THE */
/*                   SETS CONSISTING OF ALL SOLUTIONS */
/*                   OF THE EQUALITY CONSTRAINTS EX=F */
/*                   AND ALL VECTORS SATISFYING GX.GE.H */
/*                   HAVE NO POINTS ON COMMON.  (IN */
/*                   PARTICULAR THIS DOES NOT SAY THAT */
/*                   EACH INDIVIDUAL SET HAS NO POINTS */
/*                   AT ALL, ALTHOUGH THIS COULD BE THE */
/*                   CASE.) */

/*                   4  USAGE ERROR OCCURRED.  THE VALUE */
/*                      OF MDW IS .LT. ME+MA+MG, MDW IS */
/*                      .LT. N AND A COVARIANCE MATRIX IS */
/*                      REQUESTED, THE OPTION VECTOR */
/*                      PRGOPT(*) IS NOT PROPERLY DEFINED, */
/*                      OR THE LENGTHS OF THE WORKING ARRAYS */
/*                      WS(*) AND IP(*), WHEN SPECIFIED IN */
/*                      IP(1) AND IP(2) RESPECTIVELY, ARE NOT */
/*                      LONG ENOUGH. */

/*     W(*,*)        THE ARRAY W(*,*) CONTAINS THE N BY N SYMMETRIC */
/*                   COVARIANCE MATRIX OF THE SOLUTION PARAMETERS, */
/*                   PROVIDED THIS WAS REQUESTED ON INPUT WITH */
/*                   THE OPTION VECTOR PRGOPT(*) AND THE OUTPUT */
/*                   FLAG IS RETURNED WITH MODE = 0 OR 1. */

/*     IP(*)         THE INTEGER WORKING ARRAY HAS THREE ENTRIES */
/*                   THAT PROVIDE RANK AND WORKING ARRAY LENGTH */
/*                   INFORMATION AFTER COMPLETION. */

/*                      IP(1) = RANK OF EQUALITY CONSTRAINT */
/*                              MATRIX.  DEFINE THIS QUANTITY */
/*                              AS KRANKE. */

/*                      IP(2) = RANK OF REDUCED LEAST SQUARES */
/*                              PROBLEM. */

/*                      IP(3) = THE AMOUNT OF STORAGE IN THE */
/*                              WORKING ARRAY WS(*) THAT WAS */
/*                              ACTUALLY USED BY THE SUBPROGRAM. */
/*                              THE FORMULA GIVEN ABOVE FOR THE LENGTH */
/*                              OF WS(*) IS A NECESSARY OVERESTIMATE. */
/*     USER DESIGNATED */
/*     WORKING ARRAYS.. */

/*     WS(*),IP(*)              THESE ARE RESP. TYPE FLOATING POINT */
/*                              AND TYPE INTEGER WORKING ARRAYS. */
/*                              THEIR REQUIRED MINIMAL LENGTHS ARE */
/*                              GIVEN ABOVE. */


/*     SUBROUTINES CALLED */

/*     LSI           PART OF THIS PACKAGE.  SOLVES A */
/*                   CONSTRAINED LEAST SQUARES PROBLEM WITH */
/*                   INEQUALITY CONSTRAINTS. */

/* ++ */
/*     SDOT,SSCAL,   SUBROUTINES FROM THE BLAS PACKAGE. */
/*     SAXPY,SASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308. */
/*     SCOPY,SNRM2, */
/*     SSWAP */

/*     H12           SUBROUTINE TO CONSTRUCT AND APPLY A */
/*                   HOUSEHOLDER TRANSFORMATION. */

/*     XERROR        FROM SLATEC ERROR PROCESSING PACKAGE. */
/*                   THIS IS DOCUMENTED IN SANDIA TECH. REPT., */
/*                   SAND78-1189. */



/*     SUBROUTINE LSI(W,MDW,MA,MG,N,PRGOPT,X,RNORM,MODE,WS,IP) */
public:
    /* Subroutine */ int lsi_(double *w, integer *mdw, integer *ma, integer *mg, 
    integer *n, double *prgopt, double *x, double *rnorm, integer *mode, double *ws, integer *ip);

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START EDITING AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SDOT/DDOT/, */
/*     / SQRT/ DSQRT/,/AMAX1/DMAX1/,/SSWAP/DSWAP/, */
/*     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/E0/D0/,/SRELPR/DRELPR/ */

/*     THIS IS A COMPANION SUBPROGRAM TO LSEI( ). */
/*     THE DOCUMENTATION FOR LSEI( ) HAS MORE COMPLETE */
/*     USAGE INSTRUCTIONS. */
/*     WRITTEN BY R. J. HANSON, SLA. */

/*     SOLVE.. */
/*              AX = B,  A  MA BY N  (LEAST SQUARES EQUATIONS) */
/*     SUBJECT TO.. */

/*              GX.GE.H, G  MG BY N  (INEQUALITY CONSTRAINTS) */

/*     INPUT.. */

/*      W(*,*) CONTAINS  (A B) IN ROWS 1,...,MA+MG, COLS 1,...,N+1. */
/*                       (G H) */

/*     MDW,MA,MG,N */
/*              CONTAIN (RESP) VAR. DIMENSION OF W(*,*), */
/*              AND MATRIX DIMENSIONS. */

/*     PRGOPT(*), */
/*              PROGRAM OPTION VECTOR. */

/*     OUTPUT.. */

/*      X(*),RNORM */

/*              SOLUTION VECTOR(UNLESS MODE=2), LENGTH OF AX-B. */

/*      MODE */
/*              =0   INEQUALITY CONSTRAINTS ARE COMPATIBLE. */
/*              =2   INEQUALITY CONSTRAINTS CONTRADICTORY. */

/*      WS(*), */
/*              WORKING STORAGE OF DIMENSION K+N+(MG+2)*(N+7), */
/*              WHERE K=MAX(MA+MG,N). */
/*      IP(MG+2*N+1) */
/*              INTEGER WORKING STORAGE */
/*      REVISED OCT. 1, 1981. */

/*     SUBROUTINES CALLED */

/*     LPDP          THIS SUBPROGRAM MINIMIZES A SUM OF SQUARES */
/*                   OF UNKNOWNS SUBJECT TO LINEAR INEQUALITY */
/*                   CONSTRAINTS.  PART OF THIS PACKAGE. */

/* ++ */
/*     SDOT,SSCAL    SUBROUTINES FROM THE BLAS PACKAGE. */
/*     SAXPY,SASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308. */
/*     SCOPY,SSWAP */

/*     HFTI          SOLVES AN UNCONSTRAINED LINEAR LEAST SQUARES */
/*                   PROBLEM.  PART OF THIS PACKAGE. */

/*     H12           SUBROUTINE TO CONSTRUCT AND APPLY A HOUSEHOLDER */
/*                   TRANSFORMATION. */


/*     SDOT,         SUBROUTINES FROM THE BLAS PACKAGE. */
/*     SSCAL,SNRM2,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308. */
/*     SCOPY */ 

    /* Subroutine */ int lpdp_(double *a, integer *mda, integer *m, integer *n1, 
    integer *n2, double *prgopt, double *x, double *wnorm, integer *mode, double *ws, integer *is);
    
/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START EDITING AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SNRM2/DNRM2/,/SDOT/DDOT/, */
/*     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/ABS(/DABS(/, ABS/, DABS/,/E0/D0/ */

/*     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1), */
/*     WHERE N=N1+N2.  THIS IS A SLIGHT OVERESTIMATE FOR WS(*). */

/*     WRITTEN BY R. J. HANSON AND K. H. HASKELL, SANDIA LABS */
/*     REVISED OCT. 1, 1981. */

/*     DETERMINE AN N1-VECTOR W, AND */
/*               AN N2-VECTOR Z */
/*     WHICH MINIMIZES THE EUCLIDEAN LENGTH OF W */
/*     SUBJECT TO G*W+H*Z .GE. Y. */
/*     THIS IS THE LEAST PROJECTED DISTANCE PROBLEM, LPDP. */
/*     THE MATRICES G AND H ARE OF RESPECTIVE */
/*     DIMENSIONS M BY N1 AND M BY N2. */

/*     CALLED BY SUBPROGRAM LSI( ). */

/*     THE MATRIX */
/*                (G H Y) */

/*     OCCUPIES ROWS 1,...,M AND COLS 1,...,N1+N2+1 OF A(*,*). */

/*     THE SOLUTION (W) IS RETURNED IN X(*). */
/*                  (Z) */

/*     THE VALUE OF MODE INDICATES THE STATUS OF */
/*     THE COMPUTATION AFTER RETURNING TO THE USER. */

/*          MODE=1  THE SOLUTION WAS SUCCESSFULLY OBTAINED. */

/*          MODE=2  THE INEQUALITIES ARE INCONSISTENT. */

/*     SUBROUTINES CALLED */

/*     WNNLS         SOLVES A NONNEGATIVELY CONSTRAINED LINEAR LEAST */
/*                   SQUARES PROBLEM WITH LINEAR EQUALITY CONSTRAINTS. */
/*                   PART OF THIS PACKAGE. */
private:
    /* Subroutine */ int wnnls_(double *w, integer *mdw, integer *me, integer *ma, integer *n, 
    integer *l, double *prgopt, double *x, double *rnorm, integer *mode, integer *iwork, double *work);
public:
    void printErrors();

    /* Subroutine */ int wnnls_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *n, integer *l, double *prgopt, double *x, double *rnorm, integer *mode);
/*     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),IWORK(M+N),WORK(M+5*N) */
/*     ABSTRACT */
/*     THIS SUBPROGRAM SOLVES A LINEARLY CONSTRAINED LEAST SQUARES */
/*     PROBLEM.  SUPPOSE THERE ARE GIVEN MATRICES E AND A OF */
/*     RESPECTIVE DIMENSIONS ME BY N AND MA BY N, AND VECTORS F */
/*     AND B OF RESPECTIVE LENGTHS ME AND MA.  THIS SUBROUTINE */
/*     SOLVES THE PROBLEM */
/*               EX = F, (EQUATIONS TO BE EXACTLY SATISFIED) */
/*               AX = B, (EQUATIONS TO BE APPROXIMATELY SATISFIED, */
/*                        IN THE LEAST SQUARES SENSE) */
/*               SUBJECT TO COMPONENTS L+1,...,N NONNEGATIVE */
/*     ANY VALUES ME.GE.0, MA.GE.0 AND 0.LE. L .LE.N ARE PERMITTED. */
/*     THE PROBLEM IS REPOSED AS PROBLEM WNNLS */

/*               (WT*E)X = (WT*F) */
/*               (   A)    (   B), (LEAST SQUARES) */
/*               SUBJECT TO COMPONENTS L+1,...,N NONNEGATIVE. */

/*     THE SUBPROGRAM CHOOSES THE HEAVY WEIGHT (OR PENALTY PARAMETER) WT. */

/*     THE PARAMETERS FOR WNNLS ARE */

/*     INPUT.. */

/*     W(*,*),MDW,  THE ARRAY W(*,*) IS DOUBLE SUBSCRIPTED WITH FIRST */
/*     ME,MA,N,L    DIMENSIONING PARAMETER EQUAL TO MDW.  FOR THIS */
/*                  DISCUSSION LET US CALL M = ME + MA.  THEN MDW */
/*                  MUST SATISFY MDW.GE.M.  THE CONDITION MDW.LT.M */
/*                  IS AN ERROR. */

/*                  THE ARRAY W(*,*) CONTAINS THE MATRICES AND VECTORS */

/*                       (E  F) */
/*                       (A  B) */

/*                  IN ROWS AND COLUMNS 1,...,M AND 1,...,N+1 */
/*                  RESPECTIVELY.  COLUMNS 1,...,L CORRESPOND TO */
/*                  UNCONSTRAINED VARIABLES X(1),...,X(L).  THE */
/*                  REMAINING VARIABLES ARE CONSTRAINED TO BE */
/*                  NONNEGATIVE.  THE CONDITION L.LT.0 .OR. L.GT.N IS */
/*                  AN ERROR. */

/*     PRGOPT(*)    THIS ARRAY IS THE OPTION VECTOR. */
/*                  IF THE USER IS SATISFIED WITH THE NOMINAL */
/*                  SUBPROGRAM FEATURES SET */

/*                  PRGOPT(1)=1 (OR PRGOPT(1)=1.0) */

/*                  OTHERWISE PRGOPT(*) IS A LINKED LIST CONSISTING OF */
/*                  GROUPS OF DATA OF THE FOLLOWING FORM */

/*                  LINK */
/*                  KEY */
/*                  DATA SET */

/*                  THE PARAMETERS LINK AND KEY ARE EACH ONE WORD. */
/*                  THE DATA SET CAN BE COMPRISED OF SEVERAL WORDS. */
/*                  THE NUMBER OF ITEMS DEPENDS ON THE VALUE OF KEY. */
/*                  THE VALUE OF LINK POINTS TO THE FIRST */
/*                  ENTRY OF THE NEXT GROUP OF DATA WITHIN */
/*                  PRGOPT(*).  THE EXCEPTION IS WHEN THERE ARE */
/*                  NO MORE OPTIONS TO CHANGE.  IN THAT */
/*                  CASE LINK=1 AND THE VALUES KEY AND DATA SET */
/*                  ARE NOT REFERENCED. THE GENERAL LAYOUT OF */
/*                  PRGOPT(*) IS AS FOLLOWS. */

/*               ...PRGOPT(1)=LINK1 (LINK TO FIRST ENTRY OF NEXT GROUP) */
/*               .  PRGOPT(2)=KEY1 (KEY TO THE OPTION CHANGE) */
/*               .  PRGOPT(3)=DATA VALUE (DATA VALUE FOR THIS CHANGE) */
/*               .       . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK1)=LINK2 (LINK TO THE FIRST ENTRY OF */
/*               .                       NEXT GROUP) */
/*               .  PRGOPT(LINK1+1)=KEY2 (KEY TO THE OPTION CHANGE) */
/*               .  PRGOPT(LINK1+2)=DATA VALUE */
/*               ...     . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK)=1 (NO MORE OPTIONS TO CHANGE) */

/*                  VALUES OF LINK THAT ARE NONPOSITIVE ARE ERRORS. */
/*                  A VALUE OF LINK.GT.NLINK=100000 IS ALSO AN ERROR. */
/*                  THIS HELPS PREVENT USING INVALID BUT POSITIVE */
/*                  VALUES OF LINK THAT WILL PROBABLY EXTEND */
/*                  BEYOND THE PROGRAM LIMITS OF PRGOPT(*). */
/*                  UNRECOGNIZED VALUES OF KEY ARE IGNORED.  THE */
/*                  ORDER OF THE OPTIONS IS ARBITRARY AND ANY NUMBER */
/*                  OF OPTIONS CAN BE CHANGED WITH THE FOLLOWING */
/*                  RESTRICTION.  TO PREVENT CYCLING IN THE */
/*                  PROCESSING OF THE OPTION ARRAY A COUNT OF THE */
/*                  NUMBER OF OPTIONS CHANGED IS MAINTAINED. */
/*                  WHENEVER THIS COUNT EXCEEDS NOPT=1000 AN ERROR */
/*                  MESSAGE IS PRINTED AND THE SUBPROGRAM RETURNS. */

/*                  OPTIONS.. */

/*                  KEY=6 */
/*                         SCALE THE NONZERO COLUMNS OF THE */
/*                  ENTIRE DATA MATRIX */
/*                  (E) */
/*                  (A) */
/*                  TO HAVE LENGTH ONE.  THE DATA SET FOR */
/*                  THIS OPTION IS A SINGLE VALUE.  IT MUST */
/*                  BE NONZERO IF UNIT LENGTH COLUMN SCALING IS */
/*                  DESIRED. */

/*                  KEY=7 */
/*                         SCALE COLUMNS OF THE ENTIRE DATA MATRIX */
/*                  (E) */
/*                  (A) */
/*                  WITH A USER-PROVIDED DIAGONAL MATRIX. */
/*                  THE DATA SET FOR THIS OPTION CONSISTS */
/*                  OF THE N DIAGONAL SCALING FACTORS, ONE FOR */
/*                  EACH MATRIX COLUMN. */

/*                  KEY=8 */
/*                         CHANGE THE RANK DETERMINATION TOLERANCE FROM */
/*                  THE NOMINAL VALUE OF SQRT(EPS).  THIS QUANTITY CAN */
/*                  BE NO SMALLER THAN EPS, THE ARITHMETIC- */
/*                  STORAGE PRECISION.  THE QUANTITY USED */
/*                  HERE IS INTERNALLY RESTRICTED TO BE AT */
/*                  LEAST EPS.  THE DATA SET FOR THIS OPTION */
/*                  IS THE NEW TOLERANCE. */

/*                  KEY=9 */
/*                         CHANGE THE BLOW-UP PARAMETER FROM THE */
/*                  NOMINAL VALUE OF SQRT(EPS).  THE RECIPROCAL OF */
/*                  THIS PARAMETER IS USED IN REJECTING SOLUTION */
/*                  COMPONENTS AS TOO LARGE WHEN A VARIABLE IS */
/*                  FIRST BROUGHT INTO THE ACTIVE SET.  TOO LARGE */
/*                  MEANS THAT THE PROPOSED COMPONENT TIMES THE */
/*                  RECIPROCAL OF THE PARAMETERIS NOT LESS THAN */
/*                  THE RATIO OF THE NORMS OF THE RIGHT-SIDE */
/*                  VECTOR AND THE DATA MATRIX. */
/*                  THIS PARAMETER CAN BE NO SMALLER THAN EPS, */
/*                  THE ARITHMETIC-STORAGE PRECISION. */

/*                  FOR EXAMPLE, SUPPOSE WE WANT TO PROVIDE */
/*                  A DIAGONAL MATRIX TO SCALE THE PROBLEM */
/*                  MATRIX AND CHANGE THE TOLERANCE USED FOR */
/*                  DETERMINING LINEAR DEPENDENCE OF DROPPED COL */
/*                  VECTORS.  FOR THESE OPTIONS THE DIMENSIONS OF */
/*                  PRGOPT(*) MUST BE AT LEAST N+6.  THE FORTRAN */
/*                  STATEMENTS DEFINING THESE OPTIONS WOULD */
/*                  BE AS FOLLOWS. */

/*                  PRGOPT(1)=N+3 (LINK TO ENTRY N+3 IN PRGOPT(*)) */
/*                  PRGOPT(2)=7 (USER-PROVIDED SCALING KEY) */

/*                  CALL SCOPY(N,D,1,PRGOPT(3),1) (COPY THE N */
/*                  SCALING FACTORS FROM A USER ARRAY CALLED D(*) */
/*                  INTO PRGOPT(3)-PRGOPT(N+2)) */

/*                  PRGOPT(N+3)=N+6 (LINK TO ENTRY N+6 OF PRGOPT(*)) */
/*                  PRGOPT(N+4)=8 (LINEAR DEPENDENCE TOLERANCE KEY) */
/*                  PRGOPT(N+5)=... (NEW VALUE OF THE TOLERANCE) */

/*                  PRGOPT(N+6)=1 (NO MORE OPTIONS TO CHANGE) */

/*     IWORK(1),    THE AMOUNTS OF WORKING STORAGE ACTUALLY ALLOCATED */
/*     IWORK(2)     FOR THE WORKING ARRAYS WORK(*) AND IWORK(*), */
/*                  RESPECTIVELY.  THESE QUANTITIES ARE COMPARED WITH */
/*                  THE ACTUAL AMOUNTS OF STORAGE NEEDED FOR WNNLS( ). */
/*                  INSUFFICIENT STORAGE ALLOCATED FOR EITHER WORK(*) */
/*                  OR IWORK(*) IS CONSIDERED AN ERROR.  THIS FEATURE */
/*                  WAS INCLUDED IN WNNLS( ) BECAUSE MISCALCULATING */
/*                  THE STORAGE FORMULAS FOR WORK(*) AND IWORK(*) */
/*                  MIGHT VERY WELL LEAD TO SUBTLE AND HARD-TO-FIND */
/*                  EXECUTION ERRORS. */

/*                  THE LENGTH OF WORK(*) MUST BE AT LEAST */

/*                  LW = ME+MA+5*N */
/*                  THIS TEST WILL NOT BE MADE IF IWORK(1).LE.0. */

/*                  THE LENGTH OF IWORK(*) MUST BE AT LEAST */

/*                  LIW = ME+MA+N */
/*                  THIS TEST WILL NOT BE MADE IF IWORK(2).LE.0. */

/*     OUTPUT.. */

/*     X(*)         AN ARRAY DIMENSIONED AT LEAST N, WHICH WILL */
/*                  CONTAIN THE N COMPONENTS OF THE SOLUTION VECTOR */
/*                  ON OUTPUT. */

/*     RNORM        THE RESIDUAL NORM OF THE SOLUTION.  THE VALUE OF */
/*                  RNORM CONTAINS THE RESIDUAL VECTOR LENGTH OF THE */
/*                  EQUALITY CONSTRAINTS AND LEAST SQUARES EQUATIONS. */

/*     MODE         THE VALUE OF MODE INDICATES THE SUCCESS OR FAILURE */
/*                  OF THE SUBPROGRAM. */

/*                  MODE = 0  SUBPROGRAM COMPLETED SUCCESSFULLY. */

/*                       = 1  MAX. NUMBER OF ITERATIONS (EQUAL TO */
/*                            3*(N-L)) EXCEEDED. NEARLY ALL PROBLEMS */
/*                            SHOULD COMPLETE IN FEWER THAN THIS */
/*                            NUMBER OF ITERATIONS. AN APPROXIMATE */
/*                            SOLUTION AND ITS CORRESPONDING RESIDUAL */
/*                            VECTOR LENGTH ARE IN X(*) AND RNORM. */

/*                       = 2  USAGE ERROR OCCURRED.  THE OFFENDING */
/*                            CONDITION IS NOTED WITH THE ERROR */
/*                            PROCESSING SUBPROGRAM, XERROR( ). */

/*     USER-DESIGNATED */
/*     WORKING ARRAYS.. */

/*     WORK(*)      A WORKING ARRAY OF LENGTH AT LEAST */
/*                  M + 5*N. */

/*     IWORK(*)     AN INTEGER-VALUED WORKING ARRAY OF LENGTH AT LEAST */
/*                  M+N. */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/, DUMMY/,SNGL(DUMMY)/ */

/*     WRITTEN BY KAREN H. HASKELL, SANDIA LABORATORIES, */
/*     AND R.J. HANSON, SANDIA LABORATORIES. */
/*     REVISED FEB.25, 1982. */

/*     SUBROUTINES CALLED BY WNNLS( ) */

/* ++ */
/*     WNLSM         COMPANION SUBROUTINE TO WNNLS( ), WHERE */
/*                   MOST OF THE COMPUTATION TAKES PLACE. */

/*     XERROR,XERRWV FROM SLATEC ERROR PROCESSING PACKAGE. */
/*                   THIS IS DOCUMENTED IN SANDIA TECH. REPT., */
/*                   SAND78-1189. */

/*     REFERENCES */

/*     1. SOLVING LEAST SQUARES PROBLEMS, BY C.L. LAWSON */
/*        AND R.J. HANSON.  PRENTICE-HALL, INC. (1974). */

/*     2. BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE, BY */
/*        C.L. LAWSON, R.J. HANSON, D.R. KINCAID, AND F.T. KROGH. */
/*        TOMS, V. 5, NO. 3, P. 308.  ALSO AVAILABLE AS */
/*        SANDIA TECHNICAL REPORT NO. SAND77-0898. */

/*     3. AN ALGORITHM FOR LINEAR LEAST SQUARES WITH EQUALITY */
/*        AND NONNEGATIVITY CONSTRAINTS, BY K.H. HASKELL AND */
/*        R.J. HANSON.  AVAILABLE AS SANDIA TECHNICAL REPORT NO. */
/*        SAND77-0552, AND MATH. PROGRAMMING, VOL. 21, (1981), P. 98-118. */

/*     4. SLATEC COMMON MATH. LIBRARY ERROR HANDLING */
/*        PACKAGE.  BY R. E. JONES.  AVAILABLE AS SANDIA */
/*        TECHNICAL REPORT SAND78-1189. */


/*     SUBROUTINE WNLSM (W,MDW,MME,MA,N,L,PRGOPT,X,RNORM,MODE, */
/*    1                  IPIVOT,ITYPE,WD,H,SCALE,Z,TEMP,D) */
    /* Subroutine */ int wnlsm_(double *w, integer *mdw, integer *mme, integer *ma, 
    integer *n, integer *l, double *prgopt, double *x, double *rnorm, integer *
    mode, integer *ipivot, integer *itype, double *wd, double *h__, double *
    scale, double *z__, double *temp, double *d__);
    
/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START CHANGES AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SROTMG/DROTMG/, */
/*     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/SROTM/DROTM/,/AMAX1/DMAX1/, */
/*     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/E0/D0/,/SSWAP/DSWAP/, */
/*     /ISAMAX/IDAMAX/,/SRELPR/DRELPR/,/.E-/.D-/                   REMK */

/*     THIS IS A COMPANION SUBPROGRAM TO WNNLS( ). */
/*     THE DOCUMENTATION FOR WNNLS( ) HAS MORE COMPLETE */
/*     USAGE INSTRUCTIONS. */

/*     WRITTEN BY KAREN H. HASKELL, SANDIA LABORATORIES, */
/*     WITH THE HELP OF R.J. HANSON, SANDIA LABORATORIES, */
/*     DECEMBER 1976 - JANUARY 1978. */
/*     REVISED MAR. 4, 1982. */

/*     IN ADDITION TO THE PARAMETERS DISCUSSED IN THE PROLOGUE TO */
/*     SUBROUTINE WNNLS, THE FOLLOWING WORK ARRAYS ARE USED IN */
/*     SUBROUTINE WNLSM  (THEY ARE PASSED THROUGH THE CALLING */
/*     SEQUENCE FROM WNNLS FOR PURPOSES OF VARIABLE DIMENSIONING). */
/*     THEIR CONTENTS WILL IN GENERAL BE OF NO INTEREST TO THE USER. */

/*         IPIVOT(*) */
/*            AN ARRAY OF LENGTH N.  UPON COMPLETION IT CONTAINS THE */
/*         PIVOTING INFORMATION FOR THE COLS OF W(*,*). */

/*         ITYPE(*) */
/*            AN ARRAY OF LENGTH M WHICH IS USED TO KEEP TRACK */
/*         OF THE CLASSIFICATION OF THE EQUATIONS.  ITYPE(I)=0 */
/*         DENOTES EQUATION I AS AN EQUALITY CONSTRAINT. */
/*         ITYPE(I)=1 DENOTES EQUATION I AS A LEAST SQUARES */
/*         EQUATION. */

/*         WD(*) */
/*            AN ARRAY OF LENGTH N.  UPON COMPLETION IT CONTAINS THE */
/*         DUAL SOLUTION VECTOR. */

/*         H(*) */
/*            AN ARRAY OF LENGTH N.  UPON COMPLETION IT CONTAINS THE */
/*         PIVOT SCALARS OF THE HOUSEHOLDER TRANSFORMATIONS PERFORMED */
/*         IN THE CASE KRANK.LT.L. */

/*         SCALE(*) */
/*            AN ARRAY OF LENGTH M WHICH IS USED BY THE SUBROUTINE */
/*         TO STORE THE DIAGONAL MATRIX OF WEIGHTS. */
/*         THESE ARE USED TO APPLY THE MODIFIED GIVENS */
/*         TRANSFORMATIONS. */

/*         Z(*),TEMP(*) */
/*            WORKING ARRAYS OF LENGTH N. */

/*         D(*) */
/*            AN ARRAY OF LENGTH N THAT CONTAINS THE */
/*         COLUMN SCALING FOR THE MATRIX (E). */
/*                                       (A) */

    /* Subroutine */ int wnlit_(double *, integer *, integer *, integer *
        , integer *, integer *, integer *, double *, double *, double *, 
        integer *, double *, logical *);

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (BEGIN CHANGES AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/,/SROTM/DROTM/, */
/*     /SSCAL/DSCAL/,/SQRT/DSQRT/,                                  REMK */
/*     /SSWAP/DSWAP/,/AMAX1/DMAX1/,/ISAMAX/IDAMAX/,/.E-/.D-/,/E0/D0/ */

/*     THIS IS A COMPANION SUBPROGRAM TO WNNLS( ). */
/*     THE DOCUMENTATION FOR WNNLS( ) HAS MORE COMPLETE */
/*     USAGE INSTRUCTIONS. */

/*     NOTE  THE M BY (N+1) MATRIX W( , ) CONTAINS THE RT. HAND SIDE */
/*           B AS THE (N+1)ST COL. */

/*     TRIANGULARIZE L1 BY L1 SUBSYSTEM, WHERE L1=MIN(M,L), WITH */
/*     COL INTERCHANGES. */
/*     REVISED MARCH 4, 1982 */

    /* Subroutine */ int hfti_(double *a, integer *mda, integer *m, integer *n, 
    double *b, integer *mdb, integer *nb, double *tau, integer *krank, double *
    rnorm, double *h__, double *g, integer *ip);
    
    /*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (BEGIN CHANGES AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/ SQRT/ DSQRT/, */
/*     /, ABS/, DABS/,/ABS(/DABS(/,/E0/D0/ */

/*     DIMENSION A(MDA,N),(B(MDB,NB) OR B(M)),RNORM(NB),H(N),G(N),IP(N) */

/*     WRITTEN BY C. L. LAWSON AND R. J. HANSON.  FROM THE BOOK SOLVING */
/*     LEAST SQUARES PROBLEMS, PRENTICE-HALL, INC. (1974). FOR FURTHER */
/*     ALGORITHMIC DETAILS SEE ALGORITHM HFTI IN CHAPTER 14. */

/*     ABSTRACT */

/*     THIS SUBROUTINE SOLVES A LINEAR LEAST SQUARES PROBLEM OR A SET OF */
/*     LINEAR LEAST SQUARES PROBLEMS HAVING THE SAME MATRIX BUT DIFFERENT */
/*     RIGHT-SIDE VECTORS.  THE PROBLEM DATA CONSISTS OF AN M BY N MATRIX */
/*     A, AN M BY NB MATRIX B, AND AN ABSOLUTE TOLERANCE PARAMETER TAU */
/*     WHOSE USAGE IS DESCRIBED BELOW.  THE NB COLUMN VECTORS OF B */
/*     REPRESENT RIGHT-SIDE VECTORS FOR NB DISTINCT LINEAR LEAST SQUARES */
/*     PROBLEMS. */

/*     THIS SET OF PROBLEMS CAN ALSO BE WRITTEN AS THE MATRIX LEAST */
/*     SQUARES PROBLEM */

/*                       AX = B, */

/*     WHERE X IS THE N BY NB SOLUTION MATRIX. */

/*     NOTE THAT IF B IS THE M BY M IDENTITY MATRIX, THEN X WILL BE THE */
/*     PSEUDO-INVERSE OF A. */

/*     THIS SUBROUTINE FIRST TRANSFORMS THE AUGMENTED MATRIX (A B) TO A */
/*     MATRIX (R C) USING PREMULTIPLYING HOUSEHOLDER TRANSFORMATIONS WITH */
/*     COLUMN INTERCHANGES.  ALL SUBDIAGONAL ELEMENTS IN THE MATRIX R ARE */
/*     ZERO AND ITS DIAGONAL ELEMENTS SATISFY */

/*                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)), */

/*                       I = 1,...,L-1, WHERE */

/*                       L = MIN(M,N). */

/*     THE SUBROUTINE WILL COMPUTE AN INTEGER, KRANK, EQUAL TO THE NUMBER */
/*     OF DIAGONAL TERMS OF R THAT EXCEED TAU IN MAGNITUDE.  THEN A */
/*     SOLUTION OF MINIMUM EUCLIDEAN LENGTH IS COMPUTED USING THE FIRST */
/*     KRANK ROWS OF (R C). */

/*     TO BE SPECIFIC WE SUGGEST THAT THE USER CONSIDER AN EASILY */
/*     COMPUTABLE MATRIX NORM, SUCH AS, THE MAXIMUM OF ALL COLUMN SUMS OF */
/*     MAGNITUDES. */

/*     NOW IF THE RELATIVE UNCERTAINTY OF B IS EPS, (NORM OF UNCERTAINTY/ */
/*     NORM OF B), IT IS SUGGESTED THAT TAU BE SET APPROXIMATELY EQUAL TO */
/*     EPS*(NORM OF A). */

/*     THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST.. */
/*     A(MDA,N),(B(MDB,NB) OR B(M)),RNORM(NB),H(N),G(N),IP(N).  THIS */
/*     PERMITS THE SOLUTION OF A RANGE OF PROBLEMS IN THE SAME ARRAY */
/*     SPACE. */

/*     THE ENTIRE SET OF PARAMETERS FOR HFTI ARE */

/*     INPUT.. */

/*     A(*,*),MDA,M,N    THE ARRAY A(*,*) INITIALLY CONTAINS THE M BY N */
/*                       MATRIX A OF THE LEAST SQUARES PROBLEM AX = B. */
/*                       THE FIRST DIMENSIONING PARAMETER OF THE ARRAY */
/*                       A(*,*) IS MDA, WHICH MUST SATISFY MDA.GE.M */
/*                       EITHER M.GE.N OR M.LT.N IS PERMITTED.  THERE */
/*                       IS NO RESTRICTION ON THE RANK OF A.  THE */
/*                       CONDITION MDA.LT.M IS CONSIDERED AN ERROR. */

/*     B(*),MDB,NB       IF NB = 0 THE SUBROUTINE WILL PERFORM THE */
/*                       ORTHOGONAL DECOMPOSITION BUT WILL MAKE NO */
/*                       REFERENCES TO THE ARRAY B(*).  IF NB.GT.0 */
/*                       THE ARRAY B(*) MUST INITIALLY CONTAIN THE M BY */
/*                       NB MATRIX B OF THE LEAST SQUARES PROBLEM AX = */
/*                       B.  IF NB.GE.2 THE ARRAY B(*) MUST BE DOUBLY */
/*                       SUBSCRIPTED WITH FIRST DIMENSIONING PARAMETER */
/*                       MDB.GE.MAX(M,N).  IF NB = 1 THE ARRAY B(*) MAY */
/*                       BE EITHER DOUBLY OR SINGLY SUBSCRIPTED.  IN */
/*                       THE LATTER CASE THE VALUE OF MDB IS ARBITRARY */
/*                       BUT IT SHOULD BE SET TO SOME VALID INTEGER */
/*                       VALUE SUCH AS MDB = M. */

/*                       THE CONDITION OF NB.GT.1.AND.MDB.LT. MAX(M,N) */
/*                       IS CONSIDERED AN ERROR. */

/*     TAU               ABSOLUTE TOLERANCE PARAMETER PROVIDED BY USER */
/*                       FOR PSEUDORANK DETERMINATION. */

/*     H(*),G(*),IP(*)   ARRAYS OF WORKING SPACE USED BY HFTI. */

/*     OUTPUT.. */

/*     A(*,*)            THE CONTENTS OF THE ARRAY A(*,*) WILL BE */
/*                       MODIFIED BY THE SUBROUTINE.  THESE CONTENTS */
/*                       ARE NOT GENERALLY REQUIRED BY THE USER. */

/*     B(*)              ON RETURN THE ARRAY B(*) WILL CONTAIN THE N BY */
/*                       NB SOLUTION MATRIX X. */

/*     KRANK             SET BY THE SUBROUTINE TO INDICATE THE */
/*                       PSEUDORANK OF A. */

/*     RNORM(*)          ON RETURN, RNORM(J) WILL CONTAIN THE EUCLIDEAN */
/*                       NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM */
/*                       DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY */
/*                       B(*,*) FOR J = 1,...,NB. */

/*     H(*),G(*)         ON RETURN THESE ARRAYS RESPECTIVELY CONTAIN */
/*                       ELEMENTS OF THE PRE- AND POST-MULTIPLYING */
/*                       HOUSEHOLDER TRANSFORMATIONS USED TO COMPUTE */
/*                       THE MINIMUM EUCLIDEAN LENGTH SOLUTION. */

/*     IP(*)             ARRAY IN WHICH THE SUBROUTINE RECORDS INDICES */
/*                       DESCRIBING THE PERMUTATION OF COLUMN VECTORS. */
/*                       THE CONTENTS OF ARRAYS H(*),G(*) AND IP(*) */
/*                       ARE NOT GENERALLY REQUIRED BY THE USER. */

/*     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV) */
    /* Subroutine */ int h12_(integer *mode, integer *lpivot, integer *l1, 
    integer *m, double *u, integer *iue, double *up, double *c__, integer *ice, 
    integer *icv, integer *ncv);
    
/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (START CHANGES AT LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SDOT/DDOT/,/ABS,/DABS,/, */
/*     /SSWAP/DSWAP/,/SQRT/DSQRT/,/ABS(/ DABS(/,/AMAX1/DMAX1/, */
/*     /SAXPY/DAXPY/,/E0/D0/ */


/*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 */
/*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */

/*     MODIFIED AT SANDIA LABS., MAY 1977, TO -- */

/*     1)  REMOVE DOUBLE PRECISION ACCUMULATION, AND */
/*     2)  INCLUDE USAGE OF THE BASIC LINEAR ALGEBRA PACKAGE FOR */
/*         VECTORS LONGER THAN A PARTICULAR THRESHOLD. */

/*     CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
/*     HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B */

/*     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 . */
/*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
/*     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
/*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M */
/*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP    ON ENTRY TO H1 U() CONTAINS THE PIVOT VECTOR. */
/*                   IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS. */
/*                                       ON EXIT FROM H1 U() AND UP */
/*                   CONTAIN QUANTITIES DEFINING THE VECTOR U OF THE */
/*                   HOUSEHOLDER TRANSFORMATION.   ON ENTRY TO H2 U() */
/*                   AND UP SHOULD CONTAIN QUANTITIES PREVIOUSLY COMPUTED */
/*                   BY H1.  THESE WILL NOT BE MODIFIED BY H2. */
/*     C()    ON ENTRY TO H1 OR H2 C() CONTAINS A MATRIX WHICH WILL BE */
/*            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER */
/*            TRANSFORMATION IS TO BE APPLIED.  ON EXIT C() CONTAINS THE */
/*            SET OF TRANSFORMED VECTORS. */
/*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
/*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
/*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0 */
/*            NO OPERATIONS WILL BE DONE ON C(). */

    /* Subroutine */ int srotmg_(double *sd1, double *sd2, double *sx1, double *sy1, double 
    *sparam);
    
/*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/*     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)* */
/*     SY2)**T. */
/*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */

/*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */

/* END OF ABSTRACT */   

    /* Subroutine */ int scopy_(integer *n, double *sx, integer *incx, double *sy, 
    integer *incy);
    
/*     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY. */
/* END OF ABSTRACT */

    /* Subroutine */ int sswap_(integer *, double *, integer *, double *, 
        integer *);
               
/*     INTERCHANGE SINGLE PRECISION SX AND SINGLE PRECISION SY. */
/* END OF ABSTRACT */

    /* Subroutine */ int sscal_(integer *n, double *sa, double *sx, integer *incx);

/*     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX. */
/* END OF ABSTRACT */   

    /* Subroutine */ int saxpy_(integer *n, double *sa, double *sx, integer *incx, 
    double *sy, integer *incy);

/*     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY. */    

    /* Subroutine */ int srotm_(integer *n, double *sx, integer *incx, double *sy, 
    integer *incy, double *sparam);
    
/*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */

/*     (SX(1)     SX(N)) */
/*     (      ...      ) */
/*     (SY(1)     SY(N)) */

/*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */

/*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */   

    /* Subroutine */ int drotmg_(doublereal *dd1, doublereal *dd2, doublereal *
    dx1, doublereal *dy1, doublereal *dparam);
    
/* ***BEGIN PROLOGUE  DROTMG */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN */
/* ***PURPOSE */
/*    CONSTRUCT D.P. MODIFIED GIVENS TRANSFORMATION */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*      DD1  DOUBLE PRECISION SCALAR */
/*      DD2  DOUBLE PRECISION SCALAR */
/*      DX1  DOUBLE PRECISION SCALAR */
/*      DX2  DOUBLE PRECISION SCALAR */
/*   DPARAM  D.P. 5-VECTOR. DPARAM(1)=DFLAG DEFINED BELOW. */
/*             ELEMENTS 2-5  DEFINE THE TRANSFORMATION MATRIX H. */

/*     --OUTPUT-- */
/*      DD1  CHANGED TO REPRESENT THE EFFECT OF THE TRANSFORMATION */
/*      DD2  CHANGED TO REFLECT THE TRANSFORMATION */
/*      DX1  CHANGED TO REFLECT THE TRANSFORMATION */
/*      DX2  UNCHANGED */

/*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/*     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)* */
/*     DY2)**T. */
/*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */
/*     LOCATIONS 2-5 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22 */
/*     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE */
/*     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.) */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS   198   WITH     69  COMMENTS */
/* ***END PROLOGUE */   

    /* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy);
    
    /* ***BEGIN PROLOGUE  DCOPY */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*   D.P. VECTOR COPY Y = X */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */
/*       DY  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCY  STORAGE SPACING BETWEEN ELEMENTS OF DY */

/*     --OUTPUT-- */
/*       DY  COPY OF VECTOR DX (UNCHANGED IF N.LE.0) */

/*     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY. */
/*     FOR I = 0 TO N-1, COPY DX(LX+I*INCX) TO DY(LY+I*INCY), */
/*     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS */
/*     DEFINED IN A SIMILAR WAY USING INCY. */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    88   WITH     49  COMMENTS */
/* ***END PROLOGUE */

    /* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy);
    
/* ***BEGIN PROLOGUE  DSWAP */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*     INTERCHANGE D.P. VECTORS */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */
/*       DY  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCY  STORAGE SPACING BETWEEN ELEMENTS OF DY */

/*     --OUTPUT-- */
/*       DX  INPUT VECTOR DY (UNCHANGED IF N.LE.0) */
/*       DY  INPUT VECTOR DX (UNCHANGED IF N.LE.0) */

/*     INTERCHANGE DOUBLE PRECISION DX AND DOUBLE PRECISION DY. */
/*     FOR I = 0 TO N-1, INTERCHANGE  DX(LX+I*INCX) AND DY(LY+I*INCY), */
/*     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS */
/*     DEFINED IN A SIMILAR WAY USING INCY. */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    97   WITH     50  COMMENTS */
/* ***END PROLOGUE */   

    /* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx);
    
/* ***BEGIN PROLOGUE  DSCAL */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*    D.P. VECTOR SCALE X = A*X */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DA  DOUBLE PRECISION SCALE FACTOR */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */

/*     --OUTPUT-- */
/*       DX  DOUBLE PRECISION RESULT (UNCHANGED IF N.LE.0) */

/*     REPLACE DOUBLE PRECISION DX BY DOUBLE PRECISION DA*DX. */
/*     FOR I = 0 TO N-1, REPLACE DX(1+I*INCX) WITH  DA * DX(1+I*INCX) */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    68   WITH     43  COMMENTS */
/* ***END PROLOGUE */   

    /* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx, doublereal *dy, integer *incy);

/* ***BEGIN PROLOGUE  DAXPY */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*  D.P COMPUTATION Y = A*X + Y */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DA  DOUBLE PRECISION SCALAR MULTIPLIER */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */
/*       DY  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCY  STORAGE SPACING BETWEEN ELEMENTS OF DY */

/*     --OUTPUT-- */
/*       DY  DOUBLE PRECISION RESULT (UNCHANGED IF N.LE.0) */

/*     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY. */
/*     FOR I = 0 TO N-1, REPLACE  DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) + */
/*       DY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N */
/*       AND LY IS DEFINED IN A SIMILAR WAY USING INCY. */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED   NONE */
/* ***CARD COUNT IS    86   WITH     50  COMMENTS */
/* ***END PROLOGUE */

    /* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy, doublereal *dparam);
    
/* ***BEGIN PROLOGUE  DROTM */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79, LAWSON C. (JPL), HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*   APPLY D.P. MODIFIED GIVENS TRANSFORMATION */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */
/*       DY  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCY  STORAGE SPACING BETWEEN ELEMENTS OF DY */
/*   DPARAM  5-ELEMENT D.P. VECTOR. DPARAM(1) IS DFLAG DESCRIBED BELOW */
/*            ELEMENTS 2-5 FORM THE TRANSFORMATION MATRIX H. */

/*     --OUTPUT-- */
/*       DX  ROTATED VECTOR (UNCHANGED IF N.LE.0) */
/*       DY  ROTATED VECTOR (UNCHANGED IF N.LE.0) */

/*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */

/*     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN */
/*     (DY**T) */

/*     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/*     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY. */
/*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */
/*     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM. */


/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS   142   WITH     54  COMMENTS */
/* ***END PROLOGUE */   

    /* Subroutine */ int fdump_();

/*     ABSTRACT */
/*        ***NOTE*** MACHINE DEPENDENT ROUTINE */
/*        FDUMP IS INTENDED TO BE REPLACED BY A LOCALLY WRITTEN */
/*        VERSION WHICH PRODUCES A SYMBOLIC DUMP.  FAILING THIS, */
/*        IT SHOULD BE REPLACED BY A VERSION WHICH PRINTS THE */
/*        SUBPROGRAM NESTING LIST.  NOTE THAT THIS DUMP MUST BE */
/*        PRINTED ON EACH OF UP TO FIVE FILES, AS INDICATED BY THE */
/*        XGETUA ROUTINE.  SEE XSETUA AND XGETUA FOR DETAILS. */

/*     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE */
/* END OF ABSTRACT */
/*     LATEST REVISION ---  23 MAY 1979 */  


    /* Subroutine */ int xerror_(const char *messg, integer *nmessg, integer *nerr, 
    integer *level);
/*     ABSTRACT */
/*        XERROR PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER */
/*        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE */
/*        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL. */
/*        (SEE SUBROUTINE XSETF FOR DETAILS.) */

/*     DESCRIPTION OF PARAMETERS */
/*      --INPUT-- */
/*        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED, CONTAINING */
/*                NO MORE THAN 72 CHARACTERS. */
/*        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG. */
/*        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE. */
/*                NERR MUST NOT BE ZERO. */
/*        LEVEL - ERROR CATEGORY. */
/*                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR. */
/*                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS */
/*                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.) */
/*                =0 MEANS THIS IS A WARNING MESSAGE ONLY. */
/*                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE */
/*                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY */
/*                   TIMES THIS CALL IS EXECUTED. */

/*     EXAMPLES */
/*        CALL XERROR(23HSMOOTH -- NUM WAS ZERO.,23,1,2) */
/*        CALL XERROR(43HINTEG  -- LESS THAN FULL ACCURACY ACHIEVED., */
/*                    43,2,1) */
/*        CALL XERROR(65HROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL */
/*    1 FULLY COLLAPSED.,65,3,0) */
/*        CALL XERROR(39HEXP    -- UNDERFLOWS BEING SET TO ZERO.,39,1,-1) */

/*     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE */
/* END OF ABSTRACT */
/*     REVISED BY K HASKELL TO CHECK INPUT ARGS, 2/18/80 */         

    /* Subroutine */ int xerrwv_(const char *messg, integer *nmessg, integer *nerr, 
    integer *level, integer *ni, integer *i1, integer *i2, integer *nr, 
    double *r1, double *r2);        
/*     ABSTRACT */
/*        XERRWV PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER */
/*        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE */
/*        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL. */
/*        (SEE SUBROUTINE XSETF FOR DETAILS.) */
/*        IN ADDITION, UP TO TWO INTEGER VALUES AND TWO double */
/*        VALUES MAY BE PRINTED ALONG WITH THE MESSAGE. */

/*     DESCRIPTION OF PARAMETERS */
/*      --INPUT-- */
/*        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED. */
/*        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG. */
/*        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE. */
/*                NERR MUST NOT BE ZERO. */
/*        LEVEL - ERROR CATEGORY. */
/*                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR. */
/*                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS */
/*                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.) */
/*                =0 MEANS THIS IS A WARNING MESSAGE ONLY. */
/*                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE */
/*                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY */
/*                   TIMES THIS CALL IS EXECUTED. */
/*        NI    - NUMBER OF INTEGER VALUES TO BE PRINTED. (O TO 2) */
/*        I1    - FIRST INTEGER VALUE. */
/*        I2    - SECOND INTEGER VALUE. */
/*        NR    - NUMBER OF double VALUES TO BE PRINTED. (0 TO 2) */
/*        R1    - FIRST double VALUE. */
/*        R2    - SECOND double VALUE. */

/*     EXAMPLES */
/*        CALL XERRWV(29HSMOOTH -- NUM (=I1) WAS ZERO.,29,1,2, */
/*    1   1,NUM,0,0,0.,0.) */
/*        CALL XERRWV(54HQUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM */
/*    1 (R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN) */

/*     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE */
/* END OF ABSTRACT */
/*     LATEST REVISION ---  19 MAR 1980 */
/*     REVISED BY K HASKELL TO CHECK INPUT ARGS, 2/18/80 */

    /* Subroutine */ int clstp_(integer *klog, double *cond, integer *istat);

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (BEGIN THE CHANGES AT THE LINE WITH C++ IN COLS. 1-3.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/,/SDOT/DDOT/, */
/*     /SNRM2/DNRM2/,/SQRT/DSQRT/,/E0/D0/,/SSCAL/DSCAL/,/SAXPY/DAXPY/, */
/*     /SRELPR/DRELPR/,/SSWAP/DSWAP/ */

/*     REVISED 820305-2000 */
/*     REVISED YYMMDD-HHMM */

/*     THIS SUBROUTINE EXERCISES MOST OF THE MATHEMATICAL FEATURES OF THE */
/*     CONSTRAINED LEAST SQUARES SUBPROGRAMS WNNLS( ) AND LSEI( ). */
/*     THE PROBLEM THAT IS SOLVED HERE IS OF THE FORM */

/*               A*X=B (LEAST SQUARES, A MA BY N), */


/*               SUBJECT TO CONSTRAINTS */

/*               E*X=F      (CONSTRAINT EQUATIONS, ME BY N), */
/*          AND  G*X .GE. H (INEQUALITY CONSTRAINTS, MG BY N). */

/*     THE CLASS OF PROBLEMS THAT IS SOLVED HERE IS GENERATED WITH */
/*     HADAMARD MATRICES OF ORDER=POWER OF 2.  EACH OF THE MATRICES */
/*     A,E, AND G HAVE A SPECIFIED CONDITION NUMBER.  FOR EXAMPLE */
/*     A=HADAMARD MATRIX * DIAGONAL MATRIX * HADAMARD MATRIX. */
/*     DIAGONAL TERMS OF THE DIAGONAL MATRIX ARE CHOSEN SO THAT A */
/*     HAS A GIVEN CONDITION NUMBER.  THE MATRICES E AND G ARE */
/*     CONSTRUCTED IN SIMILIAR WAYS.  FURTHER, THE PROBLEM IS CONSTRUCTED */
/*     SO THAT THE TRUE SOLUTION IS X=(1,...,1) (TRANSPOSED). */
/*     THIS REQUIRES COMPUTING THE RIGHT HAND SIDE VECTORS B,F AND */
/*     H.  THE VECTOR B=A*X+COMPONENT ORTHOGONAL TO COL. SPACE OF */
/*     A, F=E*X, AND H=G*H-SLACK COMPONENTS.  THESE SLACK COMPONENTS */
/*     ARE CHOSEN SO THAT THE FIRST MI OF THE INEQUALITIES ARE */
/*     STRICT INEQUALITIES. */

/*     THE PROBLEMS DIMENSIONS ARE SPECIFIED BY */

/*                     MA = 2**KA */
/*                     ME = 2**KE */
/*                     MG = 2**KG */
/*                     MI = 2**KI */
/*                     N  = 2**KN */

/*     WHERE KA, KE, KG, KI, AND KN ARE INPUT TO THE SUBROUTINE AS */
/*     DISCUSSED BELOW. */

/*     THE SUBROUTINE ARGUMENTS ARE AS FOLLOWS */

/*     I N P U T */

/*     KLOG(*)    - AN INTEGER ARRAY WHOSE DIMENSION IS AT LEAST 5.  THE */
/*                  ENTRIES CORRESPOND TO THE POWERS OF 2 NECESSARY */
/*                  TO SPECIFY THE PROBLEM DIMENSIONS.  REFERRING TO */
/*                  THE ABOVE DISCUSSION, THE ENTRIES OF KLOG(*) */
/*                  SHOULD BE SET AS FOLLOWS */

/*                       KLOG(1) = KA */
/*                       KLOG(2) = KE */
/*                       KLOG(3) = KG */
/*                       KLOG(4) = KI */
/*                       KLOG(5) = KN */

/*                  IF KA, KE, KG, OR KI IS LESS THAN ZERO, THE */
/*                  CORRESPONDING DIMENSION WILL BE SET TO ZERO. */

/*                  KN.LT.0 WILL CAUSE THE SUBROUTINE TO SIMPLY RETURN. */

/*     COND(*)    - AN ARRAY WHOSE DIMENSION IS AT LEAST 3.  THE */
/*                  ENTRIES COND(1), COND(2), AND COND(3) RESPECTIVELY */
/*                  SPECIFY THE CONDITION NUMBER FOR THE LEAST SQUARES */
/*                  MATRIX, THE EQUALITY CONSTRAINT MATRIX, AND THE */
/*                  INEQUALITY CONSTRAINT MATRIX. */

/*     O U T P U T */

/*     ISTAT      - AN INTEGER FLAG WHICH INDICATES WHETHER THE SOLUTION */
/*                  WAS CORRECTLY COMPUTED. */
/*                  =1 NEITHER WNNLS( ) NOR LSEI( ) PASSED THE TEST. */
/*                  =2 WNNLS( ) PASSED BUT LSEI( ) FAILED THE TEST. */
/*                  =3 LSEI( ) PASSED BUT WNNLS( ) FAILED THE TEST. */
/*                  =4 BOTH WNNLS( ) AND LSEI( ) PASSED THE TEST. */

/*     THE DIMENSION STATEMENTS BELOW ARE SET UP TO SOLVE PROBLEMS FOR */
/*     WHICH NONE OF THE ABOVE LOGARITHMS IS GREATER THAN 5.  TO CHANGE */
/*     THESE DIMENSIONS TO SOLVE LARGER PROBLEMS, USE THE FOLLOWING */
/*     FORMULAS */

/*     DIMENSION W(MA+ME+MG,N+MG+1),X(N+MG),HH(MMAX,MMAX),GG(MMAX,MMAX) */
/*     DIMENSION WORK(2*(ME+N)+K+(MG+2)*(N+7)),IWORK(ME+MA+2*MG+N) */
/*     DIMENSION SA(MIN(MA,N)),SE(MIN(ME,N)),SG(MIN(MG,N)) */

/*     WHERE MMAX = MAX(MA,ME,MG,N) */
/*           K    = MAX(MA+MG,N) */

/*     NOTE THAT IF THE DIMENSIONS ARE CHANGED, THE VALUES ASSIGNED TO */
/*     MDW, MDH, AND MDG BELOW MUST BE ALTERED APPROPRIATELY.  THESE */
/*     ARE THE RESPECTIVE ROW DIMENSIONS OF THE ARRAYS W(*,*), HH(*,*), */
/*     AND GG(*,*). */  
    
	int SelfTest();
private:

    inline doublereal diff_(double *x, double *y)  { return *x - *y; } /* diff_ */

/*     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO */
/*     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES. */
/*     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/. */
/*     (APPLY CHANGES TO ENTIRE PROGRAM UNIT.) */
/*     /double (12 BLANKS)/DOUBLE PRECISION/ */

/*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUNE 7 */
/*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */    
    
    doublereal snrm2_(integer *n, double *sx, integer *incx);
/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() double AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */
/* END OF ABSTRACT */   
    
    doublereal sasum_(integer *n, double *sx, integer *incx);
    
/*     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX. */
/* END OF ABSTRACT */   
    
    integer isamax_(integer *n, double *sx, integer *incx);
    
/*     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX. */
/* END OF ABSTRACT */
    
    doublereal sdot_(integer *n, double *sx, integer *incx, double *sy, integer *incy);

/*     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY. */    
    
    doublereal dnrm2_(integer *n, doublereal *dx, integer *incx);
    
/* ***BEGIN PROLOGUE  DNRM2 */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*    EUCLIDEAN LENGTH (L2 NORM) OF D.P. VECTOR */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */

/*     --OUTPUT-- */
/*    DNRM2  DOUBLE PRECISION RESULT (ZERO IF N.LE.0) */

/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() double AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS   151   WITH    105  COMMENTS */
/* ***END PROLOGUE */   
    
    doublereal dasum_(integer *n, doublereal *dx, integer *incx);
    
/* ***BEGIN PROLOGUE  DASUM */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*    SUM OF MAGNITUDES OF D.P. VECTOR COMPONENTS */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */

/*     --OUTPUT-- */
/*    DASUM  DOUBLE PRECISION RESULT (ZERO IF N.LE.0) */

/*     RETURNS SUM OF MAGNITUDES OF DOUBLE PRECISION DX. */
/*     DASUM = SUM FROM 0 TO N-1 OF DABS(DX(1+I*INCX)) */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    65   WITH     42  COMMENTS */
/* ***END PROLOGUE */
    
    integer idamax_(integer *n, doublereal *dx, integer *incx);

/* ***BEGIN PROLOGUE  IDAMAX */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*     FIND LARGEST COMPONENT OF D.P. VECTOR */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */

/*     --OUTPUT-- */
/*   IDAMAX  SMALLEST INDEX (ZERO IF N.LE.0) */

/*     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF DOUBLE PRECISION DX. */
/*     IDAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(DX(1-INCX+I*INCX) */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    66   WITH     39  COMMENTS */
/* ***END PROLOGUE */
    
    doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
    integer *incy);

/* ***BEGIN PROLOGUE  DDOT */
/* ***REVISION 3/1/80 */
/* ***CATEGORY NO. */
/* ***KEYWORD(S) */
/* ***AUTHOR--DATE WRITTEN  10/79,LAWSON C. (JPL),HANSON R. (SLA), */
/*                            KINCAID D. (U TEXAS), KROGH F. (JPL) */
/* ***PURPOSE */
/*   D.P. INNER PRODUCT OF D.P. VECTORS */
/* ***DESCRIPTION */
/*                B L A S  SUBPROGRAM */
/*    DESCRIPTION OF PARAMETERS */

/*     --INPUT-- */
/*        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S) */
/*       DX  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCX  STORAGE SPACING BETWEEN ELEMENTS OF DX */
/*       DY  DOUBLE PRECISION VECTOR WITH N ELEMENTS */
/*     INCY  STORAGE SPACING BETWEEN ELEMENTS OF DY */

/*     --OUTPUT-- */
/*     DDOT  DOUBLE PRECISION DOT PRODUCT (ZERO IF N.LE.0) */

/*     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY. */
/*     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY) */
/*     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS */
/*     DEFINED IN A SIMILAR WAY USING INCY. */

/* ***REFERENCE(S) */
/*  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T., */
/*   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, */
/*  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE, */
/*  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 */
/* ***ROUTINES CALLED  NONE */
/* ***CARD COUNT IS    84   WITH     49  COMMENTS */
/* ***END PROLOGUE */
        
    

doublereal ran_(integer *k);

/*     RANDOM NUMBER GENERATOR - BASED ON ALGORITHM 266 */
/*      BY PIKE AND HILL (MODIFIED BY HANSSON) */
/*      COLLECTED ALG. FROM CACM. */

/*     THIS SUBPROGRAM IS INTENDED FOR USE ON COMPUTERS WITH */
/*      FIXED POINT WORDLENGTH OF AT LEAST 29 BITS.  IT IS */
/*      BEST IF THE FLOATING POINT SIGNIFICAND HAS AT MOST */
/*      29 BITS. */
};
#endif LSEI_H
