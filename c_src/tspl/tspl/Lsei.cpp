#include "Lsei.h"

/* Table of constant values */

 integer c__67 = 67;
 integer c__1 = 1;
 integer c__0 = 0;
 integer c__68 = 68;
 integer c__36 = 36;
 integer c__54 = 54;
 integer c__2 = 2;
 integer c__38 = 38;
 integer c__52 = 52;
 integer c__70 = 70;
 integer c__72 = 72;
 integer c__44 = 44;
 integer c__39 = 39;
 integer c__53 = 53;
 integer c__31 = 31;
 integer c__49 = 49;
 logical c_false = FALSE_;
 logical c_true = TRUE_;
 integer c__4 = 4;
 integer c__6 = 6;
 integer c__17 = 17;
 integer c__33 = 33;
 integer c__29 = 29;
 integer c__23 = 23;
 integer c__28 = 28;
 integer c__32 = 32;
 double c_b913 = 0;
 integer c__57 = 57;
 integer c__13 = 13;
 integer c__35 = 35;
 integer c__5 = 5;
 integer c__3 = 3;
 integer c__34 = 34;
 integer c_n1 = -1;
 integer c__4000 = 4000;

/*    extern "C" {
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *,
        ftnlen), e_rsfe(), s_wsfe(cilist *), e_wsfe(), f_clos(cllist *);
    //int s_stop(char *, ftnlen);
    }
*/

int lsei::lsei_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *mg, integer *n, double *prgopt, double *x, double *rnorme, double *
    rnorml, integer *mode)
{
    /****************** seting dimentions of matricies ***************************/
    /*     WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2) */
    /*     ABOVE, K=MAX(MA+MG,N). */    
    int k_temp =  (max(*ma+*mg, *n));   
    double *ws = new double[2*(*me+*n)+ k_temp + (*mg + 2)*(*n + 7)+1];
    integer *ip = new integer[*mg +2*(*n)+3];
    /*****************************************************************************/
    int val = lsei_(w, mdw, me, ma, mg, n , prgopt, x, rnorme, rnorml, mode, ws, ip);
    delete[] ws;
    delete[] ip;
    return val;
}

/* Subroutine */ 
int lsei::lsei_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *mg, integer *n, double *prgopt, double *x, double *rnorme, double *
    rnorml, integer *mode, double *ws, integer *ip)
    
{
    /* Initialized data */

    double zero = 0.;
    double srelpr = 0.;
    double half = .5;
    double one = 1.;

    /* Format strings */
//    char fmt_90[] = "";

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3;
    double r__1, r__2;

    /* Builtin functions */
 //   double sqrt(doublereal);

    /* Local variables */
    integer lchk, mend, link, imax, last, nerr;    
    double size;
    integer iopt, next, nopt, igo990;
    integer i__, j, k, l, m;
    double t;
    integer mdeqc;
    integer nlink;
    double enorm, fnorm, rnmax, snmax;
    double xnrme, dummy;
    integer n1, n2;
    double xnorm;
    integer mapke1;
    integer ii;
    double rb;
    integer kk, jj;
    double uj, rn, sn, vj, up;
    integer kranke, ntimes, jp1, np1;
    double gam;
    logical cov;
    double tau;
    integer key, mep1;
    /* Assigned format variables */
//    char *igo990_fmt;

    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;
    int k_temp =  (max(*ma+*mg, *n));   
    ip[1]=2*(*me+*n)+ k_temp + (*mg + 2)*(*n + 7)+1;
    ip[2]=*mg +2*(*n)+3;
//  ip[3]=

    /* Function Body */

//  goto L30;

/*     CHECK THAT ENOUGH STORAGE WAS ALLOCATED IN WS(*) AND IP(*). */
    if (! (ip[1] > 0)) {
    goto L20;
    }
/* Computing MAX */
    i__1 = *ma + *mg;
    lchk = ((*me + *n )<< 1) + max(i__1,*n) + (*mg + 2) * (*n + 7);
    if (! (ip[1] < lchk)) {
    goto L10;
    }
    *mode = 4;
    nerr = 2;
    iopt = 1;
   xerrwv_("LSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR WS(*), NEED LW=I1 BELOW", &c__67,
       &nerr, &iopt, &c__1, &lchk, &c__0, &c__0, &dummy, &dummy/*, (ftnlen)67*/);
    return 0;
L10:
L20:
    if (! (ip[2] > 0)) {
    goto L40;
    }
    lchk = *mg + (*n << 1) + 2;
    if (! (ip[2] < lchk)) {
    goto L30;
    }
    *mode = 4;
    nerr = 2;
    iopt = 1;
    xerrwv_("LSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR IP(*), NEED LIP=I1 BELOW", &c__68, &nerr, 
        &iopt, &c__1, &lchk, &c__0, &c__0 , &dummy, &dummy/*, (ftnlen)68*/);
    return 0;
L30:

/*     COMPUTE MACHINE PRECISION=SRELPR ONLY WHEN NECESSARY. */
L40:
    if (! (M_precision == zero)) goto L70;
    M_precision = one;
L50:
    if (one + M_precision == one) goto L60;
    M_precision *= half;
    goto L50;
L60:
    M_precision *= 10;
L70:
    srelpr=M_precision;
 //   srelpr *= 10;

/*     COMPUTE NUMBER OF POSSIBLE RIGHT MULTIPLYING HOUSEHOLDER */
/*     TRANSFORMATIONS. */

    m = *me + *ma + *mg;
    *mode = 0;
    if (*n <= 0 || m <= 0) {
    return 0;
    }
    if (! (*mdw < m)) {
    goto L80;
    }
    nerr = 1;
    iopt = 1;
    xerror_("LSEI( ), MDW.LT.ME+MA+MG IS AN ERROR", &c__36, &nerr, &iopt/*, (ftnlen)36*/);
    *mode = 4;
    return 0;
L80:
    np1 = *n + 1;
    kranke = min(*me,*n);
    n1 = (kranke << 1) + 1;
    n2 = n1 + *n;

/*     PROCESS-OPTION-VECTOR */
    igo990 = 0;
//    igo990_fmt = fmt_90;
    goto L480;
L90:
    if (! (cov && *mdw < *n)) {
    goto L100;
    }
    nerr = 2;
    iopt = 1;
    xerror_("LSEI( ), MDW.LT.N, WHEN COV MATRIX NEEDED, IS AN ERROR", &c__54, 
        &nerr, &iopt/*, (ftnlen)54*/);
    *mode = 4;
    return 0;
L100:
    l = kranke;

/*     COMPUTE NORM OF EQUALITY CONSTRAINT MATRIX AND RT SIDE. */
    enorm = zero;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
    r__1 = enorm; r__2 = sasum_(me, &w[j * w_dim1 + 1], &c__1); ////
    enorm = dmax(r__1,r__2);
/* L110: */
    }
    fnorm = sasum_(me, &w[np1 * w_dim1 + 1], &c__1);
    if (! (l > 0)) {
    goto L190;
    }
    snmax = zero;
    rnmax = zero;
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     COMPUTE MAXIMUM RATIO OF VECTOR LENGTHS. PARTITION */
/*     IS AT COL. I. */
    i__2 = *me;
    for (k = i__; k <= i__2; ++k) {
        i__3 = *n - i__ + 1;
        sn = sdot_(&i__3, &w[k + i__ * w_dim1], mdw, &w[k + i__ * w_dim1], mdw);
        i__3 = i__ - 1;
        rn = sdot_(&i__3, &w[k + w_dim1], mdw, &w[k + w_dim1], mdw);
        if (! (rn == zero && sn > snmax)) {
        goto L120;
        }
        snmax = sn;
        imax = k;
        goto L140;
L120:
        if (! (k == i__ || sn * rnmax > rn * snmax)) {
        goto L130;
        }
        snmax = sn;
        rnmax = rn;
        imax = k;
L130:
L140:
/* L150: */
        ;
    }

/*     INTERCHANGE ROWS IF NECESSARY. */
    if (i__ != imax) {
        sswap_(&np1, &w[i__ + w_dim1], mdw, &w[imax + w_dim1], mdw);
    }
/* Computing 2nd power */
    r__1 = tau;
    if (! (snmax > r__1 * r__1 * rnmax)) {
        goto L160;
    }

/*     ELIMINATE ELEMS I+1,...,N IN ROW I. */
    i__2 = i__ + 1;
    i__3 = m - i__;
    h12_(&c__1, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[i__], &w[i__ + 1 + w_dim1], mdw, &c__1, &i__3);
    goto L170;
L160:
    kranke = i__ - 1;
    goto L200;
L170:
/* L180: */
    ;
    }
L190:
L200:

/*     SAVE DIAG. TERMS OF LOWER TRAP. MATRIX. */
    i__1 = *mdw + 1;
    scopy_(&kranke, &w[w_offset], &i__1, &ws[kranke + 1], &c__1);

/*     USE HOUSEHOLDER TRANS FROM LEFT TO ACHIEVE KRANKE BY KRANKE UPPER */
/*     TRIANGULAR FORM. */
    if (! (kranke > 0 && kranke < *me)) {
    goto L220;
    }
    i__1 = kranke;
    for (kk = 1; kk <= i__1; ++kk) {
    k = kranke + 1 - kk;

/*     APPLY TRANFORMATION TO MATRIX COLS. 1,...,K-1. */
    i__2 = kranke + 1;
    i__3 = k - 1;
    h12_(&c__1, &k, &i__2, me, &w[k * w_dim1 + 1], &c__1, &up, &w[  w_offset], &c__1, mdw, &i__3);

/*     APPLY TO RT SIDE VECTOR. */
    i__2 = kranke + 1;
    h12_(&c__2, &k, &i__2, me, &w[k * w_dim1 + 1], &c__1, &up, &w[np1 * w_dim1 + 1], &c__1, &c__1, &c__1);
/* L210: */
    }
L220:
    if (! (kranke > 0)) {
    goto L240;
    }

/*     SOLVE FOR VARIABLES 1,...,KRANKE IN NEW COORDINATES. */
    scopy_(&kranke, &w[np1 * w_dim1 + 1], &c__1, &x[1], &c__1);
    i__1 = kranke;
    for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = i__ - 1;
    x[i__] = (x[i__] - sdot_(&i__2, &w[i__ + w_dim1], mdw, &x[1], &c__1))   / w[i__ + i__ * w_dim1];
/* L230: */
    }

/*     COMPUTE RESIDUALS FOR REDUCED PROBLEM. */
L240:
    mep1 = *me + 1;
    *rnorml = zero;
    if (! (*me < m)) {
    goto L270;
    }
    i__1 = m;
    for (i__ = mep1; i__ <= i__1; ++i__) {
    w[i__ + np1 * w_dim1] -= sdot_(&kranke, &w[i__ + w_dim1], mdw, &x[1], &c__1);
    sn = sdot_(&kranke, &w[i__ + w_dim1], mdw, &w[i__ + w_dim1], mdw);
    i__2 = *n - kranke;
    rn = sdot_(&i__2, &w[i__ + (kranke + 1) * w_dim1], mdw, &w[i__ + (kranke + 1) * w_dim1], mdw);
/* Computing 2nd power */
    r__1 = tau;
    if (! (rn <= r__1 * r__1 * sn && kranke < *n)) {
        goto L250;
    }
    w[i__ + (kranke + 1) * w_dim1] = zero;
    i__2 = *n - kranke;
    scopy_(&i__2, &w[i__ + (kranke + 1) * w_dim1], &c__0, &w[i__ + (    kranke + 1) * w_dim1], mdw);
L250:
/* L260: */
    ;
    }

/*     COMPUTE EQUAL. CONSTRAINT EQUAS. RESIDUAL LENGTH. */
L270:
    i__1 = *me - kranke;
    *rnorme = snrm2_(&i__1, &w[kranke + 1 + np1 * w_dim1], &c__1);

/*     MOVE REDUCED PROBLEM DATA UPWARD IF KRANKE.LT.ME. */
    if (! (kranke < *me)) {
    goto L290;
    }
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
    i__2 = m - *me;
    scopy_(&i__2, &w[*me + 1 + j * w_dim1], &c__1, &w[kranke + 1 + j *  w_dim1], &c__1);
/* L280: */
    }

/*     COMPUTE SOLN OF REDUCED PROBLEM. */
L290:
    i__1 = *n - kranke;
    lsi_(&w[kranke + 1 + (kranke + 1) * w_dim1], mdw, ma, mg, &i__1, &prgopt[1], &x[kranke + 1], rnorml, mode, &ws[n2], &ip[2]);
    if (! (*me > 0)) {
    goto L330;
    }

/*     TEST FOR CONSISTENCY OF EQUALITY CONSTRAINTS. */
    mdeqc = 0;
    xnrme = sasum_(&kranke, &w[np1 * w_dim1 + 1], &c__1);
    if (*rnorme > tau * (enorm * xnrme + fnorm)) {
    mdeqc = 1;
    }
    *mode += mdeqc;

/*     CHECK IF SOLN TO EQUAL. CONSTRAINTS SATISFIES INEQUAL. */
/*     CONSTRAINTS WHEN THERE ARE NO DEGREES OF FREEDOM LEFT. */
    if (! (kranke == *n && *mg > 0)) {
    goto L320;
    }
    xnorm = sasum_(n, &x[1], &c__1);
    mapke1 = *ma + kranke + 1;
    mend = *ma + kranke + *mg;
    i__1 = mend;
    for (i__ = mapke1; i__ <= i__1; ++i__) {
    size = sasum_(n, &w[i__ + w_dim1], mdw) * xnorm + (r__1 = w[i__ + np1   * w_dim1], dabs(r__1));
    if (! (w[i__ + np1 * w_dim1] > tau * size)) {
        goto L300;
    }
    *mode += 2;
    goto L450;
L300:
/* L310: */
    ;
    }
L320:
L330:
    if (! (kranke > 0)) {
    goto L420;
    }

/*     REPLACE DIAG. TERMS OF LOWER TRAP. MATRIX. */
    i__1 = *mdw + 1;
    scopy_(&kranke, &ws[kranke + 1], &c__1, &w[w_offset], &i__1);

/*     REAPPLY TRANS TO PUT SOLN IN ORIGINAL COORDINATES. */
    i__1 = kranke;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = kranke + 1 - ii;
    i__2 = i__ + 1;
    h12_(&c__2, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[i__], &x[1], &c__1, &c__1, &c__1);
/* L340: */
    }

/*     COMPUTE COV MATRIX OF EQUAL. CONSTRAINED PROBLEM. */
    if (! cov) {
    goto L410;
    }
    i__1 = kranke;
    for (jj = 1; jj <= i__1; ++jj) {
    j = kranke + 1 - jj;
    if (! (j < *n)) {
        goto L390;
    }
    rb = ws[j] * w[j + j * w_dim1];
    if (rb != zero) {
        rb = one / rb;
    }
    jp1 = j + 1;
    i__2 = *n;
    for (i__ = jp1; i__ <= i__2; ++i__) {
        i__3 = *n - j;
        w[i__ + j * w_dim1] = sdot_(&i__3, &w[i__ + jp1 * w_dim1], mdw, &w[j + jp1 * w_dim1], mdw) * rb;
/* L350: */
    }
    i__2 = *n - j;
    gam = sdot_(&i__2, &w[jp1 + j * w_dim1], &c__1, &w[j + jp1 * w_dim1], mdw) * rb;
    gam = half * gam;
    i__2 = *n - j;
    saxpy_(&i__2, &gam, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * w_dim1], &c__1);
    i__2 = *n;
    for (i__ = jp1; i__ <= i__2; ++i__) {
        i__3 = *n;
        for (k = i__; k <= i__3; ++k) {
        w[i__ + k * w_dim1] = w[i__ + k * w_dim1] + w[j + i__ * 
            w_dim1] * w[k + j * w_dim1] + w[i__ + j * w_dim1] * w[j + k * w_dim1];
        w[k + i__ * w_dim1] = w[i__ + k * w_dim1];
/* L360: */
        }
/* L370: */
    }
    uj = ws[j];
    vj = gam * uj;
    w[j + j * w_dim1] = uj * vj + uj * vj;
    i__2 = *n;
    for (i__ = jp1; i__ <= i__2; ++i__) {
        w[j + i__ * w_dim1] = uj * w[i__ + j * w_dim1] + vj * w[j + i__ * w_dim1];
/* L380: */
    }
    i__2 = *n - j;
    scopy_(&i__2, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * w_dim1], &c__1);
L390:
/* L400: */
    ;
    }
L410:

/*     APPLY THE SCALING TO THE COVARIANCE MATRIX. */
L420:
    if (! cov) {
    goto L440;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    l = n1 + i__;
    sscal_(n, &ws[l - 1], &w[i__ + w_dim1], mdw);
    sscal_(n, &ws[l - 1], &w[i__ * w_dim1 + 1], &c__1);
/* L430: */
    }
L440:
L450:

/*     RESCALE SOLN. VECTOR. */
    if (! (*mode <= 1)) {
    goto L470;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    l = n1 + j;
    x[j] *= ws[l - 1];
/* L460: */
    }
L470:
    ip[1] = kranke;
    ip[3] = ip[3] + (kranke << 1) + *n;
    return 0;
L480:
/*     TO PROCESS-OPTION-VECTOR */

/*     THE NOMINAL TOLERANCE USED IN THE CODE */
/*     FOR THE EQUALITY CONSTRAINT EQUATIONS. */
    tau = sqrt(srelpr);

/*     THE NOMINAL COLUMN SCALING USED IN THE CODE IS */
/*     THE IDENTITY SCALING. */
    ws[n1] = one;
    scopy_(n, &ws[n1], &c__0, &ws[n1], &c__1);

/*     NO COVARIANCE MATRIX IS NOMINALLY COMPUTED. */
    cov = FALSE_;

/*     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE. */
    nopt = 1000;
    ntimes = 0;

/*     DEFINE BOUND FOR POSITIVE VALUES OF LINK. */
    nlink = 100000;
    last = 1;
    link = prgopt[1];
    if (! (link <= 0 || link > nlink)) {
    goto L490;
    }
    nerr = 3;
    iopt = 1;
//    xerror_("LSEI( ) THE OPTION VECTOR IS UNDEFINED", &c__38, &nerr, &iopt, ( ftnlen)38);
    *mode = 4;
    return 0;
L490:
    if (! (link > 1)) {
    goto L540;
    }
    ++ntimes;
    if (! (ntimes > nopt)) {
    goto L500;
    }
    nerr = 3;
    iopt = 1;
    xerror_("LSEI( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING.", &c__52, &nerr, &iopt/*, (ftnlen)52*/);
    *mode = 4;
    return 0;
L500:
    key = prgopt[last + 1];
    if (key == 1) {
    cov = prgopt[last + 2] != zero;
    }
    if (! (key == 2 && prgopt[last + 2] != zero)) {
    goto L520;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    t = snrm2_(&m, &w[j * w_dim1 + 1], &c__1);
    if (t != zero) {
        t = one / t;
    }
    l = n1 + j;
    ws[l - 1] = t;
/* L510: */
    }
L520:
    if (key == 3) {
    scopy_(n, &prgopt[last + 2], &c__1, &ws[n1], &c__1);
    }
    if (key == 4) {
/* Computing MAX */
    r__1 = srelpr, r__2 = prgopt[last + 2];
    tau = dmax(r__1,r__2);
    }
    next = prgopt[link];
    if (! (next <= 0 || next > nlink)) {
    goto L530;
    }
    nerr = 3;
    iopt = 1;
    xerror_("LSEI( ) THE OPTION VECTOR IS UNDEFINED", &c__38, &nerr, &iopt/*, (ftnlen)38*/);
    *mode = 4;
    return 0;
L530:
    last = link;
    link = next;
    goto L490;
L540:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    l = n1 + j;
    sscal_(&m, &ws[l - 1], &w[j * w_dim1 + 1], &c__1);
/* L550: */
    }
    goto L560;
L560:
    switch (igo990) {
    case 0: goto L90;
    }
} /* lsei_ */

/* Subroutine */ 
int lsei::lsi_(double *w, integer *mdw, integer *ma, integer *mg, 
    integer *n, double *prgopt, double *x, double *rnorm, integer *mode, double *
    ws, integer *ip)
{
    /* Initialized data */

    double zero = 0.;
    double srelpr = 0.;
    double one = 1.;
    double half = .5;

    /* Format strings */
//    char fmt_40[] = "";
//    char fmt_60[] = "";

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3, i__4;
    double r__1, r__2;

    /* Builtin functions */
//    double sqrt(doublereal);

    /* Local variables */
    integer link;  
    integer last;
    integer next, igo990, igo994, i__, j, k, l, m;
    integer krank;
    double anorm;
    integer n1, n2;
    integer n3;
    double xnorm;
    integer ii;
    double rb;
    integer il, minman, mdlpdp, im1, ip1, np1;
    double fac, gam, tau;
    logical cov;
    integer key;
    double tol;
    integer map1, krm1, krp1;

    /* Assigned format variables */
//    char *igo994_fmt, *igo990_fmt;

    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;

    /* Function Body */

/*     COMPUTE MACHINE PRECISION=SRELPR ONLY WHEN NECESSARY. */
    if (! (M_precision == zero)) goto L30;
    M_precision = one;
L10:
    if (one + M_precision == one)   goto L20;
    M_precision *= half;
    goto L10;
L20:
    M_precision *= 10;
L30:
    srelpr = M_precision;

    *mode = 0;
    *rnorm = zero;
    m = *ma + *mg;
    np1 = *n + 1;
    krank = 0;
    if (*n <= 0 || m <= 0) {
    goto L70;
    }
    igo994 = 0;
//    igo994_fmt = fmt_40;
    goto L500;

/*     PROCESS-OPTION-VECTOR */

/*     COMPUTE MATRIX NORM OF LEAST SQUARES EQUAS. */
L40:
    anorm = zero;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
    r__1 = anorm, r__2 = sasum_(ma, &w[j * w_dim1 + 1], &c__1);
    anorm = dmax(r__1,r__2);
/* L50: */
    }

/*     SET TOL FOR HFTI( ) RANK TEST. */
    tau = tol * anorm;

/*     COMPUTE HOUSEHOLDER ORTHOGONAL DECOMP OF MATRIX. */
    if (*n > 0) {
    ws[1] = zero;
    }
    scopy_(n, &ws[1], &c__0, &ws[1], &c__1);
    scopy_(ma, &w[np1 * w_dim1 + 1], &c__1, &ws[1], &c__1);
    k = max(m,*n);
    minman = min(*ma,*n);
    n1 = k + 1;
    n2 = n1 + *n;
    hfti_(&w[w_offset], mdw, ma, n, &ws[1], &c__1, &c__1, &tau, &krank, rnorm,
         &ws[n2], &ws[n1], &ip[1]);
    fac = one;
    gam = (double) (*ma - krank);
    if (krank < *ma) {
/* Computing 2nd power */
    r__1 = *rnorm;
    fac = r__1 * r__1 / gam;
    }
    igo990 = 0;
//    igo990_fmt = fmt_60;
    goto L80;

/*     REDUCE-TO-LPDP-AND-SOLVE */
L60:
L70:
    ip[1] = krank;
    ip[2] = *n + max(m,*n) + (*mg + 2) * (*n + 7);
    return 0;
L80:

/*     TO REDUCE-TO-LPDP-AND-SOLVE */
    map1 = *ma + 1;

/*     COMPUTE INEQ. RT-HAND SIDE FOR LPDP. */
    if (! (*ma < m)) {
    goto L260;
    }
    if (! (minman > 0)) {
    goto L160;
    }
    i__1 = m;
    for (i__ = map1; i__ <= i__1; ++i__) {
    w[i__ + np1 * w_dim1] -= sdot_(n, &w[i__ + w_dim1], mdw, &ws[1], &c__1);
/* L90: */
    }
    i__1 = minman;
    for (i__ = 1; i__ <= i__1; ++i__) {
    j = ip[i__];

/*     APPLY PERMUTATIONS TO COLS OF INEQ. CONSTRAINT MATRIX. */
    sswap_(mg, &w[map1 + i__ * w_dim1], &c__1, &w[map1 + j * w_dim1], &c__1);
/* L100: */
    }

/*     APPLY HOUSEHOLDER TRANSFORMATIONS TO CONSTRAINT MATRIX. */
    if (! (0 < krank && krank < *n)) {
    goto L120;
    }
    i__1 = krank;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = krank + 1 - ii;
    l = n1 + i__;
    i__2 = krank + 1;
    h12_(&c__2, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[l - 1], &w[  map1 + w_dim1], mdw, &c__1, mg);
/* L110: */
    }

/*     COMPUTE PERMUTED INEQ. CONSTR. MATRIX TIMES R-INVERSE. */
L120:
    i__1 = m;
    for (i__ = map1; i__ <= i__1; ++i__) {
    if (! (0 < krank)) {
        goto L140;
    }
    i__2 = krank;
    for (j = 1; j <= i__2; ++j) {
        i__3 = j - 1;
        w[i__ + j * w_dim1] = (w[i__ + j * w_dim1] - sdot_(&i__3, &w[j * w_dim1 + 1], &c__1, 
            &w[i__ + w_dim1], mdw)) / w[j + j *  w_dim1];
/* L130: */
    }
L140:
/* L150: */
    ;
    }

/*     SOLVE THE REDUCED PROBLEM WITH LPDP ALGORITHM, */
/*     THE LEAST PROJECTED DISTANCE PROBLEM. */
L160:
    i__1 = *n - krank;
    lpdp_(&w[map1 + w_dim1], mdw, mg, &krank, &i__1, &prgopt[1], &x[1], &
        xnorm, &mdlpdp, &ws[n2], &ip[*n + 1]);
    if (! (mdlpdp == 1)) {
    goto L240;
    }
    if (! (krank > 0)) {
    goto L180;
    }

/*     COMPUTE SOLN IN ORIGINAL COORDINATES. */
    i__1 = krank;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = krank + 1 - ii;
    i__2 = ii - 1;
    x[i__] = (x[i__] - sdot_(&i__2, &w[i__ + (i__ + 1) * w_dim1], mdw, &x[i__ + 1], &c__1)) / w[i__ + i__ * w_dim1];
/* L170: */
    }

/*     APPLY HOUSEHOLDER TRANS. TO SOLN VECTOR. */
L180:
    if (! (0 < krank && krank < *n)) {
    goto L200;
    }
    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
    l = n1 + i__;
    i__2 = krank + 1;
    h12_(&c__2, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[l - 1], &x[1], &c__1, &c__1, &c__1);
/* L190: */
    }
L200:
    if (! (minman > 0)) {
    goto L230;
    }

/*     REPERMUTE VARIABLES TO THEIR INPUT ORDER. */
    i__1 = minman;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = minman + 1 - ii;
    j = ip[i__];
    sswap_(&c__1, &x[i__], &c__1, &x[j], &c__1);
/* L210: */
    }

/*     VARIABLES ARE NOW IN ORIG. COORDINATES. */
/*     ADD SOLN OF UNSCONSTRAINED PROB. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    x[i__] += ws[i__];
/* L220: */
    }

/*     COMPUTE THE RESIDUAL VECTOR NORM. */
/* Computing 2nd power */
    r__1 = *rnorm;
/* Computing 2nd power */
    r__2 = xnorm;
    *rnorm = sqrt(r__1 * r__1 + r__2 * r__2);
L230:
    goto L250;
L240:
    *mode = 2;
L250:
    goto L270;
L260:
    scopy_(n, &ws[1], &c__1, &x[1], &c__1);
L270:
    if (! (cov && krank > 0)) {
    goto L490;
    }

/*     COMPUTE COVARIANCE MATRIX BASED ON THE ORTHOGONAL DECOMP. */
/*     FROM HFTI( ). */

    krm1 = krank - 1;
    krp1 = krank + 1;

/*     COPY DIAG. TERMS TO WORKING ARRAY. */
    i__1 = *mdw + 1;
    scopy_(&krank, &w[w_offset], &i__1, &ws[n2], &c__1);

/*     RECIPROCATE DIAG. TERMS. */
    i__1 = krank;
    for (j = 1; j <= i__1; ++j) {
    w[j + j * w_dim1] = one / w[j + j * w_dim1];
/* L280: */
    }
    if (! (krank > 1)) {
    goto L310;
    }

/*     INVERT THE UPPER TRIANGULAR QR FACTOR ON ITSELF. */
    i__1 = krm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ip1 = i__ + 1;
    i__2 = krank;
    for (j = ip1; j <= i__2; ++j) {
        i__3 = j - i__;
        w[i__ + j * w_dim1] = -sdot_(&i__3, &w[i__ + i__ * w_dim1], mdw, &
            w[i__ + j * w_dim1], &c__1) * w[j + j * w_dim1];
/* L290: */
    }
/* L300: */
    }

/*     COMPUTE THE INVERTED FACTOR TIMES ITS TRANSPOSE. */
L310:
    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
    i__2 = krank;
    for (j = i__; j <= i__2; ++j) {
        i__3 = krank + 1 - j;
        w[i__ + j * w_dim1] = sdot_(&i__3, &w[i__ + j * w_dim1], mdw, &w[j + j * w_dim1], mdw);
/* L320: */
    }
/* L330: */
    }
    if (! (krank < *n)) {
    goto L450;
    }

/*     ZERO OUT LOWER TRAPEZOIDAL PART. */
/*     COPY UPPER TRI. TO LOWER TRI. PART. */
    i__1 = krank;
    for (j = 1; j <= i__1; ++j) {
    scopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L340: */
    }
    i__1 = *n;
    for (i__ = krp1; i__ <= i__1; ++i__) {
    w[i__ + w_dim1] = zero;
    scopy_(&i__, &w[i__ + w_dim1], &c__0, &w[i__ + w_dim1], mdw);
/* L350: */
    }

/*     APPLY RIGHT SIDE TRANSFORMATIONS TO LOWER TRI. */
    n3 = n2 + krp1;
    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
    l = n1 + i__;
    k = n2 + i__;
    rb = ws[l - 1] * ws[k - 1];
    if (! (rb < zero)) {
        goto L420;
    }

/*     IF RB.GE.ZERO, TRANSFORMATION CAN BE REGARDED AS ZERO. */
    rb = one / rb;

/*     STORE UNSCALED RANK-ONE HOUSEHOLDER UPDATE IN WORK ARRAY. */
    ws[n3] = zero;
    scopy_(n, &ws[n3], &c__0, &ws[n3], &c__1);
    l = n1 + i__;
    k = n3 + i__;
    ws[k - 1] = ws[l - 1];
    i__2 = *n;
    for (j = krp1; j <= i__2; ++j) {
        k = n3 + j;
        ws[k - 1] = w[i__ + j * w_dim1];
/* L360: */
    }
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
        l = n3 + i__;
        k = n3 + j;
        i__3 = j - i__;
        i__4 = *n - j + 1;
        ws[j] = sdot_(&i__3, &w[j + i__ * w_dim1], mdw, &ws[l - 1], &c__1)
             + sdot_(&i__4, &w[j + j * w_dim1], &c__1, &ws[k - 1], &c__1);
        ws[j] *= rb;
/* L370: */
    }
    l = n3 + i__;
    i__2 = *n - i__ + 1;
    gam = sdot_(&i__2, &ws[l - 1], &c__1, &ws[i__], &c__1) * rb;
    gam *= half;
    i__2 = *n - i__ + 1;
    saxpy_(&i__2, &gam, &ws[l - 1], &c__1, &ws[i__], &c__1);
    i__2 = *n;
    for (j = i__; j <= i__2; ++j) {
        if (! (i__ > 1)) {
        goto L390;
        }
        im1 = i__ - 1;
        k = n3 + j;
        i__3 = im1;
        for (l = 1; l <= i__3; ++l) {
        w[j + l * w_dim1] += ws[k - 1] * ws[l];
/* L380: */
        }
L390:
        k = n3 + j;
        i__3 = j;
        for (l = i__; l <= i__3; ++l) {
        il = n3 + l;
        w[j + l * w_dim1] = w[j + l * w_dim1] + ws[j] * ws[il - 1] +    ws[l] * ws[k - 1];
/* L400: */
        }
/* L410: */
    }
L420:
/* L430: */
    ;
    }

/*     COPY LOWER TRI. TO UPPER TRI. TO SYMMETRIZE THE COVARIANCE MATRIX. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    scopy_(&i__, &w[i__ + w_dim1], mdw, &w[i__ * w_dim1 + 1], &c__1);
/* L440: */
    }

/*     REPERMUTE ROWS AND COLS. */
L450:
    i__1 = minman;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = minman + 1 - ii;
    k = ip[i__];
    if (! (i__ != k)) {
        goto L460;
    }
    sswap_(&c__1, &w[i__ + i__ * w_dim1], &c__1, &w[k + k * w_dim1], &c__1);
    i__2 = i__ - 1;
    sswap_(&i__2, &w[i__ * w_dim1 + 1], &c__1, &w[k * w_dim1 + 1], &c__1);
    i__2 = k - i__ - 1;
    sswap_(&i__2, &w[i__ + (i__ + 1) * w_dim1], mdw, &w[i__ + 1 + k * w_dim1], &c__1);
    i__2 = *n - k;
    sswap_(&i__2, &w[i__ + (k + 1) * w_dim1], mdw, &w[k + (k + 1) * w_dim1], mdw);
L460:
/* L470: */
    ;
    }

/*     PUT IN NORMALIZED RESIDUAL SUM OF SQUARES SCALE FACTOR */
/*     AND SYMMETRIZE THE RESULTING COVARIANCE MARIX. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    sscal_(&j, &fac, &w[j * w_dim1 + 1], &c__1);
    scopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L480: */
    }
L490:
    goto L540;
L500:

/*     TO PROCESS-OPTION-VECTOR */

/*     THE NOMINAL TOLERANCE USED IN THE CODE, */
    tol = sqrt(srelpr);
    cov = FALSE_;
    last = 1;
    link = prgopt[1];
L510:
    if (! (link > 1)) {
    goto L520;
    }
    key = prgopt[last + 1];
    if (key == 1) {
    cov = prgopt[last + 2] != zero;
    }
    if (key == 5) {
/* Computing MAX */
    r__1 = srelpr, r__2 = prgopt[last + 2];
    tol = dmax(r__1,r__2);
    }
    next = prgopt[link];
    last = link;
    link = next;
    goto L510;
L520:
    goto L530;
L530:
    switch (igo994) {
    case 0: goto L40;
    }
L540:
    switch (igo990) {
    case 0: goto L60;
    }
} /* lsi_ */

/* Subroutine */ int lsei::lpdp_(double *a, integer *mda, integer *m, integer *n1, 
    integer *n2, double *prgopt, double *x, double *wnorm, integer *mode, double *ws, integer *is)
{
    /* Initialized data */

    double zero = 0.;
    double one = 1.;
    double fac = .1;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer i__, j, l, n;
    integer modew;
    double rnorm;
    double ynorm, sc;
    integer iw, ix, np1;

    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --prgopt;
    --x;
    --ws;
    --is;

    /* Function Body */
    n = *n1 + *n2;
    *mode = 1;
    if (! (*m <= 0)) {
    goto L20;
    }
    if (! (n > 0)) {
    goto L10;
    }
    x[1] = zero;
    scopy_(&n, &x[1], &c__0, &x[1], &c__1);
L10:
    *wnorm = zero;
    return 0;
L20:
    np1 = n + 1;

/*     SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    sc = snrm2_(&n, &a[i__ + a_dim1], mda);
    if (! (sc != zero)) {
        goto L30;
    }
    sc = one / sc;
    sscal_(&np1, &sc, &a[i__ + a_dim1], mda);
L30:
/* L40: */
    ;
    }

/*     SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO). */
    ynorm = snrm2_(m, &a[np1 * a_dim1 + 1], &c__1);
    if (! (ynorm != zero)) {
    goto L50;
    }
    sc = one / ynorm;
    sscal_(m, &sc, &a[np1 * a_dim1 + 1], &c__1);

/*     SCALE COLS OF MATRIX H. */
L50:
    j = *n1 + 1;
L60:
    if (! (j <= n)) {
    goto L70;
    }
    sc = snrm2_(m, &a[j * a_dim1 + 1], &c__1);
    if (sc != zero) {
    sc = one / sc;
    }
    sscal_(m, &sc, &a[j * a_dim1 + 1], &c__1);
    x[j] = sc;
    ++j;
    goto L60;
L70:
    if (! (*n1 > 0)) {
    goto L130;
    }

/*     COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY. */
    scopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
    iw += *n2;

/*     MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY. */
    scopy_(n1, &a[i__ + a_dim1], mda, &ws[iw + 1], &c__1);
    iw += *n1;

/*     MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY. */
    ws[iw + 1] = a[i__ + np1 * a_dim1];
    ++iw;
/* L80: */
    }
    ws[iw + 1] = zero;
    scopy_(&n, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += n;
    ws[iw + 1] = one;
    ++iw;

/*     SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE */
/*     MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR */
/*     F = TRANSPOSE OF (0,...,0,1). */
    ix = iw + 1;
    iw += *m;

/*     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = np1 - *n2;
    wnnls_(&ws[1], &np1, n2, &i__1, m, &c__0, &prgopt[1], &ws[ix], &rnorm, &modew, &is[1], &ws[iw + 1]);

/*     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W. */
    sc = one - sdot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (! (one + fac * dabs(sc) != one && rnorm > zero)) {
    goto L110;
    }
    sc = one / sc;
    i__1 = *n1;
    for (j = 1; j <= i__1; ++j) {
    x[j] = sc * sdot_(m, &a[j * a_dim1 + 1], &c__1, &ws[ix], &c__1);
/* L90: */
    }

/*     COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS VECTOR. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    a[i__ + np1 * a_dim1] -= sdot_(n1, &a[i__ + a_dim1], mda, &x[1], &c__1);
/* L100: */
    }
    goto L120;
L110:
    *mode = 2;
    return 0;
L120:
L130:
    if (! (*n2 > 0)) {
    goto L180;
    }

/*     COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    scopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
    iw += *n2;
    ws[iw + 1] = a[i__ + np1 * a_dim1];
    ++iw;
/* L140: */
    }
    ws[iw + 1] = zero;
    scopy_(n2, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += *n2;
    ws[iw + 1] = one;
    ++iw;
    ix = iw + 1;
    iw += *m;

/*     SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE */
/*     OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE */
/*     OF (0,...,0,1)). */

/*     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = *n2 + 1;
    i__2 = *n2 + 1;
    wnnls_(&ws[1], &i__1, &c__0, &i__2, m, &c__0, &prgopt[1], &ws[ix], &rnorm,&modew, &is[1], &ws[iw + 1]);

/*     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z. */
    sc = one - sdot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (! (one + fac * dabs(sc) != one && rnorm > zero)) {
    goto L160;
    }
    sc = one / sc;
    i__1 = *n2;
    for (j = 1; j <= i__1; ++j) {
    l = *n1 + j;
    x[l] = sc * sdot_(m, &a[l * a_dim1 + 1], &c__1, &ws[ix], &c__1) * x[l];
/* L150: */
    }
    goto L170;
L160:
    *mode = 2;
    return 0;
L170:

/*     ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION. */
L180:
    sscal_(&n, &ynorm, &x[1], &c__1);
    *wnorm = snrm2_(n1, &x[1], &c__1);
    return 0;
} /* lpdp_ */

/* Subroutine */ int lsei::wnnls_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *n, integer *l, double *prgopt, double *x, double *rnorm, integer *mode)
{
    /****************** seting dimentions of working arrays************************/
    
    int k_temp =  *ma+*me +5*(*n) +1;       // max(*ma+*mg, *n));
    double *work = new double[k_temp/*2*(*me+*n)+ k_temp + (*mg + 2)*(*n + 7)*/]; 
    integer *iwork = new integer[*ma +*me + *n +1/*+2*(*n)+2*/];

    /*****************************************************************************/
    int val = wnnls_(w, mdw, me, ma, n, l, prgopt, x, rnorm, mode, iwork, work);
    delete[] work;
    delete[] iwork;
    return val;
}

/* Subroutine */ int lsei::wnnls_(double *w, integer *mdw, integer *me, integer *ma, 
    integer *n, integer *l, double *prgopt, double *x, double *rnorm, integer *
    mode, integer *iwork, double *work)
{
    /* System generated locals */
    integer w_dim1, w_offset;

    /* Local variables */
    integer nerr, iopt;
    double dummy;
    integer l1, l2, l3, l4, l5;
    integer lw;
    integer liw;

    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --prgopt;
    --x;
    --iwork;
    --work;
    /* Function Body */
    *mode = 0;
    if (*ma + *me <= 0 || *n <= 0) {return 0;}
    goto L30;
    if (! (iwork[1] > 0)) {
    goto L20;
    }
    lw = *me + *ma + *n * 5;
    if (! (iwork[1] < lw)) {
    goto L10;
    }
    nerr = 2;
    iopt = 1;
    xerrwv_("WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR WORK(*), NEED LW=I1 BELOW", &c__70, 
        &nerr, &iopt, &c__1, &lw, &c__0, &c__0, &dummy, &dummy/*, (ftnlen)70*/);
    *mode = 2;
    return 0;
L10:
L20:
    if (! (iwork[2] > 0)) {
    goto L40;
    }
    liw = *me + *ma + *n;
    if (! (iwork[2] < liw)) {
    goto L30;
    }
    nerr = 2;
    iopt = 1;
    xerrwv_("WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR IWORK(*), NEED LIW=I1 BELOW", &c__72, 
        &nerr, &iopt, &c__1, &liw, &c__0, &c__0, &dummy, &dummy/*, ( ftnlen)72*/);
    *mode = 2;
    return 0;
L30:
L40:
    if (! (*mdw < *me + *ma)) {
    goto L50;
    }
    nerr = 1;
    iopt = 1;
    xerror_("WNNLS( ), THE VALUE MDW.LT.ME+MA IS AN ERROR", &c__44, &nerr, &iopt/*, (ftnlen)44*/);
    *mode = 2;
    return 0;
L50:
    if (0 <= *l && *l <= *n) {
    goto L60;
    }
    nerr = 2;
    iopt = 1;
    xerror_("WNNLS( ), L.LE.0.AND.L.LE.N IS REQUIRED", &c__39, &nerr, &iopt/*, (        ftnlen)39*/);
    *mode = 2;
    return 0;

/*     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS */
/*     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS */
/*     REQUIRED BY THE MAIN SUBROUTINE WNLSM( ). */

L60:
    l1 = *n + 1;
    l2 = l1 + *n;
    l3 = l2 + *me + *ma;
    l4 = l3 + *n;
    l5 = l4 + *n;

    wnlsm_(&w[w_offset], mdw, me, ma, n, l, &prgopt[1], &x[1], rnorm, mode, &
        iwork[1], &iwork[l1], &work[1], &work[l1], &work[l2], &work[l3], &work[l4], &work[l5]);

    return 0;
} /* wnnls_ */

/* Subroutine */ 
int lsei::wnlsm_(double *w, integer *mdw, integer *mme, integer *ma, 
    integer *n, integer *l, double *prgopt, double *x, double *rnorm, integer *
    mode, integer *ipivot, integer *itype, double *wd, double *h__, double *
    scale, double *z__, double *temp, double *d__)
{
    /* Initialized data */

    double zero = 0.;
    double one = 1.;
    double two = 2.;
    double srelpr = 0.;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    double r__1, r__2;

    /* Builtin functions */
   // double sqrt(doublereal);

    /* Local variables */
    logical done;
    double amax, dope[4];
    integer jcon, link, imax;
    double alsq;
    integer last, iter, nerr, isol, iopt;
    double wmax;
    integer next, nopt, igo980, igo991, igo983, igo938, igo958, igo995,
         igo986, igo977, igo998, igo897;
    integer i__, j, m;
    double t, alpha;
    integer idope[8];    
    integer krank, nlink;
    double bnorm;
    integer itemp, itmax, iwmax;
    integer nsoln;
    integer l1;
    double z2;    
    double alamda;
    integer me, jj, jp;
    logical feasbl;
    double sm, eanorm;    
    double sparam[5], zz;
    logical hitcon;
    integer ntimes;
    double blowup;
    integer jm1, nm1, lp1;    
    integer np1;    
    double fac;
    integer key;
    double tau;
    integer niv;
    logical pos= FALSE_;
    integer mep1, krp1, niv1, nsp1;

    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --prgopt;
    --x;
    --ipivot;
    --itype;
    --wd;
    --h__;
    --scale;
    --z__;
    --temp;
    --d__;

    /* Function Body */

/*     INITIALIZE-VARIABLES */
    igo998 = 0;
//    igo998_fmt = fmt_10;
    goto L180;

/*     PERFORM INITIAL TRIANGULARIZATION IN THE SUBMATRIX */
/*     CORRESPONDING TO THE UNCONSTRAINED VARIABLES USING */
/*     THE PROCEDURE INITIALLY-TRIANGULARIZE. */
L10:
    igo995 = 0;
//    igo995_fmt = fmt_20;
    goto L280;

/*     PERFORM WNNLS ALGORITHM USING THE FOLLOWING STEPS. */

/*     UNTIL(DONE) */

/*        COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT */

/*        WHEN (HITCON) ADD-CONSTRAINTS */

/*        ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT */

/*        FIN */

/*     COMPUTE-FINAL-SOLUTION */

L20:
    if (done) {
    goto L80;
    }

    igo991 = 0;
//    igo991_fmt = fmt_30;
    goto L300;

/*     COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT */

L30:
    if (! hitcon) {
    goto L50;
    }
    igo986 = 0;
//    igo986_fmt = fmt_40;
    goto L370;
L40:
    goto L70;

/*     WHEN (HITCON) ADD-CONSTRAINTS */

L50:
    igo983 = 0;
//    igo983_fmt = fmt_60;
    goto L640;
L60:

/*     ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT */

L70:
    goto L20;

L80:
    igo980 = 0;
//    igo980_fmt = fmt_90;
    goto L1000;

/*     COMPUTE-FINAL-SOLUTION */

L90:
    return 0;
L100:

/*     TO PROCESS-OPTION-VECTOR */
    fac = 1e-4;

/*     THE NOMINAL TOLERANCE USED IN THE CODE, */
    tau = sqrt(srelpr);

/*     THE NOMINAL BLOW-UP FACTOR USED IN THE CODE. */
    blowup = tau;

/*     THE NOMINAL COLUMN SCALING USED IN THE CODE IS */
/*     THE IDENTITY SCALING. */
    d__[1] = one;
    scopy_(n, &d__[1], &c__0, &d__[1], &c__1);

/*     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE. */
    nopt = 1000;

/*     DEFINE BOUND FOR POSITIVE VALUE OF LINK. */
    nlink = 100000;
    ntimes = 0;
    last = 1;
    link = prgopt[1];
    if (! (link <= 0 || link > nlink)) {
    goto L110;
    }
    nerr = 3;
    iopt = 1;
    xerror_("WNNLS( ) THE OPTION VECTOR IS UNDEFINED", &c__39, &nerr, &iopt/*, (    ftnlen)39*/);
    *mode = 2;
    return 0;
L110:
    if (! (link > 1)) {
    goto L160;
    }
    ++ntimes;
    if (! (ntimes > nopt)) {
    goto L120;
    }
    nerr = 3;
    iopt = 1;
    xerror_("WNNLS( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING.", &c__53, &nerr, &iopt/*, (ftnlen)53*/);
    *mode = 2;
    return 0;
L120:
    key = prgopt[last + 1];
    if (! (key == 6 && prgopt[last + 2] != zero)) {
    goto L140;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    t = snrm2_(&m, &w[j * w_dim1 + 1], &c__1);
    if (t != zero) {
        t = one / t;
    }
    d__[j] = t;
/* L130: */
    }
L140:
    if (key == 7) {
    scopy_(n, &prgopt[last + 2], &c__1, &d__[1], &c__1);
    }
    if (key == 8) {
/* Computing MAX */
    r__1 = srelpr, r__2 = prgopt[last + 2];
    tau = dmax(r__1,r__2);
    }
    if (key == 9) {
/* Computing MAX */
    r__1 = srelpr, r__2 = prgopt[last + 2];
    blowup = dmax(r__1,r__2);
    }
    next = prgopt[link];
    if (! (next <= 0 || next > nlink)) {
    goto L150;
    }
    nerr = 3;
    iopt = 1;
    xerror_("WNNLS( ) THE OPTION VECTOR IS UNDEFINED", &c__39, &nerr, &iopt/*, ( ftnlen)39*/);
    *mode = 2;
    return 0;
L150:
    last = link;
    link = next;
    goto L110;
L160:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    sscal_(&m, &d__[j], &w[j * w_dim1 + 1], &c__1);
/* L170: */
    }
    goto L1260;
L180:

/*     TO INITIALIZE-VARIABLES */

/*     SRELPR IS THE PRECISION FOR THE PARTICULAR MACHINE */
/*     BEING USED.  THIS LOGIC AVOIDS RECOMPUTING IT EVERY ENTRY. */
    if (! (M_precision == zero))    goto L210;
    M_precision = one;
L190:
    if (one + M_precision == one)   goto L200;
    M_precision /= two;
    goto L190;
L200:
    M_precision *=10;
L210:
    srelpr = M_precision;

    m = *ma + *mme;
    me = *mme;
    mep1 = me + 1;
    igo977 = 0;
//    igo977_fmt = fmt_220;
    goto L100;

/*     PROCESS-OPTION-VECTOR */
L220:
    done = FALSE_;
    iter = 0;
    itmax = (*n - *l) * 10;//gleb, was 3
    *mode = 0;
    lp1 = *l + 1;
    nsoln = *l;
    nsp1 = nsoln + 1;
    np1 = *n + 1;
    nm1 = *n - 1;
    l1 = min(m,*l);

/*     COMPUTE SCALE FACTOR TO APPLY TO EQUAL. CONSTRAINT EQUAS. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    wd[j] = sasum_(&m, &w[j * w_dim1 + 1], &c__1);
/* L230: */
    }
    imax = isamax_(n, &wd[1], &c__1);
    eanorm = wd[imax];
    bnorm = sasum_(&m, &w[np1 * w_dim1 + 1], &c__1);
    alamda = eanorm / (srelpr * fac);

/*     DEFINE SCALING DIAG MATRIX FOR MOD GIVENS USAGE AND */
/*     CLASSIFY EQUATION TYPES. */
/* Computing 2nd power */
    r__1 = alamda;
    alsq = r__1 * r__1;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     WHEN EQU I IS HEAVILY WEIGHTED ITYPE(I)=0, ELSE ITYPE(I)=1. */
    if (! (i__ <= me)) {
        goto L240;
    }
    t = alsq;
    itemp = 0;
    goto L250;
L240:
    t = one;
    itemp = 1;
L250:
    scale[i__] = t;
    itype[i__] = itemp;
/* L260: */
    }

/*     SET THE SOLN VECTOR X(*) TO ZERO AND THE COL INTERCHANGE */
/*     MATRIX TO THE IDENTITY. */
    x[1] = zero;
    scopy_(n, &x[1], &c__0, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ipivot[i__] = i__;
/* L270: */
    }
    goto L1230;
L280:

/*     TO INITIALLY-TRIANGULARIZE */

/*     SET FIRST L COMPS. OF DUAL VECTOR TO ZERO BECAUSE */
/*     THESE CORRESPOND TO THE UNCONSTRAINED VARIABLES. */
    if (! (*l > 0)) {
    goto L290;
    }
    wd[1] = zero;
    scopy_(l, &wd[1], &c__0, &wd[1], &c__1);

/*     THE ARRAYS IDOPE(*) AND DOPE(*) ARE USED TO PASS */
/*     INFORMATION TO WNLIT().  THIS WAS DONE TO AVOID */
/*     A LONG CALLING SEQUENCE OR THE USE OF COMMON. */
L290:
    idope[0] = me;
    idope[1] = mep1;
    idope[2] = 0;
    idope[3] = 1;
    idope[4] = nsoln;
    idope[5] = 0;
    idope[6] = 1;
    idope[7] = l1;

    dope[0] = alsq;
    dope[1] = eanorm;
    dope[2] = fac;
    dope[3] = tau;
    wnlit_(&w[w_offset], mdw, &m, n, l, &ipivot[1], &itype[1], &h__[1], &
        scale[1], rnorm, idope, dope, &done);
    me = idope[0];
    mep1 = idope[1];
    krank = idope[2];
    krp1 = idope[3];
    nsoln = idope[4];
    niv = idope[5];
    niv1 = idope[6];
    l1 = idope[7];
    goto L1240;
L300:

/*     TO COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT */

/*     SOLVE THE TRIANGULAR SYSTEM OF CURRENTLY NON-ACTIVE */
/*     VARIABLES AND STORE THE SOLUTION IN Z(*). */

/*     SOLVE-SYSTEM */
    igo958 = 0;
//    igo958_fmt = fmt_310;
    goto L1110;

/*     INCREMENT ITERATION COUNTER AND CHECK AGAINST MAX. NUMBER */
/*     OF ITERATIONS. */
L310:
    ++iter;
    if (! (iter > itmax)) {
    goto L320;
    }
    *mode = 1;
    done = TRUE_;

/*     CHECK TO SEE IF ANY CONSTRAINTS HAVE BECOME ACTIVE. */
/*     IF SO, CALCULATE AN INTERPOLATION FACTOR SO THAT ALL */
/*     ACTIVE CONSTRAINTS ARE REMOVED FROM THE BASIS. */
L320:
    alpha = two;
    hitcon = FALSE_;
    if (! (*l < nsoln)) {
    goto L360;
    }
    i__1 = nsoln;
    for (j = lp1; j <= i__1; ++j) {
    zz = z__[j];
    if (! (zz <= zero)) {
        goto L340;
    }
    t = x[j] / (x[j] - zz);
    if (! (t < alpha)) {
        goto L330;
    }
    alpha = t;
    jcon = j;
L330:
    hitcon = TRUE_;
L340:
/* L350: */
    ;
    }
L360:
    goto L1220;
L370:

/*     TO ADD-CONSTRAINTS */

/*     USE COMPUTED ALPHA TO INTERPOLATE BETWEEN LAST */
/*     FEASIBLE SOLUTION X(*) AND CURRENT UNCONSTRAINED */
/*     (AND INFEASIBLE) SOLUTION Z(*). */
    if (! (lp1 <= nsoln)) {
    goto L390;
    }
    i__1 = nsoln;
    for (j = lp1; j <= i__1; ++j) {
    x[j] += alpha * (z__[j] - x[j]);
/* L380: */
    }
L390:
    feasbl = FALSE_;
    goto L410;
L400:
    if (feasbl) {
    goto L610;
    }

/*     REMOVE COL JCON AND SHIFT COLS JCON+1 THROUGH N TO THE */
/*     LEFT. SWAP COL JCON INTO THE N-TH POSITION.  THIS ACHIEVES */
/*     UPPER HESSENBERG FORM FOR THE NONACTIVE CONSTRAINTS AND */
/*     LEAVES AN UPPER HESSENBERG MATRIX TO RETRIANGULARIZE. */
L410:
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    t = w[i__ + jcon * w_dim1];
    i__2 = *n - jcon;
    scopy_(&i__2, &w[i__ + (jcon + 1) * w_dim1], mdw, &w[i__ + jcon *   w_dim1], mdw);
    w[i__ + *n * w_dim1] = t;
/* L420: */
    }

/*     UPDATE PERMUTED INDEX VECTOR TO REFLECT THIS SHIFT AND SWAP. */
    itemp = ipivot[jcon];
    if (! (jcon < *n)) {
    goto L440;
    }
    i__1 = nm1;
    for (i__ = jcon; i__ <= i__1; ++i__) {
    ipivot[i__] = ipivot[i__ + 1];
/* L430: */
    }
L440:
    ipivot[*n] = itemp;

/*     SIMILARLY REPERMUTE X(*) VECTOR. */
    i__1 = *n - jcon;
    scopy_(&i__1, &x[jcon + 1], &c__1, &x[jcon], &c__1);
    x[*n] = zero;
    nsp1 = nsoln;
    --nsoln;
    niv1 = niv;
    --niv;

/*     RETRIANGULARIZE UPPER HESSENBERG MATRIX AFTER ADDING CONSTRAINTS. */
    j = jcon;
    i__ = krank + jcon - *l;
L450:
    if (! (j <= nsoln)) {
    goto L570;
    }
    if (! (itype[i__] == 0 && itype[i__ + 1] == 0)) {
    goto L470;
    }
    igo938 = 0;
//    igo938_fmt = fmt_460;
    goto L620;

/*     (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) ZERO-IP1-TO-I-IN-COL-J */
L460:
    goto L560;
L470:
    if (! (itype[i__] == 1 && itype[i__ + 1] == 1)) {
    goto L490;
    }
    igo938 = 1;
//    igo938_fmt = fmt_480;
    goto L620;

/*     (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) ZERO-IP1-TO-I-IN-COL-J */
L480:
    goto L560;
L490:
    if (! (itype[i__] == 1 && itype[i__ + 1] == 0)) {
    goto L510;
    }
    sswap_(&np1, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1], mdw);
    sswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
    itemp = itype[i__ + 1];
    itype[i__ + 1] = itype[i__];
    itype[i__] = itemp;

/*     SWAPPED ROW WAS FORMERLY A PIVOT ELT., SO IT WILL */
/*     BE LARGE ENOUGH TO PERFORM ELIM. */
    igo938 = 2;
//    igo938_fmt = fmt_500;
    goto L620;

/*     ZERO-IP1-TO-I-IN-COL-J */
L500:
    goto L560;
L510:
    if (! (itype[i__] == 0 && itype[i__ + 1] == 1)) {
    goto L550;
    }
/* Computing 2nd power */
    r__1 = w[i__ + j * w_dim1];
    t = scale[i__] * (r__1 * r__1) / alsq;
/* Computing 2nd power */
    r__1 = tau;
/* Computing 2nd power */
    r__2 = eanorm;
    if (! (t > r__1 * r__1 * (r__2 * r__2))) {
    goto L530;
    }
    igo938 = 3;
//    igo938_fmt = fmt_520;
    goto L620;
L520:
    goto L540;
L530:
    sswap_(&np1, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1], mdw);
    sswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
    itemp = itype[i__ + 1];
    itype[i__ + 1] = itype[i__];
    itype[i__] = itemp;
    w[i__ + 1 + j * w_dim1] = zero;
L540:
L550:
L560:
    ++i__;
    ++j;
    goto L450;

/*     SEE IF THE REMAINING COEFFS IN THE SOLN SET ARE FEASIBLE.  THEY */
/*     SHOULD BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.  IF ANY ARE */
/*     INFEASIBLE IT IS DUE TO ROUNDOFF ERROR.  ANY THAT ARE NON- */
/*     POSITIVE WILL BE SET TO ZERO AND REMOVED FROM THE SOLN SET. */
L570:
    if (! (lp1 <= nsoln)) {
    goto L590;
    }
    i__1 = nsoln;
    for (jcon = lp1; jcon <= i__1; ++jcon) {
    if (x[jcon] <= zero) {
        goto L600;
    }
/* L580: */
    }
L590:
    feasbl = TRUE_;
L600:
    goto L400;
L610:
    goto L1200;
L620:

/*     TO ZERO-IP1-TO-I-IN-COL-J */
    if (! (w[i__ + 1 + j * w_dim1] != zero)) {
    goto L630;
    }
    srotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1], &w[i__ + 1 +  j * w_dim1], sparam);
    w[i__ + 1 + j * w_dim1] = zero;
    i__1 = np1 - j;
    srotm_(&i__1, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 + (j + 1) *  w_dim1], mdw, sparam);
L630:
    goto L1290;
L640:

/*     TO PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT */
    scopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);
    if (! (nsoln < *n)) {
    goto L650;
    }
    x[nsp1] = zero;
    i__1 = *n - nsoln;
    scopy_(&i__1, &x[nsp1], &c__0, &x[nsp1], &c__1);
L650:
    i__ = niv1;
L660:
    if (! (i__ <= me)) {
    goto L690;
    }

/*     RECLASSIFY LEAST SQUARES EQATIONS AS EQUALITIES AS */
/*     NECESSARY. */
    if (! (itype[i__] == 0)) {
    goto L670;
    }
    ++i__;
    goto L680;
L670:
    sswap_(&np1, &w[i__ + w_dim1], mdw, &w[me + w_dim1], mdw);
    sswap_(&c__1, &scale[i__], &c__1, &scale[me], &c__1);
    itemp = itype[i__];
    itype[i__] = itype[me];
    itype[me] = itemp;
    mep1 = me;
    --me;
L680:
    goto L660;

/*     FORM INNER PRODUCT VECTOR WD(*) OF DUAL COEFFS. */
L690:
    if (! (nsp1 <= *n)) {
    goto L730;
    }
    i__1 = *n;
    for (j = nsp1; j <= i__1; ++j) {
    sm = zero;
    if (! (nsoln < m)) {
        goto L710;
    }
    i__2 = m;
    for (i__ = nsp1; i__ <= i__2; ++i__) {
        sm += scale[i__] * w[i__ + j * w_dim1] * w[i__ + np1 * w_dim1];
/* L700: */
    }
L710:
    wd[j] = sm;
/* L720: */
    }
L730:
    goto L750;
L740:
    if (pos || done) {
    goto L970;
    }

/*     FIND J SUCH THAT WD(J)=WMAX IS MAXIMUM.  THIS DETERMINES */
/*     THAT THE INCOMING COL J WILL REDUCE THE RESIDUAL VECTOR */
/*     AND BE POSITIVE. */
L750:
    wmax = zero;
    iwmax = nsp1;
    if (! (nsp1 <= *n)) {
    goto L780;
    }
    i__1 = *n;
    for (j = nsp1; j <= i__1; ++j) {
    if (! (wd[j] > wmax)) {
        goto L760;
    }
    wmax = wd[j];
    iwmax = j;
L760:
/* L770: */
    ;
    }
L780:
    if (! (wmax <= zero)) {
    goto L790;
    }
    done = TRUE_;
    goto L960;

/*     SET DUAL COEFF TO ZERO FOR INCOMING COL. */
L790:
    wd[iwmax] = zero;

/*     WMAX .GT. ZERO, SO OKAY TO MOVE COL IWMAX TO SOLN SET. */
/*     PERFORM TRANSFORMATION TO RETRIANGULARIZE, AND TEST */
/*     FOR NEAR LINEAR DEPENDENCE. */
/*     SWAP COL IWMAX INTO NSOLN-TH POSITION TO MAINTAIN UPPER */
/*     HESSENBERG FORM OF ADJACENT COLS, AND ADD NEW COL TO */
/*     TRIANGULAR DECOMPOSITION. */
    nsoln = nsp1;
    nsp1 = nsoln + 1;
    niv = niv1;
    niv1 = niv + 1;
    if (! (nsoln != iwmax)) {
    goto L800;
    }
    sswap_(&m, &w[nsoln * w_dim1 + 1], &c__1, &w[iwmax * w_dim1 + 1], &c__1);
    wd[iwmax] = wd[nsoln];
    wd[nsoln] = zero;
    itemp = ipivot[nsoln];
    ipivot[nsoln] = ipivot[iwmax];
    ipivot[iwmax] = itemp;

/*     REDUCE COL NSOLN SO THAT THE MATRIX OF NONACTIVE */
/*     CONSTRAINTS VARIABLES IS TRIANGULAR. */
L800:
    j = m;
L810:
    if (! (j > niv)) {
    goto L870;
    }
    jm1 = j - 1;
    jp = jm1;

/*     WHEN OPERATING NEAR THE ME LINE, TEST TO SEE IF THE PIVOT ELT. */
/*     IS NEAR ZERO.  IF SO, USE THE LARGEST ELT. ABOVE IT AS THE PIVOT. */
/*     THIS IS TO MAINTAIN THE SHARP INTERFACE BETWEEN WEIGHTED AND */
/*     NON-WEIGHTED ROWS IN ALL CASES. */
    if (! (j == mep1)) {
    goto L850;
    }
    imax = me;
/* Computing 2nd power */
    r__1 = w[me + nsoln * w_dim1];
    amax = scale[me] * (r__1 * r__1);
L820:
    if (! (jp >= niv)) {
    goto L840;
    }
/* Computing 2nd power */
    r__1 = w[jp + nsoln * w_dim1];
    t = scale[jp] * (r__1 * r__1);
    if (! (t > amax)) {
    goto L830;
    }
    imax = jp;
    amax = t;
L830:
    --jp;
    goto L820;
L840:
    jp = imax;
L850:
    if (! (w[j + nsoln * w_dim1] != zero)) {
    goto L860;
    }
    srotmg_(&scale[jp], &scale[j], &w[jp + nsoln * w_dim1], &w[j + nsoln *  w_dim1], sparam);
    w[j + nsoln * w_dim1] = zero;
    i__1 = np1 - nsoln;
    srotm_(&i__1, &w[jp + nsp1 * w_dim1], mdw, &w[j + nsp1 * w_dim1], mdw, sparam);
L860:
    j = jm1;
    goto L810;

/*     SOLVE FOR Z(NSOLN)=PROPOSED NEW VALUE FOR X(NSOLN). */
/*     TEST IF THIS IS NONPOSITIVE OR TOO LARGE. */
/*     IF THIS WAS TRUE OR IF THE PIVOT TERM WAS ZERO REJECT */
/*     THE COL AS DEPENDENT. */
L870:
    if (! (w[niv + nsoln * w_dim1] != zero)) {
    goto L890;
    }
    isol = niv;
    igo897 = 0;
//    igo897_fmt = fmt_880;
    goto L980;

/*     TEST-PROPOSED-NEW-COMPONENT */
L880:
    goto L940;
L890:
    if (! (niv <= me && w[mep1 + nsoln * w_dim1] != zero)) {
    goto L920;
    }

/*     TRY TO ADD ROW MEP1 AS AN ADDITIONAL EQUALITY CONSTRAINT. */
/*     CHECK SIZE OF PROPOSED NEW SOLN COMPONENT. */
/*     REJECT IT IF IT IS TOO LARGE. */
    isol = mep1;
    igo897 = 1;
//    igo897_fmt = fmt_900;
    goto L980;

/*     TEST-PROPOSED-NEW-COMPONENT */
L900:
    if (! pos) {
    goto L910;
    }

/*     SWAP ROWS MEP1 AND NIV, AND SCALE FACTORS FOR THESE ROWS. */
    sswap_(&np1, &w[mep1 + w_dim1], mdw, &w[niv + w_dim1], mdw);
    sswap_(&c__1, &scale[mep1], &c__1, &scale[niv], &c__1);
    itemp = itype[mep1];
    itype[mep1] = itype[niv];
    itype[niv] = itemp;
    me = mep1;
    mep1 = me + 1;
L910:
    goto L930;
L920:
    pos = FALSE_;
L930:
L940:
    if (pos) {
    goto L950;
    }
    nsp1 = nsoln;
    --nsoln;
    niv1 = niv;
    --niv;
L950:
L960:
    goto L740;
L970:
    goto L1250;
L980:

/*     TO TEST-PROPOSED-NEW-COMPONENT */
    z2 = w[isol + np1 * w_dim1] / w[isol + nsoln * w_dim1];
    z__[nsoln] = z2;
    pos = z2 > zero;
    if (! (z2 * eanorm >= bnorm && pos)) {
    goto L990;
    }
    pos = ! (blowup * z2 * eanorm >= bnorm);
L990:
    goto L1280;
L1000:
/*     TO COMPUTE-FINAL-SOLUTION */

/*     SOLVE SYSTEM, STORE RESULTS IN X(*). */

    igo958 = 1;
//    igo958_fmt = fmt_1010;
    goto L1110;
/*     SOLVE-SYSTEM */
L1010:
    scopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);

/*     APPLY HOUSEHOLDER TRANSFORMATIONS TO X(*) IF KRANK.LT.L */
    if (! (0 < krank && krank < *l)) {
    goto L1030;
    }
    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
    h12_(&c__2, &i__, &krp1, l, &w[i__ + w_dim1], mdw, &h__[i__], &x[1], &c__1, &c__1, &c__1);
/* L1020: */
    }

/*     FILL IN TRAILING ZEROES FOR CONSTRAINED VARIABLES NOT IN SOLN. */
L1030:
    if (! (nsoln < *n)) {
    goto L1040;
    }
    x[nsp1] = zero;
    i__1 = *n - nsoln;
    scopy_(&i__1, &x[nsp1], &c__0, &x[nsp1], &c__1);

/*     REPERMUTE SOLN VECTOR TO NATURAL ORDER. */
L1040:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    j = i__;
L1050:
    if (ipivot[j] == i__) {
        goto L1060;
    }
    ++j;
    goto L1050;
L1060:
    ipivot[j] = ipivot[i__];
    ipivot[i__] = j;
    sswap_(&c__1, &x[j], &c__1, &x[i__], &c__1);
/* L1070: */
    }

/*     RESCALE THE SOLN USING THE COL SCALING. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
    x[j] *= d__[j];
/* L1080: */
    }
/*     IF (.NOT.(NSOLN.LT.M)) GO TO 1100                           REMK */
/*     DO 1090 I=NSP1,M                                            REMK */
    if (! (niv < m)) {
    goto L1100;
    }
    i__1 = m;
    for (i__ = niv1; i__ <= i__1; ++i__) {
    t = w[i__ + np1 * w_dim1];
    if (i__ <= me) {
        t /= alamda;
    }
    t = scale[i__] * t * t;
    *rnorm += t;
/* L1090: */
    }
L1100:
    *rnorm = sqrt(*rnorm);
    goto L1210;

/*     TO SOLVE-SYSTEM */

L1110:
    if (! done) {
    goto L1120;
    }
    isol = 1;
    goto L1130;
L1120:
    isol = lp1;
L1130:
    if (! (nsoln >= isol)) {
    goto L1190;
    }

/*     COPY RT. HAND SIDE INTO TEMP VECTOR TO USE OVERWRITING METHOD. */
    scopy_(&niv, &w[np1 * w_dim1 + 1], &c__1, &temp[1], &c__1);
    i__1 = nsoln;
    for (jj = isol; jj <= i__1; ++jj) {
    j = nsoln - jj + isol;
    if (! (j > krank)) {
        goto L1140;
    }
    i__ = niv - jj + isol;
    goto L1150;
L1140:
    i__ = j;
L1150:
    if (! (j > krank && j <= *l)) {
        goto L1160;
    }
    z__[j] = zero;
    goto L1170;
L1160:
    z__[j] = temp[i__] / w[i__ + j * w_dim1];
    i__2 = i__ - 1;
    r__1 = -z__[j];
    saxpy_(&i__2, &r__1, &w[j * w_dim1 + 1], &c__1, &temp[1], &c__1);
L1170:
/* L1180: */
    ;
    }
L1190:
    goto L1270;
L1200:
    switch (igo986) {
    case 0: goto L40;
    }
L1210:
    switch (igo980) {
    case 0: goto L90;
    }
L1220:
    switch (igo991) {
    case 0: goto L30;
    }
L1230:
    switch (igo998) {
    case 0: goto L10;
    }
L1240:
    switch (igo995) {
    case 0: goto L20;
    }
L1250:
    switch (igo983) {
    case 0: goto L60;
    }
L1260:
    switch (igo977) {
    case 0: goto L220;
    }
L1270:
    switch (igo958) {
    case 0: goto L310;
    case 1: goto L1010;
    }
L1280:
    switch (igo897) {
    case 0: goto L880;
    case 1: goto L900;
    }
L1290:
    switch (igo938) {
    case 0: goto L460;
    case 1: goto L480;
    case 2: goto L500;
    case 3: goto L520;
    }
} /* wnlsm_ */

/* Subroutine */ 
int lsei::wnlit_(double *w, integer *mdw, integer *m, integer *n, 
    integer *l, integer *ipivot, integer *itype, double *h__, double *scale, 
    double *rnorm, integer *idope, double *dope, logical *done)
{
    /* Initialized data */

    double tenm3 = .001;
    double zero = 0.;
    double one = 1.;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    double r__1, r__2, r__3;

    /* Builtin functions */
 //   double sqrt(doublereal);

    /* Local variables */
    double hbar;
    integer lend, mend;
    double alsq;
    integer igo990, igo993, igo996, i__, j, k;
    double t;
    logical indep;
    integer krank, itemp, nsoln, i1, j1, l1;
    integer ic, lb, me, jj, kk;
    logical recalc;
    integer jp, ir;
    double rn, sn, factor, eanorm;    
    double sparam[5];
    integer jm1, ip1;    
    integer lp1, np1, max__;
    double tau;
    integer niv, mep1, irp1, krp1, niv1;

/*     double             ALSQ, AMAX, EANORM, FAC, FACTOR, HBAR, ONE, RN */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --ipivot;
    --itype;
    --h__;
    --scale;
    --idope;
    --dope;

    /* Function Body */

    me = idope[1];
    mep1 = idope[2];
    krank = idope[3];
    krp1 = idope[4];
    nsoln = idope[5];
    niv = idope[6];
    niv1 = idope[7];
    l1 = idope[8];

    alsq = dope[1];
    eanorm = dope[2];
/*     FAC = DOPE(3)                                               REMK */
    tau = dope[4];
    np1 = *n + 1;
/* Computing MIN */
    i__1 = *m - 1;
    lb = min(i__1,*l);
    recalc = TRUE_;
    *rnorm = zero;
    krank = 0;
/*     WE SET FACTOR=1.E0 SO THAT THE HEAVY WEIGHT ALAMDA WILL BE */
/*     INCLUDED IN THE TEST FOR COL INDEPENDENCE. */
    factor = 1.;
    i__ = 1;
    ip1 = 2;
    lend = *l;
L10:
    if (! (i__ <= lb)) {
    goto L150;
    }
    if (! (i__ <= me)) {
    goto L130;
    }

/*     SET IR TO POINT TO THE I-TH ROW. */
    ir = i__;
    mend = *m;
    igo996 = 0;
//    igo996_fmt = fmt_20;
    goto L460;

/*     UPDATE-COL-SS-AND-FIND-PIVOT-COL */
L20:
    igo993 = 0;
//    igo993_fmt = fmt_30;
    goto L560;

/*     PERFORM-COL-INTERCHANGE */

/*     SET IC TO POINT TO I-TH COL. */
L30:
    ic = i__;
    igo990 = 0;
//    igo990_fmt = fmt_40;
    goto L520;

/*     TEST-INDEP-OF-INCOMING-COL */
L40:
    if (! indep) {
    goto L110;
    }

/*     ELIMINATE I-TH COL BELOW DIAG. USING MOD. GIVENS TRANSFORMATIONS */
/*     APPLIED TO (A B). */
    j = *m;
    i__1 = *m;
    for (jj = ip1; jj <= i__1; ++jj) {
    jm1 = j - 1;
    jp = jm1;
/*     WHEN OPERATING NEAR THE ME LINE, USE THE LARGEST ELT.        REMK */
/*     ABOVE IT AS THE PIVOT.                                       REMK */
/*       IF (.NOT.(J.EQ.MEP1)) GO TO 80                             REMK */
/*       IMAX = ME                                                  REMK */
/*       AMAX = SCALE(ME)*W(ME,I)**2                                REMK */
/*  50   IF (.NOT.(JP.GE.I)) GO TO 70                               REMK */
/*       T = SCALE(JP)*W(JP,I)**2                                   REMK */
/*       IF (.NOT.(T.GT.AMAX)) GO TO 60                             REMK */
/*       IMAX = JP                                                  REMK */
/*       AMAX = T                                                   REMK */
/*  60   JP = JP - 1                                                REMK */
/*       GO TO 50                                                   REMK */
/*  70   JP = IMAX                                                  REMK */
    if (! (jj == *m)) {
        goto L70;
    }
    if (! (i__ < mep1)) {
        goto L80;
    }
    j = mep1;
    jp = i__;
/* Computing 2nd power */
    r__1 = w[jp + i__ * w_dim1];
/* Computing 2nd power */
    r__2 = tau;
    t = scale[jp] * (r__1 * r__1) * (r__2 * r__2);
/* Computing 2nd power */
    r__1 = w[j + i__ * w_dim1];
    if (! (t > scale[j] * (r__1 * r__1))) {
        goto L130;
    }
    goto L80;
L70:
    if (! (j == mep1)) {
        goto L80;
    }
    j = jm1;
    jm1 = j - 1;
    jp = jm1;
L80:
    if (! (w[j + i__ * w_dim1] != zero)) {
        goto L90;
    }
    srotmg_(&scale[jp], &scale[j], &w[jp + i__ * w_dim1], &w[j + i__ *  w_dim1], sparam);
    w[j + i__ * w_dim1] = zero;
    i__2 = np1 - i__;
    srotm_(&i__2, &w[jp + ip1 * w_dim1], mdw, &w[j + ip1 * w_dim1], mdw,    sparam);
L90:
    j = jm1;
/* L100: */
    }
    goto L140;
L110:
    if (! (lend > i__)) {
    goto L130;
    }

/*     COL I IS DEPENDENT. SWAP WITH COL LEND. */
    max__ = lend;

/*     PERFORM-COL-INTERCHANGE */
    igo993 = 1;
//    igo993_fmt = fmt_120;
    goto L560;
L120:
    --lend;

/*     FIND COL IN REMAINING SET WITH LARGEST SS. */
    i__1 = lend - i__ + 1;
    max__ = isamax_(&i__1, &h__[i__], &c__1) + i__ - 1;
    hbar = h__[max__];
    goto L30;
L130:
    krank = i__ - 1;
    goto L160;
L140:
    i__ = ip1;
    ++ip1;
    goto L10;
L150:
    krank = l1;
L160:
    krp1 = krank + 1;
    if (! (krank < me)) {
    goto L290;
    }
    factor = alsq;
    i__1 = me;
    for (i__ = krp1; i__ <= i__1; ++i__) {
    if (*l > 0) {
        w[i__ + w_dim1] = zero;
    }
    scopy_(l, &w[i__ + w_dim1], &c__0, &w[i__ + w_dim1], mdw);
/* L170: */
    }

/*     DETERMINE THE RANK OF THE REMAINING EQUALITY CONSTRAINT */
/*     EQUATIONS BY ELIMINATING WITHIN THE BLOCK OF CONSTRAINED */
/*     VARIABLES.  REMOVE ANY REDUNDANT CONSTRAINTS. */
    ir = krp1;
    if (! (*l < *n)) {
    goto L245;
    }
    lp1 = *l + 1;
    recalc = TRUE_;
/* Computing MIN */
    i__1 = *l + me - krank;
    lb = min(i__1,*n);
    i__ = lp1;
    ip1 = i__ + 1;
L180:
    if (! (i__ <= lb)) {
    goto L280;
    }
    ir = krank + i__ - *l;
    lend = *n;
    mend = me;
    igo996 = 1;
//    igo996_fmt = fmt_190;
    goto L460;

/*     UPDATE-COL-SS-AND-FIND-PIVOT-COL */
L190:
    igo993 = 2;
//    igo993_fmt = fmt_200;
    goto L560;

/*     PERFORM-COL-INTERCHANGE */

/*     ELIMINATE ELEMENTS IN THE I-TH COL. */
L200:
    j = me;
L210:
    if (! (j > ir)) {
    goto L230;
    }
    jm1 = j - 1;
    if (! (w[j + i__ * w_dim1] != zero)) {
    goto L220;
    }
    srotmg_(&scale[jm1], &scale[j], &w[jm1 + i__ * w_dim1], &w[j + i__ *  w_dim1], sparam);
    w[j + i__ * w_dim1] = zero;
    i__1 = np1 - i__;
    srotm_(&i__1, &w[jm1 + ip1 * w_dim1], mdw, &w[j + ip1 * w_dim1], mdw,   sparam);
L220:
    j = jm1;
    goto L210;

/*     SET IC=I=COL BEING ELIMINATED */
L230:
    ic = i__;
    igo990 = 1;
//    igo990_fmt = fmt_240;
    goto L520;

/*     TEST-INDEP-OF-INCOMING-COL */
L240:
    if (indep) {
    goto L270;
    }

/*     REMOVE ANY REDUNDANT OR DEPENDENT EQUALITY CONSTRAINTS. */
L245:
    jj = ir;
L250:
    if (! (ir <= me)) {
    goto L260;
    }
    w[ir + w_dim1] = zero;
    scopy_(n, &w[ir + w_dim1], &c__0, &w[ir + w_dim1], mdw);
    *rnorm += scale[ir] * w[ir + np1 * w_dim1] / alsq * w[ir + np1 * w_dim1];
    w[ir + np1 * w_dim1] = zero;
    scale[ir] = one;
/*     RECLASSIFY THE ZEROED ROW AS A LEAST SQUARES EQUATION. */
    itype[ir] = 1;
    ++ir;
    goto L250;

/*     REDUCE ME TO REFLECT ANY DISCOVERED DEPENDENT EQUALITY */
/*     CONSTRAINTS. */
L260:
    me = jj - 1;
    mep1 = me + 1;
    goto L300;
L270:
    i__ = ip1;
    ++ip1;
    goto L180;
L280:
L290:
L300:
    if (! (krank < l1)) {
    goto L420;
    }

/*     TRY TO DETERMINE THE VARIABLES KRANK+1 THROUGH L1 FROM THE */
/*     LEAST SQUARES EQUATIONS.  CONTINUE THE TRIANGULARIZATION WITH */
/*     PIVOT ELEMENT W(MEP1,I). */

    recalc = TRUE_;

/*     SET FACTOR=ALSQ TO REMOVE EFFECT OF HEAVY WEIGHT FROM */
/*     TEST FOR COL INDEPENDENCE. */
    factor = alsq;
    kk = krp1;
    i__ = kk;
    ip1 = i__ + 1;
L310:
    if (! (i__ <= l1)) {
    goto L410;
    }

/*     SET IR TO POINT TO THE MEP1-ST ROW. */
    ir = mep1;
    lend = *l;
    mend = *m;
    igo996 = 2;
//    igo996_fmt = fmt_320;
    goto L460;

/*     UPDATE-COL-SS-AND-FIND-PIVOT-COL */
L320:
    igo993 = 3;
//    igo993_fmt = fmt_330;
    goto L560;

/*     PERFORM-COL-INTERCHANGE */

/*     ELIMINATE I-TH COL BELOW THE IR-TH ELEMENT. */
L330:
    irp1 = ir + 1;
    if (! (irp1 <= *m)) {
    goto L355;
    }
    j = *m;
    i__1 = *m;
    for (jj = irp1; jj <= i__1; ++jj) {
    jm1 = j - 1;
    if (! (w[j + i__ * w_dim1] != zero)) {
        goto L340;
    }
    srotmg_(&scale[jm1], &scale[j], &w[jm1 + i__ * w_dim1], &w[j + i__ *    w_dim1], sparam);
    w[j + i__ * w_dim1] = zero;
    i__2 = np1 - i__;
    srotm_(&i__2, &w[jm1 + ip1 * w_dim1], mdw, &w[j + ip1 * w_dim1], mdw, sparam);
L340:
    j = jm1;
/* L350: */
    }
L355:

/*     TEST IF NEW PIVOT ELEMENT IS NEAR ZERO. IF SO, THE COL IS */
/*     DEPENDENT. */
/* Computing 2nd power */
    r__1 = w[ir + i__ * w_dim1];
    t = scale[ir] * (r__1 * r__1);
/* Computing 2nd power */
    r__1 = tau;
/* Computing 2nd power */
    r__2 = eanorm;
    indep = t > r__1 * r__1 * (r__2 * r__2);
    if (! indep) {
    goto L380;
    }

/*     COL TEST PASSED. NOW MUST PASS ROW NORM TEST TO BE CLASSIFIED */
/*     AS INDEPENDENT. */
    rn = zero;
    i__1 = *m;
    for (i1 = ir; i1 <= i__1; ++i1) {
    i__2 = *n;
    for (j1 = ip1; j1 <= i__2; ++j1) {
/* Computing MAX */
/* Computing 2nd power */
        r__3 = w[i1 + j1 * w_dim1];
        r__1 = rn, r__2 = scale[i1] * (r__3 * r__3);
        rn = dmax(r__1,r__2);
/* L360: */
    }
/* L370: */
    }
/* Computing 2nd power */
    r__1 = tau;
    indep = t > r__1 * r__1 * rn;

/*     IF INDEPENDENT, SWAP THE IR-TH AND KRP1-ST ROWS TO MAINTAIN THE */
/*     TRIANGULAR FORM.  UPDATE THE RANK INDICATOR KRANK AND THE */
/*     EQUALITY CONSTRAINT POINTER ME. */
L380:
    if (! indep) {
    goto L390;
    }
    sswap_(&np1, &w[krp1 + w_dim1], mdw, &w[ir + w_dim1], mdw);
    sswap_(&c__1, &scale[krp1], &c__1, &scale[ir], &c__1);
/*     RECLASSIFY THE LEAST SQ. EQUATION AS AN EQUALITY CONSTRAINT AND */
/*     RESCALE IT. */
    itype[ir] = 0;
    t = sqrt(scale[krp1]);
    sscal_(&np1, &t, &w[krp1 + w_dim1], mdw);
    scale[krp1] = alsq;
    me = mep1;
    mep1 = me + 1;
    krank = krp1;
    krp1 = krank + 1;
    goto L400;
L390:
    goto L430;
L400:
    i__ = ip1;
    ++ip1;
    goto L310;
L410:
L420:
L430:

/*     IF PSEUDORANK IS LESS THAN L, APPLY HOUSEHOLDER TRANS. */
/*     FROM RIGHT. */
    if (! (krank < *l)) {
    goto L450;
    }
    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
    j = krp1 - i__;
    i__2 = j - 1;
    h12_(&c__1, &j, &krp1, l, &w[j + w_dim1], mdw, &h__[j], &w[w_offset], 
        mdw, &c__1, &i__2);
/* L440: */
    }
L450:
    niv = krank + nsoln - *l;
    niv1 = niv + 1;
    if (*l == *n) {
    *done = TRUE_;
    }

/*  END OF INITIAL TRIANGULARIZATION. */
    idope[1] = me;
    idope[2] = mep1;
    idope[3] = krank;
    idope[4] = krp1;
    idope[5] = nsoln;
    idope[6] = niv;
    idope[7] = niv1;
    idope[8] = l1;
    return 0;
L460:

/*     TO UPDATE-COL-SS-AND-FIND-PIVOT-COL */

/*     THE COL SS VECTOR WILL BE UPDATED AT EACH STEP. WHEN */
/*     NUMERICALLY NECESSARY, THESE VALUES WILL BE RECOMPUTED. */

    if (! (ir != 1 && ! recalc)) {
    goto L480;
    }
/*     UPDATE COL SS =SUM OF SQUARES. */
    i__1 = lend;
    for (j = i__; j <= i__1; ++j) {
/* Computing 2nd power */
    r__1 = w[ir - 1 + j * w_dim1];
    h__[j] -= scale[ir - 1] * (r__1 * r__1);
/* L470: */
    }

/*     TEST FOR NUMERICAL ACCURACY. */
    i__1 = lend - i__ + 1;
    max__ = isamax_(&i__1, &h__[i__], &c__1) + i__ - 1;
    recalc = hbar + tenm3 * h__[max__] == hbar;

/*     IF REQUIRED, RECALCULATE COL SS, USING ROWS IR THROUGH MEND. */
L480:
    if (! recalc) {
    goto L510;
    }
    i__1 = lend;
    for (j = i__; j <= i__1; ++j) {
    h__[j] = zero;
    i__2 = mend;
    for (k = ir; k <= i__2; ++k) {
/* Computing 2nd power */
        r__1 = w[k + j * w_dim1];
        h__[j] += scale[k] * (r__1 * r__1);
/* L490: */
    }
/* L500: */
    }

/*     FIND COL WITH LARGEST SS. */
    i__1 = lend - i__ + 1;
    max__ = isamax_(&i__1, &h__[i__], &c__1) + i__ - 1;
    hbar = h__[max__];
L510:
    goto L600;
L520:

/*     TO TEST-INDEP-OF-INCOMING-COL */

/*     TEST THE COL IC TO DETERMINE IF IT IS LINEARLY INDEPENDENT */
/*     OF THE COLS ALREADY IN THE BASIS.  IN THE INIT TRI */
/*     STEP, WE USUALLY WANT THE HEAVY WEIGHT ALAMDA TO */
/*     BE INCLUDED IN THE TEST FOR INDEPENDENCE.  IN THIS CASE THE */
/*     VALUE OF FACTOR WILL HAVE BEEN SET TO 1.E0 BEFORE THIS */
/*     PROCEDURE IS INVOKED.  IN THE POTENTIALLY RANK DEFICIENT */
/*     PROBLEM, THE VALUE OF FACTOR WILL HAVE BEEN */
/*     SET TO ALSQ=ALAMDA**2 TO REMOVE THE EFFECT OF THE HEAVY WEIGHT */
/*     FROM THE TEST FOR INDEPENDENCE. */

/*     WRITE NEW COL AS PARTITIONED VECTOR */
/*             (A1)  NUMBER OF COMPONENTS IN SOLN SO FAR = NIV */
/*             (A2)  M-NIV COMPONENTS */
/*     AND COMPUTE  SN = INVERSE WEIGHTED LENGTH OF A1 */
/*                  RN = INVERSE WEIGHTED LENGTH OF A2 */
/*     CALL THE COL INDEPENDENT WHEN RN .GT. TAU*SN */
    sn = zero;
    rn = zero;
    i__1 = mend;
    for (j = 1; j <= i__1; ++j) {
    t = scale[j];
    if (j <= me) {
        t /= factor;
    }
/* Computing 2nd power */
    r__1 = w[j + ic * w_dim1];
    t *= r__1 * r__1;
    if (! (j < ir)) {
        goto L530;
    }
    sn += t;
    goto L540;
L530:
    rn += t;
L540:
/* L550: */
    ;
    }
/* Computing 2nd power */
    r__1 = tau;
    indep = rn > r__1 * r__1 * sn;
    goto L590;
L560:

/*     TO PERFORM-COL-INTERCHANGE */

    if (! (max__ != i__)) {
    goto L570;
    }
/*     EXCHANGE ELEMENTS OF PERMUTED INDEX VECTOR AND PERFORM COL */
/*     INTERCHANGES. */
    itemp = ipivot[i__];
    ipivot[i__] = ipivot[max__];
    ipivot[max__] = itemp;
    sswap_(m, &w[max__ * w_dim1 + 1], &c__1, &w[i__ * w_dim1 + 1], &c__1);
    t = h__[max__];
    h__[max__] = h__[i__];
    h__[i__] = t;
L570:
    goto L580;
L580:
    switch (igo993) {
    case 0: goto L30;
    case 1: goto L120;
    case 2: goto L200;
    case 3: goto L330;
    }
L590:
    switch (igo990) {
    case 0: goto L40;
    case 1: goto L240;
    }
L600:
    switch (igo996) {
    case 0: goto L20;
    case 1: goto L190;
    case 2: goto L320;
    }
} /* wnlit_ */

/* Subroutine */ 
int lsei::hfti_(double *a, integer *mda, integer *m, integer *n, 
    double *b, integer *mdb, integer *nb, double *tau, integer *krank, double *
    rnorm, double *h__, double *g, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    double r__1;

    /* Builtin functions */
//    double sqrt(doublereal);

    /* Local variables */
    
    integer jcol;
    double hmax;
    integer lmax, nerr, iopt;
    double zero;
    integer i__, j, k, l, ldiag;
    integer jb, ii, jj;
    double sm, factor;
    integer ip1, kp1;
    double sm1;
    double tmp;

    /* Parameter adjustments */
    --ip;
    --g;
    --h__;
    a_dim1 = *mda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --rnorm;

    /* Function Body */
    zero = 0.;
    factor = .001;

    k = 0;
    ldiag = min(*m,*n);
    if (ldiag <= 0) {
    goto L310;
    }
    if (! (*mda < *m)) {
    goto L10;
    }
    nerr = 2;
    iopt = 2;
    xerror_("HFTI MDA.LT.M.. PROBABLE ERROR.", &c__31, &nerr, &iopt/*, (ftnlen)    31*/);
    return 0;
L10:

    if (! (*nb > 1 && max(*m,*n) > *mdb)) {
    goto L20;
    }
    nerr = 2;
    iopt = 2;
    xerror_("HFTI MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.", &c__49, &   nerr, &iopt/*, (ftnlen)49*/);
    return 0;
L20:

    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
    if (j == 1) {
        goto L40;
    }

/*     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*    .. */
    lmax = j;
    i__2 = *n;
    for (l = j; l <= i__2; ++l) {
/* Computing 2nd power */
        r__1 = a[j - 1 + l * a_dim1];
        h__[l] -= r__1 * r__1;
        if (h__[l] > h__[lmax]) {
        lmax = l;
        }
/* L30: */
    }
    r__1 = hmax + factor * h__[lmax];
    if (diff_(&r__1, &hmax) <= 0.) {
        goto L40;
    } else {
        goto L70;
    }

/*     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*    .. */
L40:
    lmax = j;
    i__2 = *n;
    for (l = j; l <= i__2; ++l) {
        h__[l] = zero;
        i__3 = *m;
        for (i__ = j; i__ <= i__3; ++i__) {
/* Computing 2nd power */
        r__1 = a[i__ + l * a_dim1];
        h__[l] += r__1 * r__1;
/* L50: */
        }
        if (h__[l] > h__[lmax]) {
        lmax = l;
        }
/* L60: */
    }
    hmax = h__[lmax];
/*    .. */
/*     LMAX HAS BEEN DETERMINED */

/*     DO COLUMN INTERCHANGES IF NEEDED. */
/*    .. */
L70:
    ip[j] = lmax;
    if (ip[j] == j) {
        goto L90;
    }
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
        tmp = a[i__ + j * a_dim1];
        a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
        a[i__ + lmax * a_dim1] = tmp;
/* L80: */
    }
    h__[lmax] = h__[j];
L90:
/* Computing MIN */
    i__2 = j + 1;
    jcol = min(i__2,*n);

/*     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B. */
/*    .. */
    i__2 = j + 1;
    i__3 = *n - j;
    h12_(&c__1, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &a[jcol * a_dim1 + 1], &c__1, mda, &i__3);
    i__2 = j + 1;
    h12_(&c__2, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &b[b_offset], &c__1, mdb, nb);
/* L100: */
    }

/*     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU. */
/*    .. */
    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
    if ((r__1 = a[j + j * a_dim1], dabs(r__1)) <= *tau) {
        goto L120;
    }
/* L110: */
    }
    k = ldiag;
    goto L130;
L120:
    k = j - 1;
L130:
    kp1 = k + 1;

/*     COMPUTE THE NORMS OF THE RESIDUAL VECTORS. */

    if (*nb <= 0) {
    goto L170;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
    tmp = zero;
    if (kp1 > *m) {
        goto L150;
    }
    i__2 = *m;
    for (i__ = kp1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
        r__1 = b[i__ + jb * b_dim1];
        tmp += r__1 * r__1;
/* L140: */
    }
L150:
    rnorm[jb] = sqrt(tmp);
/* L160: */
    }
L170:
/*                                           SPECIAL FOR PSEUDORANK = 0 */
    if (k > 0) {
    goto L200;
    }
    if (*nb <= 0) {
    goto L310;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + jb * b_dim1] = zero;
/* L180: */
    }
/* L190: */
    }
    goto L310;

/*     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER */
/*     DECOMPOSITION OF FIRST K ROWS. */
/*    .. */
L200:
    if (k == *n) {
    goto L220;
    }
    i__1 = k;
    for (ii = 1; ii <= i__1; ++ii) {
    i__ = kp1 - ii;
    i__2 = i__ - 1;
    h12_(&c__1, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &a[
        a_offset], mda, &c__1, &i__2);
/* L210: */
    }
L220:


    if (*nb <= 0) {
    goto L310;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {

/*     SOLVE THE K BY K TRIANGULAR SYSTEM. */
/*    .. */
    i__2 = k;
    for (l = 1; l <= i__2; ++l) {
        sm = zero;
        i__ = kp1 - l;
        if (i__ == k) {
        goto L240;
        }
        ip1 = i__ + 1;
        i__3 = k;
        for (j = ip1; j <= i__3; ++j) {
        sm += a[i__ + j * a_dim1] * b[j + jb * b_dim1];
/* L230: */
        }
L240:
        sm1 = sm;
        b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1] - sm1) / a[i__ + i__ * a_dim1];
/* L250: */
    }

/*     COMPLETE COMPUTATION OF SOLUTION VECTOR. */
/*    .. */
    if (k == *n) {
        goto L280;
    }
    i__2 = *n;
    for (j = kp1; j <= i__2; ++j) {
        b[j + jb * b_dim1] = zero;
/* L260: */
    }
    i__2 = k;
    for (i__ = 1; i__ <= i__2; ++i__) {
        h12_(&c__2, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &b[jb * b_dim1 + 1], &c__1, mdb, &c__1);
/* L270: */
    }

/*      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE */
/*      COLUMN INTERCHANGES. */
/*    .. */
L280:
    i__2 = ldiag;
    for (jj = 1; jj <= i__2; ++jj) {
        j = ldiag + 1 - jj;
        if (ip[j] == j) {
        goto L290;
        }
        l = ip[j];
        tmp = b[l + jb * b_dim1];
        b[l + jb * b_dim1] = b[j + jb * b_dim1];
        b[j + jb * b_dim1] = tmp;
L290:
        ;
    }
/* L300: */
    }
/*    .. */
/*     THE SOLUTION VECTORS, X, ARE NOW */
/*     IN THE FIRST  N  ROWS OF THE ARRAY B(,). */

L310:
    *krank = k;
    return 0;
} /* hfti_ */

/* Subroutine */ int lsei::h12_(integer *mode, integer *lpivot, integer *l1, 
    integer *m, double *u, integer *iue, double *up, double *c__, integer *ice, 
    integer *icv, integer *ncv)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;
    double r__1, r__2;

    /* Builtin functions */
 //   double sqrt(doublereal);

    /* Local variables */
    integer incr;
    double ul1m1;
    double b;
    integer i__, j, mml1p2;
    double clinv;
    integer i2, i3, i4;
    double cl, sm;
    integer kl1, kl2, l1m1;
    double one;
    integer klp;

    /* Parameter adjustments */
    u_dim1 = *iue;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --c__;

    /* Function Body */
    one = 1.;

    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
    return 0;
    }
    cl = (r__1 = u[*lpivot * u_dim1 + 1], dabs(r__1));
    if (*mode == 2) {
    goto L60;
    }
/*                            ****** CONSTRUCT THE TRANSFORMATION. ****** */
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L10: */
/* Computing MAX */
    r__2 = (r__1 = u[j * u_dim1 + 1], dabs(r__1));
    cl = dmax(r__2,cl);
    }
    if (cl <= 0.) {
    goto L130;
    } else {
    goto L20;
    }
L20:
    clinv = one / cl;
/* Computing 2nd power */
    r__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = r__1 * r__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L30: */
/* Computing 2nd power */
    r__1 = u[j * u_dim1 + 1] * clinv;
    sm += r__1 * r__1;
    }
    cl *= sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] <= 0.) {
    goto L50;
    } else {
    goto L40;
    }
L40:
    cl = -cl;
L50:
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L70;
/*            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** */

L60:
    if (cl <= 0.) {
    goto L130;
    } else {
    goto L70;
    }
L70:
    if (*ncv <= 0) {
    return 0;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
/*                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. */

    if (b >= 0.) {
    goto L130;
    } else {
    goto L80;
    }
L80:
    b = one / b;
    mml1p2 = *m - *l1 + 2;
    if (mml1p2 > 20) {
    goto L140;
    }
    i2 = 1 - *icv + *ice * (*lpivot - 1);
    incr = *ice * (*l1 - *lpivot);
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
    i2 += *icv;
    i3 = i2 + incr;
    i4 = i3;
    sm = c__[i2] * *up;
    i__2 = *m;
    for (i__ = *l1; i__ <= i__2; ++i__) {
        sm += c__[i3] * u[i__ * u_dim1 + 1];
/* L90: */
        i3 += *ice;
    }
    if (sm != 0.) {
        goto L100;
    } else {
        goto L120;
    }
L100:
    sm *= b;
    c__[i2] += sm * *up;
    i__2 = *m;
    for (i__ = *l1; i__ <= i__2; ++i__) {
        c__[i4] += sm * u[i__ * u_dim1 + 1];
/* L110: */
        i4 += *ice;
    }
L120:
    ;
    }
L130:
    return 0;
L140:
    l1m1 = *l1 - 1;
    kl1 = (l1m1 - 1) * *ice + 1;
    kl2 = kl1;
    klp = (*lpivot - 1) * *ice + 1;
    ul1m1 = u[l1m1 * u_dim1 + 1];
    u[l1m1 * u_dim1 + 1] = *up;
    if (*lpivot == l1m1) {
    goto L150;
    }
    sswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
L150:
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
    sm = sdot_(&mml1p2, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
    sm *= b;
    saxpy_(&mml1p2, &sm, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
    kl1 += *icv;
/* L160: */
    }
    u[l1m1 * u_dim1 + 1] = ul1m1;
    if (*lpivot == l1m1) {
    return 0;
    }
    kl1 = kl2;
    sswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
    return 0;
} /* h12_ */




/* Subroutine */ 
int lsei::srotmg_(double *sd1, double *sd2, double *sx1, double *sy1, double *sparam)
{
    /* Initialized data */

    double zero = 0.;
    double one = 1.;
    double two = 2.;
    integer iflag = 1;
    double gam = 4096.;
    double gamsq = 1.678e7;
    double rgam = 2.441e-4;
    double rgamsq = 5.96e-8;

    /* System generated locals */
    double r__1;

    /* Local variables */
    double sflag, stemp, su, sp1, sp2, sq2, sq1, sh11, sh21, sh12, sh22;
    integer igo;

    /* Parameter adjustments */
    --sparam;

    /* Function Body */

    if (! (*sd1 < zero)) {
    goto L10;
    }
/*       GO ZERO-H-D-AND-SX1.. */
    goto L60;
L10:
/*     CASE-SD1-NONNEGATIVE */
    sp2 = *sd2 * *sy1;
    if (! (sp2 == zero)) {
    goto L20;
    }
    sflag = -two;
    goto L260;
/*     REGULAR-CASE.. */
L20:
    sp1 = *sd1 * *sx1;
    sq2 = sp2 * *sy1;
    sq1 = sp1 * *sx1;

    if (! (dabs(sq1) > dabs(sq2))) {
    goto L40;
    }
    sh21 = -(*sy1) / *sx1;
    sh12 = sp2 / sp1;

    su = one - sh12 * sh21;

    if (! (su <= zero)) {
    goto L30;
    }
/*         GO ZERO-H-D-AND-SX1.. */
    goto L60;
L30:
    sflag = zero;
    *sd1 /= su;
    *sd2 /= su;
    *sx1 *= su;
/*         GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (sq2 < zero)) {
    goto L50;
    }
/*         GO ZERO-H-D-AND-SX1.. */
    goto L60;
L50:
    sflag = one;
    sh11 = sp1 / sp2;
    sh22 = *sx1 / *sy1;
    su = one + sh11 * sh22;
    stemp = *sd2 / su;
    *sd2 = *sd1 / su;
    *sd1 = stemp;
    *sx1 = *sy1 * su;
/*         GO SCALE-CHECK */
    goto L100;
/*     PROCEDURE..ZERO-H-D-AND-SX1.. */
L60:
    sflag = -one;
    sh11 = zero;
    sh12 = zero;
    sh21 = zero;
    sh22 = zero;

    *sd1 = zero;
    *sd2 = zero;
    *sx1 = zero;
/*         RETURN.. */
    goto L220;
/*     PROCEDURE..FIX-H.. */
L70:
    if (! (sflag >= zero)) {
    goto L90;
    }

    if (! (sflag == zero)) {
    goto L80;
    }
    sh11 = one;
    sh22 = one;
    sflag = -one;
    goto L90;
L80:
    sh21 = -one;
    sh12 = one;
    sflag = -one;
L90:
    switch (igo) {
    case 0: goto L120;
    case 1: goto L150;
    case 2: goto L180;
    case 3: goto L210;
    }
/*     PROCEDURE..SCALE-CHECK */
L100:
    if (! (iflag == 1)) {
    goto L105;
    }

/*                   RECOMPUTE RESCALING PARAMETERS */
/*                   MORE ACCURATELY.. */

    rgam = one / gam;
/* Computing 2nd power */
    r__1 = gam;
    gamsq = r__1 * r__1;
/* Computing 2nd power */
    r__1 = rgam;
    rgamsq = r__1 * r__1;
    iflag = 2;
L105:
L110:
    if (! (*sd1 <= rgamsq)) {
    goto L130;
    }
    if (*sd1 == zero) {
    goto L160;
    }
    igo = 0;
//    igo_fmt = fmt_120;
/*              FIX-H.. */
    goto L70;
L120:
    *sd1 *= gamsq;
    *sx1 *= rgam;
    sh11 *= rgam;
    sh12 *= rgam;
    goto L110;
L130:
L140:
    if (! (*sd1 >= gamsq)) {
    goto L160;
    }
    igo = 1;
//    igo_fmt = fmt_150;
/*              FIX-H.. */
    goto L70;
L150:
    *sd1 *= rgamsq;
    *sx1 *= gam;
    sh11 *= gam;
    sh12 *= gam;
    goto L140;
L160:
L170:
    if (! (dabs(*sd2) <= rgamsq)) {
    goto L190;
    }
    if (*sd2 == zero) {
    goto L220;
    }
    igo = 2;
//    igo_fmt = fmt_180;
/*              FIX-H.. */
    goto L70;
L180:
    *sd2 *= gamsq;
    sh21 *= rgam;
    sh22 *= rgam;
    goto L170;
L190:
L200:
    if (! (dabs(*sd2) >= gamsq)) {
    goto L220;
    }
    igo = 3;
//    igo_fmt = fmt_210;
/*              FIX-H.. */
    goto L70;
L210:
    *sd2 *= rgamsq;
    sh21 *= gam;
    sh22 *= gam;
    goto L200;
L220:
    if (sflag < 0.) {
    goto L250;
    } else if (sflag == 0) {
    goto L230;
    } else {
    goto L240;
    }
L230:
    sparam[3] = sh21;
    sparam[4] = sh12;
    goto L260;
L240:
    sparam[2] = sh11;
    sparam[5] = sh22;
    goto L260;
L250:
    sparam[2] = sh11;
    sparam[3] = sh21;
    sparam[4] = sh12;
    sparam[5] = sh22;
L260:
    sparam[1] = sflag;
    return 0;
} /* srotmg_ */

/* Subroutine */ int lsei::scopy_(integer *n, double *sx, integer *incx, double *sy, 
    integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    sy[iy] = sx[ix];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7. */

L20:
    m = *n - *n / 7 * 7;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    sy[i__] = sx[i__];
/* L30: */
    }
    if (*n < 7) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
    sy[i__] = sx[i__];
    sy[i__ + 1] = sx[i__ + 1];
    sy[i__ + 2] = sx[i__ + 2];
    sy[i__ + 3] = sx[i__ + 3];
    sy[i__ + 4] = sx[i__ + 4];
    sy[i__ + 5] = sx[i__ + 5];
    sy[i__ + 6] = sx[i__ + 6];
/* L50: */
    }
    return 0;

/*        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    sy[i__] = sx[i__];
/* L70: */
    }
    return 0;
} /* scopy_ */

/* Subroutine */ int lsei::sswap_(integer *n, double *sx, integer *incx, double *sy, 
    integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m;
    double stemp1, stemp2, stemp3;
    integer ix, iy, ns, mp1;

    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    stemp1 = sx[ix];
    sx[ix] = sy[iy];
    sy[iy] = stemp1;
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3. */

L20:
    m = *n - *n / 3 * 3;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    stemp1 = sx[i__];
    sx[i__] = sy[i__];
    sy[i__] = stemp1;
/* L30: */
    }
    if (*n < 3) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
    stemp1 = sx[i__];
    stemp2 = sx[i__ + 1];
    stemp3 = sx[i__ + 2];
    sx[i__] = sy[i__];
    sx[i__ + 1] = sy[i__ + 1];
    sx[i__ + 2] = sy[i__ + 2];
    sy[i__] = stemp1;
    sy[i__ + 1] = stemp2;
    sy[i__ + 2] = stemp3;
/* L50: */
    }
    return 0;
L60:

/*     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    stemp1 = sx[i__];
    sx[i__] = sy[i__];
    sy[i__] = stemp1;
/* L70: */
    }
    return 0;
} /* sswap_ */

doublereal lsei::snrm2_(integer *n, double *sx, integer *incx)
{
    /* Initialized data */

    double zero = 0.;
    double one = 1.;
    double cutlo = 4.441e-16;
    double cuthi = 1.304e19;

    /* System generated locals */
    integer i__1, i__2;
    double ret_val, r__1;

    /* Builtin functions */
//    double sqrt(doublereal);

    /* Local variables */
    double xmax;
    integer next, i__, j, nn;
    double hitest, sum;

    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n > 0) {
    goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
//    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 BEGIN MAIN LOOP */
    i__ = 1;
L20:
    switch (next) {
    case 0: goto L30;
    case 1: goto L50;
    case 2: goto L70;
    case 3: goto L110;
    }
L30:
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
    goto L85;
    }
    next = 1;
//    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (sx[i__] == zero) {
    goto L200;
    }
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
    goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
    next = 2;
//    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i__ = j;
    next = 3;
//    next_fmt = fmt_110;
    sum = sum / sx[i__] / sx[i__];
L105:
    xmax = (r__1 = sx[i__], dabs(r__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
    goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
    if ((r__1 = sx[i__], dabs(r__1)) <= xmax) {
    goto L115;
    }
/* Computing 2nd power */
    r__1 = xmax / sx[i__];
    sum = one + sum * (r__1 * r__1);
    xmax = (r__1 = sx[i__], dabs(r__1));
    goto L200;

L115:
/* Computing 2nd power */
    r__1 = sx[i__] / xmax;
    sum += r__1 * r__1;
    goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;


/*     FOR double OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / (double) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
    if ((r__1 = sx[j], dabs(r__1)) >= hitest) {
        goto L100;
    }
/* L95: */
/* Computing 2nd power */
    r__1 = sx[j];
    sum += r__1 * r__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
    goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* snrm2_ */

doublereal lsei::sasum_(integer *n, double *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    double ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    integer i__, m, ns, mp1;

    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
    return ret_val;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    ret_val += (r__1 = sx[i__], dabs(r__1));
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6. */

L20:
    m = *n - *n / 6 * 6;
    if (m == 0) {
    goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
    ret_val += (r__1 = sx[i__], dabs(r__1));
/* L30: */
    }
    if (*n < 6) {
    return ret_val;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
    ret_val = ret_val + (r__1 = sx[i__], dabs(r__1)) + (r__2 = sx[i__ + 1]
        , dabs(r__2)) + (r__3 = sx[i__ + 2], dabs(r__3)) + (r__4 = sx[
        i__ + 3], dabs(r__4)) + (r__5 = sx[i__ + 4], dabs(r__5)) + (
        r__6 = sx[i__ + 5], dabs(r__6));
/* L50: */
    }
    return ret_val;
} /* sasum_ */

/* Subroutine */ int lsei::sscal_(integer *n, double *sa, double *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ns, mp1;

    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    sx[i__] = *sa * sx[i__];
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. */

L20:
    m = *n - *n / 5 * 5;
    if (m == 0) {
    goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
    sx[i__] = *sa * sx[i__];
/* L30: */
    }
    if (*n < 5) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
    sx[i__] = *sa * sx[i__];
    sx[i__ + 1] = *sa * sx[i__ + 1];
    sx[i__ + 2] = *sa * sx[i__ + 2];
    sx[i__ + 3] = *sa * sx[i__ + 3];
    sx[i__ + 4] = *sa * sx[i__ + 4];
/* L50: */
    }
    return 0;
} /* sscal_ */

integer lsei::isamax_(integer *n, double *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    double r__1;

    /* Local variables */
    double xmag, smax;
    integer i__, ii, ns;

    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n <= 0) {
    return ret_val;
    }
    ret_val = 1;
    if (*n <= 1) {
    return ret_val;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    smax = dabs(sx[1]);
    ns = *n * *incx;
    ii = 1;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    xmag = (r__1 = sx[i__], dabs(r__1));
    if (xmag <= smax) {
        goto L5;
    }
    ret_val = ii;
    smax = xmag;
L5:
    ++ii;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */

L20:
    smax = dabs(sx[1]);
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
    xmag = (r__1 = sx[i__], dabs(r__1));
    if (xmag <= smax) {
        goto L30;
    }
    ret_val = i__;
    smax = xmag;
L30:
    ;
    }
    return ret_val;
} /* isamax_ */

doublereal lsei::sdot_(integer *n, double *sx, integer *incx, double *sy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    double ret_val;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
    return ret_val;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ret_val += sx[ix] * sy[iy];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. */

L20:
    m = *n - *n / 5 * 5 ;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ret_val += sx[i__] * sy[i__];
/* L30: */
    }
    if (*n < 5) {
    return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
    ret_val = ret_val + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] + 
        sx[i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] + sx[
        i__ + 4] * sy[i__ + 4];
/* L50: */
    }
    return ret_val;

/*        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    ret_val += sx[i__] * sy[i__];
/* L70: */
    }
    return ret_val;
} /* sdot_ */

/* Subroutine */ int lsei::saxpy_(integer *n, double *sa, double *sx, integer *incx, 
    double *sy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0 || *sa == 0.) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    sy[iy] += *sa * sx[ix];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4. */

L20:
    m = *n - (*n / 4 << 2);
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    sy[i__] += *sa * sx[i__];
/* L30: */
    }
    if (*n < 4) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
    sy[i__] += *sa * sx[i__];
    sy[i__ + 1] += *sa * sx[i__ + 1];
    sy[i__ + 2] += *sa * sx[i__ + 2];
    sy[i__ + 3] += *sa * sx[i__ + 3];
/* L50: */
    }
    return 0;

/*        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    sy[i__] = *sa * sx[i__] + sy[i__];
/* L70: */
    }
    return 0;
} /* saxpy_ */

/* Subroutine */ int lsei::srotm_(integer *n, double *sx, integer *incx, double *sy, 
    integer *incy, double *sparam)
{
    /* Initialized data */

    double zero = 0.;
    double two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__;
    double w, z__, sflag;
    integer kx, ky, nsteps;
    double sh11, sh12, sh21, sh22;

    /* Parameter adjustments */
    --sparam;
    --sy;
    --sx;

    /* Function Body */

    sflag = sparam[1];
    if (*n <= 0 || sflag + two == zero) {
    goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
    goto L70;
    }

    nsteps = *n * *incx;
    if (sflag < 0.) {
    goto L50;
    } else if (sflag == 0) {
    goto L10;
    } else {
    goto L30;
    }
L10:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    w = sx[i__];
    z__ = sy[i__];
    sx[i__] = w + z__ * sh12;
    sy[i__] = w * sh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
    w = sx[i__];
    z__ = sy[i__];
    sx[i__] = w * sh11 + z__;
    sy[i__] = -w + sh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    w = sx[i__];
    z__ = sy[i__];
    sx[i__] = w * sh11 + z__ * sh12;
    sy[i__] = w * sh21 + z__ * sh22;
/* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0) {
    kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
    ky = (1 - *n) * *incy + 1;
    }

    if (sflag < 0.) {
    goto L120;
    } else if (sflag == 0) {
    goto L80;
    } else {
    goto L100;
    }
L80:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = sx[kx];
    z__ = sy[ky];
    sx[kx] = w + z__ * sh12;
    sy[ky] = w * sh21 + z__;
    kx += *incx;
    ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = sx[kx];
    z__ = sy[ky];
    sx[kx] = w * sh11 + z__;
    sy[ky] = -w + sh22 * z__;
    kx += *incx;
    ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = sx[kx];
    z__ = sy[ky];
    sx[kx] = w * sh11 + z__ * sh12;
    sy[ky] = w * sh21 + z__ * sh22;
    kx += *incx;
    ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* srotm_ */

/* Subroutine */ int lsei::drotmg_(doublereal *dd1, doublereal *dd2, doublereal *
    dx1, doublereal *dy1, doublereal *dparam)
{
    /* Initialized data */

    doublereal zero = 0.;
    doublereal one = 1.;
    doublereal two = 2.;
    doublereal gam = 4096.;
    doublereal gamsq = 16777216.;
    doublereal rgamsq = 5.9604645e-8;

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    doublereal dflag, dtemp, du, dp1, dp2, dq2, dq1, dh11, dh21, dh12, 
        dh22;
    integer igo;

    /* Parameter adjustments */
    --dparam;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT */
    if (! (*dd1 < zero)) {
    goto L10;
    }
/*       GO ZERO-H-D-AND-DX1.. */
    goto L60;
L10:
/*     CASE-DD1-NONNEGATIVE */
    dp2 = *dd2 * *dy1;
    if (! (dp2 == zero)) {
    goto L20;
    }
    dflag = -two;
    goto L260;
/*     REGULAR-CASE.. */
L20:
    dp1 = *dd1 * *dx1;
    dq2 = dp2 * *dy1;
    dq1 = dp1 * *dx1;

    if (! (abs(dq1) > abs(dq2))) {
    goto L40;
    }
    dh21 = -(*dy1) / *dx1;
    dh12 = dp2 / dp1;

    du = one - dh12 * dh21;

    if (! (du <= zero)) {
    goto L30;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L30:
    dflag = zero;
    *dd1 /= du;
    *dd2 /= du;
    *dx1 *= du;
/*         GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (dq2 < zero)) {
    goto L50;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L50:
    dflag = one;
    dh11 = dp1 / dp2;
    dh22 = *dx1 / *dy1;
    du = one + dh11 * dh22;
    dtemp = *dd2 / du;
    *dd2 = *dd1 / du;
    *dd1 = dtemp;
    *dx1 = *dy1 * du;
/*         GO SCALE-CHECK */
    goto L100;
/*     PROCEDURE..ZERO-H-D-AND-DX1.. */
L60:
    dflag = -one;
    dh11 = zero;
    dh12 = zero;
    dh21 = zero;
    dh22 = zero;

    *dd1 = zero;
    *dd2 = zero;
    *dx1 = zero;
/*         RETURN.. */
    goto L220;
/*     PROCEDURE..FIX-H.. */
L70:
    if (! (dflag >= zero)) {
    goto L90;
    }

    if (! (dflag == zero)) {
    goto L80;
    }
    dh11 = one;
    dh22 = one;
    dflag = -one;
    goto L90;
L80:
    dh21 = -one;
    dh12 = one;
    dflag = -one;
L90:
    switch (igo) {
    case 0: goto L120;
    case 1: goto L150;
    case 2: goto L180;
    case 3: goto L210;
    }
/*     PROCEDURE..SCALE-CHECK */
L100:
L110:
    if (! (*dd1 <= rgamsq)) {
    goto L130;
    }
    if (*dd1 == zero) {
    goto L160;
    }
    igo = 0;
//    igo_fmt = fmt_120;
/*              FIX-H.. */
    goto L70;
L120:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 *= d__1 * d__1;
    *dx1 /= gam;
    dh11 /= gam;
    dh12 /= gam;
    goto L110;
L130:
L140:
    if (! (*dd1 >= gamsq)) {
    goto L160;
    }
    igo = 1;
//    igo_fmt = fmt_150;
/*              FIX-H.. */
    goto L70;
L150:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 /= d__1 * d__1;
    *dx1 *= gam;
    dh11 *= gam;
    dh12 *= gam;
    goto L140;
L160:
L170:
    if (! (abs(*dd2) <= rgamsq)) {
    goto L190;
    }
    if (*dd2 == zero) {
    goto L220;
    }
    igo = 2;
//    igo_fmt = fmt_180;
/*              FIX-H.. */
    goto L70;
L180:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 *= d__1 * d__1;
    dh21 /= gam;
    dh22 /= gam;
    goto L170;
L190:
L200:
    if (! (abs(*dd2) >= gamsq)) {
    goto L220;
    }
    igo = 3;
//    igo_fmt = fmt_210;
/*              FIX-H.. */
    goto L70;
L210:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 /= d__1 * d__1;
    dh21 *= gam;
    dh22 *= gam;
    goto L200;
L220:
    if (dflag < 0.) {
    goto L250;
    } else if (dflag == 0) {
    goto L230;
    } else {
    goto L240;
    }
L230:
    dparam[3] = dh21;
    dparam[4] = dh12;
    goto L260;
L240:
    dparam[2] = dh11;
    dparam[5] = dh22;
    goto L260;
L250:
    dparam[2] = dh11;
    dparam[3] = dh21;
    dparam[4] = dh12;
    dparam[5] = dh22;
L260:
    dparam[1] = dflag;
    return 0;
} /* drotmg_ */

/* Subroutine */ int lsei::dcopy_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dy[iy] = dx[ix];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7. */

L20:
    m = *n % 7;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
    dy[i__] = dx[i__];
    dy[i__ + 1] = dx[i__ + 1];
    dy[i__ + 2] = dx[i__ + 2];
    dy[i__ + 3] = dx[i__ + 3];
    dy[i__ + 4] = dx[i__ + 4];
    dy[i__ + 5] = dx[i__ + 5];
    dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;

/*        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dy[i__] = dx[i__];
/* L70: */
    }
    return 0;
} /* dcopy_ */

/* Subroutine */ int lsei::dswap_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m;
    doublereal dtemp1, dtemp2, dtemp3;
    integer ix, iy, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp1 = dx[ix];
    dx[ix] = dy[iy];
    dy[iy] = dtemp1;
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3. */

L20:
    m = *n % 3;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dtemp1 = dx[i__];
    dx[i__] = dy[i__];
    dy[i__] = dtemp1;
/* L30: */
    }
    if (*n < 3) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
    dtemp1 = dx[i__];
    dtemp2 = dx[i__ + 1];
    dtemp3 = dx[i__ + 2];
    dx[i__] = dy[i__];
    dx[i__ + 1] = dy[i__ + 1];
    dx[i__ + 2] = dy[i__ + 2];
    dy[i__] = dtemp1;
    dy[i__ + 1] = dtemp2;
    dy[i__ + 2] = dtemp3;
/* L50: */
    }
    return 0;
L60:

/*     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dtemp1 = dx[i__];
    dx[i__] = dy[i__];
    dy[i__] = dtemp1;
/* L70: */
    }
    return 0;
} /* dswap_ */

doublereal lsei::dnrm2_(integer *n, doublereal *dx, integer *incx)
{
    /* Initialized data */

    doublereal zero = 0.;
    doublereal one = 1.;
    doublereal cutlo = 8.232e-11;
    doublereal cuthi = 1.304e19;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
 //   double sqrt(doublereal);

    /* Local variables */
    doublereal xmax;
    integer next, i__, j, nn;
    doublereal hitest, sum;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT */
    if (*n > 0) {
    goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
//    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 BEGIN MAIN LOOP */
    i__ = 1;
L20:
    switch (next) {
    case 0: goto L30;
    case 1: goto L50;
    case 2: goto L70;
    case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
    goto L85;
    }
    next = 1;
//    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (dx[i__] == zero) {
    goto L200;
    }
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
    goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
    next = 2;
//    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i__ = j;
    next = 3;
//    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
    goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
    if ((d__1 = dx[i__], abs(d__1)) <= xmax) {
    goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;


/*     FOR double OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / (double) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
    if ((d__1 = dx[j], abs(d__1)) >= hitest) {
        goto L100;
    }
/* L95: */
/* Computing 2nd power */
    d__1 = dx[j];
    sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
    goto L20;
    }

/*          TERM. OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */

doublereal lsei::dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    integer i__, m, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
    return ret_val;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    ret_val += (d__1 = dx[i__], abs(d__1));
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6. */

L20:
    m = *n % 6;
    if (m == 0) {
    goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
    ret_val += (d__1 = dx[i__], abs(d__1));
/* L30: */
    }
    if (*n < 6) {
    return ret_val;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
    ret_val = ret_val + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1],
         abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[
        i__ + 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) + (
        d__6 = dx[i__ + 5], abs(d__6));
/* L50: */
    }
    return ret_val;
} /* dasum_ */

/* Subroutine */ int lsei::dscal_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
    return 0;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. */

L20:
    m = *n % 5;
    if (m == 0) {
    goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
    dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
    dx[i__] = *da * dx[i__];
    dx[i__ + 1] = *da * dx[i__ + 1];
    dx[i__ + 2] = *da * dx[i__ + 2];
    dx[i__ + 3] = *da * dx[i__ + 3];
    dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

    integer lsei::idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    doublereal dmax__, xmag;
    integer i__, ii, ns;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n <= 0) {
    return ret_val;
    }
    ret_val = 1;
    if (*n <= 1) {
    return ret_val;
    }
    if (*incx == 1) {
    goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    dmax__ = abs(dx[1]);
    ns = *n * *incx;
    ii = 1;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    xmag = (d__1 = dx[i__], abs(d__1));
    if (xmag <= dmax__) {
        goto L5;
    }
    ret_val = ii;
    dmax__ = xmag;
L5:
    ++ii;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */

L20:
    dmax__ = abs(dx[1]);
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
    xmag = (d__1 = dx[i__], abs(d__1));
    if (xmag <= dmax__) {
        goto L30;
    }
    ret_val = i__;
    dmax__ = xmag;
L30:
    ;
    }
    return ret_val;
} /* idamax_ */

    doublereal lsei::ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy,   integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
    return ret_val;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ret_val += dx[ix] * dy[iy];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1. */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. */

L20:
    m = *n % 5;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    ret_val += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
    return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
    ret_val = ret_val + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
        dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[
        i__ + 4] * dy[i__ + 4];
/* L50: */
    }
    return ret_val;

/*         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    ret_val += dx[i__] * dy[i__];
/* L70: */
    }
    return ret_val;
} /* ddot_ */

/* Subroutine */ int lsei::daxpy_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, m, ix, iy, ns, mp1;

/* ***FIRST EXECUTABLE STATEMENT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0 || *da == 0.) {
    return 0;
    }
    if (*incx == *incy) {
    if ((i__1 = *incx - 1) < 0) {
        goto L5;
    } else if (i__1 == 0) {
        goto L20;
    } else {
        goto L60;
    }
    }
L5:

/*        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS. */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dy[iy] += *da * dx[ix];
    ix += *incx;
    iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4. */

L20:
    m = *n % 4;
    if (m == 0) {
    goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
    dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
    return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
    dy[i__] += *da * dx[i__];
    dy[i__ + 1] += *da * dx[i__ + 1];
    dy[i__ + 2] += *da * dx[i__ + 2];
    dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;

/*        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    dy[i__] = *da * dx[i__] + dy[i__];
/* L70: */
    }
    return 0;
} /* daxpy_ */

/* Subroutine */ int lsei::drotm_(integer *n, doublereal *dx, integer *incx, 
    doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    doublereal zero = 0.;
    doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__;
    doublereal dflag, w, z__;
    integer kx, ky, nsteps;
    doublereal dh11, dh12, dh22, dh21;

    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT */
    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero) {
    goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
    goto L70;
    }

    nsteps = *n * *incx;
    if (dflag < 0.) {
    goto L50;
    } else if (dflag == 0) {
    goto L10;
    } else {
    goto L30;
    }
L10:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    w = dx[i__];
    z__ = dy[i__];
    dx[i__] = w + z__ * dh12;
    dy[i__] = w * dh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
    w = dx[i__];
    z__ = dy[i__];
    dx[i__] = w * dh11 + z__;
    dy[i__] = -w + dh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
    w = dx[i__];
    z__ = dy[i__];
    dx[i__] = w * dh11 + z__ * dh12;
    dy[i__] = w * dh21 + z__ * dh22;
/* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0) {
    kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
    ky = (1 - *n) * *incy + 1;
    }

    if (dflag < 0.) {
    goto L120;
    } else if (dflag == 0) {
    goto L80;
    } else {
    goto L100;
    }
L80:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = dx[kx];
    z__ = dy[ky];
    dx[kx] = w + z__ * dh12;
    dy[ky] = w * dh21 + z__;
    kx += *incx;
    ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = dx[kx];
    z__ = dy[ky];
    dx[kx] = w * dh11 + z__;
    dy[ky] = -w + dh22 * z__;
    kx += *incx;
    ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
    w = dx[kx];
    z__ = dy[ky];
    dx[kx] = w * dh11 + z__ * dh12;
    dy[ky] = w * dh21 + z__ * dh22;
    kx += *incx;
    ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* drotm_ */


/*************************error functions*********************************************/
/* Subroutine */ int lsei::xerror_(const char *messg, integer *nmessg, integer *nerr,
    integer *level)
{
    err_message temp(messg,nerr,level);
    errors.insert(temp);
    return 0;
} /* xerror_ */

/* Subroutine */ int lsei::xerrwv_(const char *messg, integer *nmessg, integer *nerr, 
    integer *level, integer *ni, integer *i1, integer *i2, integer *nr, 
    double *r1, double *r2)
{    /* Format strings */
    err_message temp(messg,nerr,level);
    errors.insert(temp);
    return 0;
} /* xerrwv_ */

/*************************error functions*********************************************/
int lsei::SelfTest()
{
	integer klog[5];
	double cond[3];
	integer istat;
	int stat=4;

	printf("istat=1 NEITHER WNNLS( ) NOR LSEI( ) PASSED THE TEST. \n=2 WNNLS( ) PASSED BUT LSEI( ) FAILED THE TEST. \n=3 LSEI( ) PASSED BUT WNNLS( ) FAILED THE TEST. \n=4 BOTH WNNLS( ) AND LSEI( ) PASSED THE TEST. \n");

    klog[0]=3;    klog[1]=1;    klog[2]=2;    klog[3]=0;    klog[4]=2 ;
	cond[0]=1; cond[1]=cond[2]=1000.;
	clstp_(klog,cond,&istat);
	printf("Stat: %d\n",istat);
	if(stat>istat) stat=istat;

    klog[0]=3;    klog[1]=1;    klog[2]=2;    klog[3]=1;    klog[4]=2 ;
	cond[0]=1; cond[1]=cond[2]=1000.;
	clstp_(klog,cond,&istat);
	printf("Stat: %d\n",istat);
	if(stat>istat) stat=istat;

    klog[0]=3;    klog[1]=1;    klog[2]=2;    klog[3]=2;    klog[4]=2 ;
	cond[0]=1000.; cond[1]=cond[2]=1000.;
	clstp_(klog,cond,&istat);
	printf("Stat: %d\n",istat);
	if(stat>istat) stat=istat;

    klog[0]=4;    klog[1]=2;    klog[2]=2;    klog[3]=2;    klog[4]=3 ;
	cond[0]=1000.; cond[1]=cond[2]=1000.;
	clstp_(klog,cond,&istat);
	printf("Stat: %d\n",istat);
	if(stat>istat) stat=istat;

    klog[0]=5;    klog[1]=3;    klog[2]=2;    klog[3]=2;    klog[4]=4 ;
	cond[0]=1000.; cond[1]=cond[2]=1000.;
	clstp_(klog,cond,&istat);
	printf("Stat: %d\n",istat);
	if(stat>istat) stat=istat;

	return stat;
}

/* Subroutine */ int lsei::clstp_(integer *klog, double *cond, integer *istat)
{    /* Format strings */


    /* System generated locals */
    integer i__1, i__2;
    double r__1;

    /* Builtin functions */
//    double sqrt(doublereal);

    /* Local variables */

    double beta;
    integer mode;
    logical done;
    integer jcol, mepi;
    
    double ansr;
    integer jpnn;    

    double zero;
    integer irow;
    integer i__, j, k, l, n;
    double conda, t;
    integer icase;
    double w[6240]	/* was [96][65] */, conde, x[64], condg;
    integer iseed, mepma;
    double bnorm, gnorm;
    integer iwork[640];
    double dxbyx;    
    integer ka;
    double gg[1024]	/* was [32][32] */;
    integer ma;
    double hh[1024]	/* was [32][32] */;
    integer ke, kg, me;
    double sa[32];
    integer ki, mg, mi;
    double se[32];
    integer kn;
    double sg[32];
    integer nn;
    double rnorme, solerr, srelpr, rnorml, prgopt[4];
    integer np1;
    double gam;
    integer mdg, mdh, mna, mne;
    double phi, one;
    
    integer ngo, mng, mdw;
    double rho, two;
    integer n20100, n20011, n20020, n20130, n20032, n20024, n20016, 
	    n20042, n20050, n20028, n20037, n20046, n20054, n20058, n20063, 
	    n20068, n20073, n20077, n20081, n20085, n20089, n20096, n20104, 
	    n20114, n20118, n20122, n20126;


/*     THE FOLLOWING DIMENSION STATEMENTS NEED NOT BE ALTERED TO */
/*     SOLVE LARGER PROBLEMS. */
    /* Parameter adjustments */
    --cond;
    --klog;

    /* Function Body */
    mdw = 96;
    mdh = 32;
    mdg = 32;
    zero = 0.;
    one = 1.;
    two = 2.;
    *istat = 1;

/*     COMPUTE THE RELATIVE MACHINE PRECISION. */

    M_precision = one;
L10:
    if (one + M_precision == one) 	goto L20;
    M_precision /= two;
    goto L10;
L20:

// fix
	M_precision *= 10;
    srelpr = M_precision;

/*     SET THE OUTPUT UNIT TO WHICH ERROR MESSAGES WILL BE PRINTED. */
/*     SET UP THE PROBLEM DIMENSIONS */
    ka = klog[1];
    ke = klog[2];
    kg = klog[3];
    ki = klog[4];
    kn = klog[5];
    conda = one;
    conde = one;
    condg = one;
    done = kn < 0;
    if (! done) {	goto L30; }
    return 0;

L30:
    ma = 0;
    me = 0;
    mg = 0;
    n = 0;

/*     SET NOISE-TO-SIGNAL RATIO PARAMETER FOR LEAST SQUARES EQUAS. */
/*     THIS ESSENTIALLY SPECIFIES THE RATIO */

/*                   NORM(RESIDUAL VECTOR) */
/*                   --------------------- */
/*                       NORM(A*X) */
    ansr = .01;
/*     SET UP THE CONDITION NUMBERS FOR THE MATRICES A, E, AND G. */
    if (ka >= 0) {
	conda = cond[1];
    }
    if (ke >= 0) {
	conde = cond[2];
    }

    if (kg >= 0) {
	condg = cond[3];
    }

/*     CHECK THE VALIDITY OF THE INPUT */

    if (ka <= 5 && ke <= 5 && kg <= 5 && ki <= 5 && kn <= 5) { goto L40;    }


	xerror_("KA, KE, KG, KI, AND KN MUST ALL BE .LE . 5 AS REQUIRED BY THE CURRENT SUBPROGRAM DIMENSION STATEMENTS.\n",0,0,0);
    return 0;

L40:
    if (! (conda < one || conde < one || condg < one)) {
	goto L50;
    }

    xerror_("CONDA, CONDE, AND CONDG MUST ALL BE .GE. ONE.\n",0,0,0);;

    return 0;

L50:
    icase = 1;
L60:
    iseed = 100001;
    t = ran_(&iseed);
    iseed = 0;

/*     COMPUTE THE PRE-MULTIPLYING HADAMARD MATRIX FOR E. */
    k = ke;
    ngo = 0;
    goto L900;
L70:
    me = nn;

/*     SAVE THE HADAMARD MATRIX. */
    j = 1;
    n20011 = me;
    goto L90;
L80:
    ++j;
L90:
    if (n20011 - j < 0) {
	goto L100;
    }
    scopy_(&me, &hh[(j << 5) - 32], &c__1, &gg[(j << 5) - 32], &c__1);
    goto L80;

/*     NOW FORM THE POST-MULTIPLYING HADAMARD MATRIX. */
L100:
    k = kn;
    ngo = 1;
    goto L900;

L110:
    n = nn;

/*     COMPUTE THE SINGULAR VALUES OF THE MATRIX E. */
/*     DISTRIBUTE THEM UNIFORMLY BETWEEN 1. AND CONDE. */
    se[0] = conde;
/* Computing MAX */
    i__1 = 1, i__2 = min(me,n);
    mne = max(i__1,i__2);
    se[mne - 1] = one;
    i__ = 2;
    n20016 = mne - 1;
    goto L130;
L120:
    ++i__;
L130:
    if (n20016 - i__ < 0) {
	goto L140;
    }
    se[i__ - 1] = one + ran_(&iseed) * (conde - one);
    goto L120;
L140:
    j = 1;
    n20020 = mne;
    goto L160;
L150:
    ++j;
L160:
    if (n20020 - j < 0) {
	goto L170;
    }
    sscal_(&me, &se[j - 1], &gg[(j << 5) - 32], &c__1);
    goto L150;
L170:
    j = 1;
    n20024 = n;
    goto L190;
L180:
    ++j;
L190:
    if (n20024 - j < 0) {
	goto L230;
    }
    if (me > 0) {
	w[j * 96 - 96] = zero;
    }
    scopy_(&me, &w[j * 96 - 96], &c__0, &w[j * 96 - 96], &c__1);
    i__ = 1;
    n20028 = mne;
    goto L210;
L200:
    ++i__;
L210:
    if (n20028 - i__ < 0) {
	goto L220;
    }
    saxpy_(&me, &hh[i__ + (j << 5) - 33], &gg[(i__ << 5) - 32], &c__1, &w[j * 96 - 96], &c__1);
    goto L200;

L220:
    goto L180;

/*     COMPUTE E*X AND STORE IN W(*,N+1). */
L230:
    i__ = 1;
    n20032 = me;
    goto L250;
L240:
    ++i__;
L250:
    if (n20032 - i__ < 0) {
	goto L260;
    }


    x[0] = one;
    w[i__ + (n + 1) * 96 - 97] = sdot_(&n, x, &c__0, &w[i__ - 1], &mdw);
    goto L240;

/*     COMPUTE THE PRE-MULTIPLYING HADAMARD MATRIX FOR A. */
L260:
    k = ka;
    ngo = 2;
    goto L900;

L270:
    ma = nn;

/*     SAVE THE HADAMARD MATRIX. */
    j = 1;
    n20037 = ma;
    goto L290;
L280:
    ++j;
L290:
    if (n20037 - j < 0) {
	goto L300;
    }
    scopy_(&ma, &hh[(j << 5) - 32], &c__1, &gg[(j << 5) - 32], &c__1);
    goto L280;

/*     NOW FORM THE POST-MULTIPLYING HADAMARD MATRIX. */
L300:
    k = kn;
    ngo = 3;
    goto L900;

L310:
    n = nn;

/*     COMPUTE THE SINGULAR VALUES OF THE MATRIX A. */
/*     DISTRUBUTE THEM UNIFORMLY BETWEEN 1. AND CONDA. */
    sa[0] = conda;
/* Computing MAX */
    i__1 = 1, i__2 = min(ma,n);
    mna = max(i__1,i__2);
    sa[mna - 1] = one;
    i__ = 2;
    n20042 = mna - 1;
    goto L330;
L320:
    ++i__;
L330:
    if (n20042 - i__ < 0) {
	goto L340;
    }
    sa[i__ - 1] = one + ran_(&iseed) * (conda - one);
    goto L320;
L340:
    j = 1;
    n20046 = mna;
    goto L360;
L350:
    ++j;
L360:
    if (n20046 - j < 0) {
	goto L370;
    }
    sscal_(&ma, &sa[j - 1], &gg[(j << 5) - 32], &c__1);
    goto L350;
L370:
    j = 1;
    n20050 = n;
    goto L390;
L380:
    ++j;
L390:
    if (n20050 - j < 0) {
	goto L430;
    }

/*     COMPUTE THE PRODUCT IN PLACE INTO W(*,*). */
    if (ma > 0) {	w[me + 1 + j * 96 - 97] = zero;    }
    scopy_(&ma, &w[me + 1 + j * 96 - 97], &c__0, &w[me + 1 + j * 96 - 97], &c__1);
    i__ = 1;
    n20054 = mna;
    goto L410;
L400:
    ++i__;
L410:
    if (n20054 - i__ < 0) {
	goto L420;
    }
    saxpy_(&ma, &hh[i__ + (j << 5) - 33], &gg[(i__ << 5) - 32], &c__1, &w[me  + 1 + j * 96 - 97], &c__1);
    goto L400;
L420:
    goto L380;

/*     COMPUTE A*X AND STORE IN W(*,N+1). */
L430:
    i__ = 1;
    n20058 = ma;
    goto L450;
L440:
    ++i__;
L450:
    if (n20058 - i__ < 0) {
	goto L460;
    }
    mepi = me + i__;
    x[0] = one;
    w[mepi + (n + 1) * 96 - 97] = sdot_(&n, x, &c__0, &w[mepi - 1], &mdw);
    goto L440;
L460:
    bnorm = snrm2_(&ma, &w[me + 1 + (n + 1) * 96 - 97], &c__1);

/*     ADD COMPONENTS TO RIGHT SIDE THAT ARE ORTHOGONAL TO COL. */
/*     SPACE OF A. */
    k = ka;
    ngo = 4;
    goto L900;
L470:
    ma = nn;
    i__ = n + 1;
    n20063 = ma;
    goto L490;
L480:
    ++i__;
L490:
    if (n20063 - i__ < 0) {
	goto L500;
    }
    t = ran_(&iseed) * bnorm * ansr;
    saxpy_(&ma, &t, &hh[(i__ << 5) - 32], &c__1, &w[me + 1 + (n + 1) * 96 -  97], &c__1);
    goto L480;

/*     COMPUTE THE PRE-MULTIPLYING HADAMARD MATRIX FOR G. */
L500:
    k = kg;
    ngo = 5;
    goto L900;
L510:
    mg = nn;

/*     SAVE THE HADAMARD MATRIX. */
    j = 1;
    n20068 = mg;
    goto L530;
L520:
    ++j;
L530:
    if (n20068 - j < 0) {
	goto L540;
    }
    scopy_(&mg, &hh[(j << 5) - 32], &c__1, &gg[(j << 5) - 32], &c__1);
    goto L520;

/*     NOW FORM THE POST-MULTIPLYING HADAMARD MATRIX. */
L540:
    k = kn;
    ngo = 6;
    goto L900;
L550:
    n = nn;

/*     COMPUTE THE SINGULAR VALUES OF G. */
/*     DISTRIBUTE THEM UNIFORMLY BETWEEN 1. AND CONDG. */
    sg[0] = condg;
/* Computing MAX */
    i__1 = 1, i__2 = min(mg,n);
    mng = max(i__1,i__2);
    sg[mng - 1] = one;
    i__ = 2;
    n20073 = mng - 1;
    goto L570;
L560:
    ++i__;
L570:
    if (n20073 - i__ < 0) {
	goto L580;
    }
    sg[i__ - 1] = one + ran_(&iseed) * (condg - one);
    goto L560;
L580:
    j = 1;
    n20077 = mng;
    goto L600;
L590:
    ++j;
L600:
    if (n20077 - j < 0) {
	goto L610;
    }
    sscal_(&mg, &sg[j - 1], &gg[(j << 5) - 32], &c__1);
    goto L590;
L610:
    j = 1;
    n20081 = n;
    goto L630;
L620:
    ++j;
L630:
    if (n20081 - j < 0) {
	goto L670;
    }
    mepma = me + ma;
    if (mg > 0) {
	w[mepma + 1 + j * 96 - 97] = zero;
    }
    scopy_(&mg, &w[mepma + 1 + j * 96 - 97], &c__0, &w[mepma + 1 + j * 96 -   97], &c__1);
    i__ = 1;
    n20085 = mng;
    goto L650;

L640:
    ++i__;
L650:
    if (n20085 - i__ < 0) {
	goto L660;
    }
    saxpy_(&mg, &hh[i__ + (j << 5) - 33], &gg[(i__ << 5) - 32], &c__1, &w[ mepma + 1 + j * 96 - 97], &c__1);
    goto L640;

L660:
    goto L620;

/*     COMPUTE G*X AND STORE IN W(*,N+1). */
L670:
    i__ = 1;
    n20089 = mg;
    goto L690;
L680:
    ++i__;
L690:
    if (n20089 - i__ < 0) {
	goto L700;
    }
    irow = i__ + mepma;
    x[0] = one;
    w[irow + (n + 1) * 96 - 97] = sdot_(&n, x, &c__0, &w[irow - 1], &mdw);
    goto L680;

/*     MAKE FIRST MI=2**KI OF THE INEQUALITIES STRICT. */
L700:
    if (! (ki >= 0)) {
	goto L710;
    }
    mi = 1;
    goto L720;
L710:
    mi = 0;
L720:
    k = 1;
    n20096 = ki;
    goto L740;
L730:
    ++k;
L740:
    if (n20096 - k < 0) {
	goto L750;
    }
    mi += mi;
    goto L730;
L750:
    gnorm = snrm2_(&mg, &w[mepma + 1 + (n + 1) * 96 - 97], &c__1);
    i__ = 1;
    n20100 = min(mi,mg);
    goto L770;
L760:
    ++i__;
L770:
    if (n20100 - i__ < 0) {
	goto L780;
    }
    irow = i__ + mepma;
    w[irow + (n + 1) * 96 - 97] -= ran_(&iseed) * gnorm;
    goto L760;

/*     OBTAIN THE CONSTRAINED LEAST SQUARES SOLUTION. */
/*     NOTE THE LENGTHS OF THE WORK ARRAYS IN IWORK(*). */
L780:
    iwork[0] = 1518;
    iwork[1] = 640;
    if (! (icase == 1)) {
	goto L810;
    }

/*     EXCHANGE POSITIONS OF THE ROWS (A B) AND (G H). */
    np1 = n + 1;
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	irow = me + (ma + mg + 2) / 2;
	i__2 = (ma + mg + 1) / 2;
	sswap_(&i__2, &w[me + 1 + j * 96 - 97], &c__1, &w[irow + j * 96 - 97],&c_n1);

/* L790: */
    }

/*     MOVE RT-SIDE TO W(*,N+MG+1). */
    jcol = n + mg + 1;
    i__1 = me + ma + mg;
    scopy_(&i__1, &w[(n + 1) * 96 - 96], &c__1, &w[jcol * 96 - 96], &c__1);

/*     PUT IN SLACK VARIABLE COLS. AS REQUIRED. */
    if (! (mg > 0)) {
	goto L800;
    }
    w[(n + 1) * 96 - 96] = zero;
    i__1 = mdw * mg;
    scopy_(&i__1, &w[(n + 1) * 96 - 96], &c__0, &w[(n + 1) * 96 - 96], &c__1);
    w[me + 1 + (n + 1) * 96 - 97] = -one;
    i__1 = mdw + 1;
    scopy_(&mg, &w[me + 1 + (n + 1) * 96 - 97], &c__0, &w[me + 1 + (n + 1) * 96 - 97], &i__1);
L800:

/*    SET THE OPTION (NUMBER 6) FOR WNNLS( ) TO SCALE THE */
/*     COLUMNS OF THE MATRIX TO HAVE UNIT LENGTH. */
    prgopt[0] = 4.;
    prgopt[1] = 6.;
    prgopt[2] = 1.;
    prgopt[3] = 1.;
    i__1 = me + mg;
    i__2 = n + mg;
    wnnls_(w, &mdw, &i__1, &ma, &i__2, &n, prgopt, x, &rnorml, &mode/*, iwork,   work*/);
    goto L820;

L810:

/*     SET THE OPTION (NUMBER 2) FOR LSEI( ) TO SCALE THE */
/*     COLUMNS OF THE MATRIX TO HAVE UNIT LENGTH. */
    prgopt[0] = 4.;
    prgopt[1] = 2.;
    prgopt[2] = 1.;
    prgopt[3] = 1.;
    lsei_(w, &mdw, &me, &ma, &mg, &n , prgopt, x, &rnorme, &rnorml, &mode/*, work, iwork*/);
L820:

/*     COMPUTE THE RESIDUAL SUM OF SQUARES OF ERRORS. */
    solerr = zero;
    i__ = 1;
    n20104 = n;
    goto L840;
L830:
    ++i__;
L840:
    if (n20104 - i__ < 0) {
	goto L850;
    }
    x[i__ - 1] = one - x[i__ - 1];
/* Computing 2nd power */

    r__1 = x[i__ - 1];
    solerr += r__1 * r__1;
    goto L830;
L850:
    solerr = sqrt(solerr);

/*     TEST SIZE OF ERROR (REF. LAWSON-HANSON, PAGE 51 AND CH. 16) */
    phi = 100.;
    t = (double) n;
    dxbyx = solerr / sqrt(t);
    rho = one;
    if (bnorm != zero) {
	rho = rnorml / bnorm;
    }

    gam = (double) ((max(ma,n) * 6 - min(ma,n) * 3) * min(ma,n));
    beta = conda * (one + conda * rho) * gam * phi;
    if (! (dxbyx + beta == beta)) {
	goto L860;
    }
    *istat += icase;
    goto L870;
L860:
    if (ma == 0 && mode == 0) {
	*istat += icase;
    }
L870:
    if (! (icase == 1)) {
	goto L880;
    }

	printf("LEAST SQS. RESID. %f  ,FOR WNNLS( )\n ERRORS, 1-X(I), FOR WNNLS( ).\n",rnorml);

    i__1 = n;

    for (i__ = 1; i__ <= i__1; ++i__) {
		printf("%d %g ",i__,x[i__-1]);
    }
	printf("\n");

    goto L890;

L880:

    printf( "LEAST SQS. RESID. %f ,FOR LSEI( ) COMP.RANK OF E %d, COMP.RANK OF REDUCED A %d \nERRORS, 1-X(I), FOR LSEI( ).\n",rnorml,iwork[0],iwork[1]);

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		printf("%d %g ",i__,x[i__-1]);
    }

	printf("\n");
L890:
    if (icase == 2) {
	return 0;
    }
    icase = 2;
    goto L60;

/*     PROCEDURE (GET HADAMARD MATRIX) */
L900:
    nn = 0;
    if (! (k >= 0)) {
	goto L940;
    }
    nn = 1;
    i__ = 1;
    n20114 = k;
    goto L920;
L910:
    ++i__;
L920:
    if (n20114 - i__ < 0) {
	goto L930;
    }
    nn += nn;
    goto L910;
L930:
    goto L950;
L940:
    nn = 0;
    goto L1080;

/*     PLACE THE SYMMETRIC HADAMARD MATRIX IN THE ARRAY HH(*,*). */
L950:
    hh[0] = one;
    nn = 1;
    l = 1;
    n20118 = k;
    goto L970;
L960:
    ++l;
L970:
    if (n20118 - l < 0) {
	goto L1040;
    }
    j = 1;
    n20122 = nn;
    goto L990;
L980:
    ++j;
L990:
    if (n20122 - j < 0) {
	goto L1000;
    }
    scopy_(&nn, &hh[(j << 5) - 32], &c__1, &hh[nn + 1 + (j << 5) - 33], &c__1) ;
    goto L980;
L1000:
   j = 1;
    n20126 = nn;
    goto L1020;
L1010:
    ++j;
L1020:
    if (n20126 - j < 0) {
	goto L1030;
    }
    jpnn = j + nn;
    i__1 = nn << 1;
    scopy_(&i__1, &hh[(j << 5) - 32], &c__1, &hh[(jpnn << 5) - 32], &c__1);
    r__1 = -one;
    sscal_(&nn, &r__1, &hh[nn + 1 + (jpnn << 5) - 33], &c__1);
    goto L1010;
L1030:
    nn += nn;
    goto L960;


/*     MAKE THE MATRIX ORTHOGONAL BY SCALING THE ENTRIES. */
L1040:
    t = (double) nn;
    t = one / sqrt(t);
    j = 1;
    n20130 = nn;
    goto L1060;
L1050:
    ++j;
L1060:
    if (n20130 - j < 0) {
	goto L1070;
    }
    sscal_(&nn, &t, &hh[(j << 5) - 32], &c__1);
    goto L1050;
L1070:
L1080:
    switch (ngo) {
	case 0: goto L70;
	case 1: goto L110;
	case 2: goto L270;
	case 3: goto L310;
	case 4: goto L470;
	case 5: goto L510;
	case 6: goto L550;
    }
 return 0;
} /* clstp_ */


/*global function*/
doublereal lsei::ran_(integer *k)
{
    /* Initialized data */

    integer iy = 100001;

    /* System generated locals */
    double ret_val;

    if (*k > 0) {
    iy = *k;
    }
    iy *= 125;
    iy -= iy / 2796203 * 2796203;
    ret_val = (double) iy / 2796203.;
    return ret_val;
/*     ---------- LAST CARD OF RAN ---------- */
} /* ran_ */
/************************Prints all the errors recorded*********************/
void lsei::printErrors()
{
    errors.printErrors();
}


// err_message1.cpp: implementation of the err_message class.
//
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

err_message::err_message()
{
	//messg = 0;
	nerr = 0;
	level = 0;
}

err_message::err_message(const char *messge, integer *nerror, integer *errlevel)
{
	//messg = (char*)malloc((strlen(messge))*sizeof(char));
	//strcpy(messg, messge);
	
	string *temp = new string(messge);
	messg = *temp;
	delete temp;

	nerr = *nerror;
	level = *errlevel;
}

err_message::~err_message()
{}

const char* err_message::Getmessage()
{
	 return messg.c_str();	 
}

int err_message::GetErrNo()
{
	return nerr;
}

int err_message::GetErrLevel()
{
	return level;
}


err_list::err_list()
{}

err_list::~err_list()
{}

bool err_list::insert(err_message& n)
{
	v.push_back(n);
	return true;
}

void err_list::printErrors()
{
	for(int i = 0; i<v.size(); i++)
	{
		printf("Message: %s, errono: %d, level: %d\n",v[i].Getmessage(),v[i].GetErrNo(),v[i].GetErrLevel());
		//cout << v[i].Getmessage() << ", " << v[i].GetErrNo() << ", " << v[i].GetErrLevel() << endl;
	}
}
