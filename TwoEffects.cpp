/* 

This program estimates the features of a multiple-allele locus in a finite population, allowing for exponential or quadratic selection and reversible mutation.
The goal is to estimate the steady-state distribution of alternative allele types over time. The equilibrium long-term result is a balance between the forces of mutation, selection, and random genetic drift. 

***** This version has two types of sites simultaneously present: major with large fitness effects, and minor with small fitness effects.

***** This version allows the use of arbitrary numbers of sites and stength of selection, and also an arbitrary fitness function of the form exp(-s * (n^powerexp)).
      Powerexp = 1.0 is the negative exponential model; Poweexp = 2.0 is the half-Gaussian model. 

The actual population size is assumed to be equal to the effective size (N), which is constant in time.

The population experiences sequential episodes of mutation, selection, and random genetic drit.

Haploidy is assumed, and there is no recombination.

Allele designation:	the general code applies to any trait with a series of sites, each with biallelic states -/+.
	-----...., for the most deleterious allele.
	+++++...., for the most fit allele.

But because there is complete linkage, order does not matter, and a fraction of sites is allocated to the minor vs. major component.

The population is regularly censused after reproduction, at which point the population has size N.

All mutation and selection processes are treated deterministically, prior to random sampling.

Mutation:
	There are just two types of mutations at each type of site: - to +, u01 (beneficial); and + to -, u10 (deleterious); the same rates are assumed at both loci.
	
Selection:
	Fitness is determined by a function that has to be set internally. 
	
Running averages of the population features are kept track of: probabilities of being fixed in the alternative monomorphic states, and the allele-frequency distribution conditional on being polymorphic.

The run starts with an allele frequency distribution equal to the expectation under the sequential model.

After a burnin, statistics are then taken at intervals of Ne/xinc generations, with a total of ngen sampling intervals.

The rates of mutation and strength of selection are defined by input parameters entered immediately below.

NOTE THAT NEAR THE TOP OF THE PROGRAM, A SCALING FACTOR (scalef) CAN BE UTILIZED TO SCALE UP THE SELECTION COEFFICIENTS. IF THE NE IS SCALED DOWN BY THE SAME FACTOR AND THE MUTATION RATES UP BY THE SAME FACTOR,
THIS ENABLES THE OVERALL PROGRAM TO RUN FASTER, AS ALL IS A FUNCTION OF THE PRODUCTS NE*S AND NE*U.

A SERIES OF 21 COMBINATIONS OF POPULATION SIZE (efpopn) AND MUTATION RATE (delmutrate) (USING SCALING FROM LYNCH ET AL. 2016), COVERING THE FULL RANGE OF NATURAL POPULATION VALUES, IS GIVEN INTERNALLY. THE LOOP IMMEDIATELY
BELOW THIS (itera) NEEDS TO BE SET TO DETERMINE THE RANGE OF POPULATION SIZES TO RUN. 

*/



/* ********************************************************************************************************************** */

#define mutbeta		0.1					/* ratio of beneficial to deleterious mutation rates */

#define ellmaj		100					/* total number of sites for the trait with major effects */

#define ellmin		1000					/* total number of sites for the trait with minor effects */

#define scomaj		0.00001					/* selection coefficient for major-effect alleles */

#define scomin		0.0000001				/* selection coefficient for minor-effect alleles */

#define powerexp	1.0						/* 1.0 means an exponential function; 2.0 means a half-Gaussian. */

#define xinc		50						/* statistics to be recorded every ne/xinc generations */

#define ngens		4000000					/* number of sampling increments in run; each increment is (ne/10) generations */

#define burnin		200000					/* number of initial burn-in sampling increments, ignored in the statistics */

#define tintprint	50000					/* number of sampling increments between screen printing */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>




/* ********************************************************************************************************** */


extern __int64 ignbin(__int64 n,double pp);


/* Definitions for the binomial generator. */

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


/* Definitions for the random number generator. */

#define TRUE -1
#define FALSE 0
#define boolean int


static double u[98], c, cd, cm;
static __int64 i97, j97;
static boolean test = FALSE;


double ranmar(void)
/*
C This is the random number generator proposed by George Marsaglia in 
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.
*/
{
        double uni;
        
        if (test==FALSE) {
                puts("Call the init routine rmarin() before calling ranmar().");
                exit(2);
        }
	uni = u[i97] - u[j97];
	if (uni < 0.0) uni += 1.0;
	u[i97] = uni;
	i97--;
	if (i97==0) i97 = 97;
	j97--;
	if (j97==0) j97 = 97;
	c -= cd;
	if (c<0.0) c += cm;
	uni -= c;
	if (uni<0.0) uni += 1.0;
	return uni;
}




/* Seed for the random number generator. */

void rmarin(int ij,int kl) 
{
/*
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient 
C length to complete an entire calculation with. For example, if sveral 
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random 
C number generator can create 900 million different subsequences -- with 
C each subsequence having a length of approximately 10^30.
C 
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0
*/
        int i, j, k, l, ii, jj, m;
        double s, t;
        
        if (ij<0 || ij>31328 || kl<0 || kl>30081) {
                puts("The first random number seed must have a value between 0 and 31328.");
                puts("The second seed must have a value between 0 and 30081.");
                exit(1);
        }
        
        i = (ij/177)%177 + 2;
        j = ij%177 + 2;
        k = (kl/169)%178 + 1;
        l = kl%169;
        
        for (ii=1; ii<=97; ii++) {
                s = 0.0;
                t = 0.5;
                for (jj=1; jj<=24; jj++) {
                        m = (((i*j)%179)*k) % 179;
                        i = j;
                        j = k;
                        k = m;
                        l = (53*l + 1) % 169;
                        if ((l*m)%64 >= 32) s += t;
                        t *= 0.5;
                }
                u[ii] = s;
        }
        
        c = 362436.0 / 16777216.0;
        cd = 7654321.0 / 16777216.0;
        cm = 16777213.0 / 16777216.0;
        
        i97 = 97;
        j97 = 33;
        
        test = TRUE;
}




/* ****************************************************************************************** */

/* binomial random number generator */

__int64 ignbin(__int64 n,double pp)
/*
**********************************************************************
     long ignbin(long n,double pp)
                    GENerate BINomial random deviate
                              Function
     Generates a single random deviate from a binomial
     distribution whose number of trials is N and whose
     probability of an event in each trial is P.
                              Arguments
     n  --> The number of trials in the binomial distribution
            from which a random deviate is to be generated.
     p  --> The probability of an event in each trial of the
            binomial distribution from which a random deviate
            is to be generated.
     ignbin <-- A random deviate yielding the number of events
                from N independent trials, each of which has
                a probability of event P.
                              Method
     This is algorithm BTPE from:
         Kachitvichyanukul, V. and Schmeiser, B. W.
         Binomial Random Variate Generation.
         Communications of the ACM, 31, 2
         (February, 1988) 216.
**********************************************************************
     SUBROUTINE BTPEC(N,PP,ISEED,JX)
     BINOMIAL RANDOM VARIATE GENERATOR
     MEAN .LT. 30 -- INVERSE CDF
       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
       USABLE ALGORITHM.
     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
       "BINOMIAL RANDOM VARIATE GENERATION,"
       COMMUNICATIONS OF THE ACM, FORTHCOMING
     WRITTEN:  SEPTEMBER 1980.
       LAST REVISED:  MAY 1985, JULY 1987
     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
                           GENERATOR
     ARGUMENTS
       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
     VARIABLES
       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
       M:  INTEGER VALUE OF THE CURRENT MODE
       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
       P1:  AREA OF THE TRIANGLE
       C:  HEIGHT OF THE PARALLELOGRAMS
       XM:  CENTER OF THE TRIANGLE
       XL:  LEFT END OF THE TRIANGLE
       XR:  RIGHT END OF THE TRIANGLE
       AL:  TEMPORARY VARIABLE
       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
       P2:  AREA OF THE PARALLELOGRAMS
       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
           FROM THE REGION
       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
           REJECT THE CANDIDATE VALUE
       IX:  INTEGER CANDIDATE VALUE
       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
       K:  ABSOLUTE VALUE OF (IX-M)
       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
           ALSO USED IN THE INVERSE TRANSFORMATION
       R: THE RATIO P/Q
       G: CONSTANT USED IN CALCULATION OF PROBABILITY
       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
            OF F WHEN IX IS GREATER THAN M
       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
             CALCULATION OF F WHEN IX IS LESS THAN M
       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
       YNORM: LOGARITHM OF NORMAL BOUND
       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
       USED IN THE FINAL ACCEPT/REJECT TEST
       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
     REMARK
       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
       ARE NOT INVOLVED.
     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/

{
static float psave = -1.0;
static __int64 nsave = -1;


static __int64 ignbin,i,ix,ix1,k,m,mp,T1;
static double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;



    if(pp != psave) goto S10;
    if(n != nsave) goto S20;
    if(xnp < 30.0) goto S150;
    goto S30;
S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
*/
    psave = pp;
    p = min(psave,1.0-psave);
    q = 1.0-p;
S20:
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = ffm;
    fm = m;
    xnpq = xnp*q;
    p1 = (__int64) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = ranmar()*p4;
    v = ranmar();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = xm-p1*v+u;
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = xl+log(v)/xll;
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = xr-log(v)/xlr;
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = powl(q,double(n));
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = ranmar();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    ignbin = ix;
    return ignbin;
}



/* point to the output file */

FILE *stream;


void main()
{




/* NOTE: The seed variables can have values between:    0 <= IJ <= 31328 */
/*                                                      0 <= KL <= 30081 */
/* Default random seeds for the ranmar() random number generator:        */

int ij = 1802;
int kl = 9373;

/* Uncomment the next line for system clock seeding of ranmar(): */
 ij = time ((long) 0) % 31329; 

/* Initialization of ranmar() */
/* i.e., must always call rmarin(ij,kl) before first usage of ranmar */

rmarin(ij,kl);





/* ***************************************************************************************** */

/* MAIN BODY OF PROGRAM. */


int ig, igmaj, igmin;									/* counters for the classes, with 0 meaning no + alleles, and ell meaning all + alleles */

int igen;												/* generation counter */

long tint;												/* interval number for statistics */

long double difmin, difmaj;								/* deviation between number of + alleles and the optimum, for minor- and major-effect sites */

long double meanfit;									/* mean fitness */

long double gmaxrate;									/* maximum growth rate under additive model */

long double phi;										/* ratio of beneficial and deleterious mutation rates per site */

long double selcomaj, selcomin;							/* selection coefficients associated with major- and minor-effect loci */

long double mut0min, mutdown0min, mutemin, mutupemin;	/* parameters for mutational transition frequencies between classes */
long double mut0maj, mutdown0maj, mutemaj, mutupemaj;

long efpopn[40];										/* effective population sizes to run -- NOTE THESE ARE SCALED DOWN ACCORDING TO SCALEF TO INCREASE RUN SPEED */
long double delmutrate[40];								/* deleterious mutation rates to run -- NOTE THESE ARE SCALED UP ACCORDING TO SCALEF TO INCREASE RUN SPEED */
int scalef[40];											/* SCALING FACTORS FOR SPEEDING UP RUNS */

int itera;												/* counter for population-size/mutation-rate iterations */

long ne;												/* effective population size, from efpopn[] */
long double u10, u01;									/* deleterious and beneficial mutation rates, from delmutrate[] */
int kfac;												/* scaling factor for speeding up runs, from scalef[] */

int newlowmin, newhighmin, newlowmaj, newhighmaj;		/* running settings for upper and lower genoytpic states for minor and major allele counts */
int lowigmin, highigmin, lowigmaj, highigmaj;
int lowmutmin, highmutmin, lowmutmaj, highmutmaj;

long double sump;										/* sum of frequencies */

int startcheckmin, startcheckmaj;						/* checkers for new low and high frequency classes */

long double pp;											/* probability associated with the binomial for drift */
long ntot;												/* integer associated with the binomial for drift */
long draw;												/* drift drawn from the binomial */
long double epoi, rnum;									/* terms for Poisson draws */

long increment;											/* increment between surveys */
long tcount;											/* counter for printing to screen to monitor simulation progress */
long counter;											/* initial number of surveys to be skipped */

long double meanw, meanigmin, meanigmaj;				/* generational means for fitness and allelic class */

long double totw, totwsq;								/* summations for grand means and variances */
long double totigmin, totigsqmin, totwvarigmin;			/* summations for grand means and variances */
long double totigmaj, totigsqmaj, totwvarigmaj;

long double grandmeanw, varw;							/* mean and variance of fitness */

long double grandmeanigmin, varigmin;					/* mean and variances of numbers of minor and major alleles*/
long double grandmeanigmaj, varigmaj;

long double ssqigmin, wvarigmin;						/* sum of squares and within-generation variance of trait */
long double ssqigmaj, wvarigmaj;

long double meanperform;								/* mean relative performance under additive model */

long double exterm, pseqmin, pseqmaj;					/* expected frequencies under the sequential model (assumes no interference) */

long double totmin0, meanminhet, neuthet;				/* statistics for computing the long-term mean heterozygosity at a single neutral minor site (ellmin must be set to 1) */
long double neestimate;


long double *mutransupmin = (long double*)malloc((ellmin+2) * sizeof(long double));
long double *mutransdownmin = (long double*)malloc((ellmin+2) * sizeof(long double));
long double *mutrans0min = (long double*)malloc((ellmin+2) * sizeof(long double));

long double *mutransupmaj = (long double*)malloc((ellmaj + 2) * sizeof(long double));
long double *mutransdownmaj = (long double*)malloc((ellmaj + 2) * sizeof(long double));
long double *mutrans0maj = (long double*)malloc((ellmaj + 2) * sizeof(long double));

long double *sumfreqmaj = (long double*)malloc((ellmaj + 2) * sizeof(long double)); 
long double *sumfreqmin = (long double*)malloc((ellmin + 2) * sizeof(long double)); 
long double *meanpfreqmaj = (long double*)malloc((ellmaj + 2) * sizeof(long double)); 
long double *meanpfreqmin = (long double*)malloc((ellmin + 2) * sizeof(long double));


double **wfit = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ wfit[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **relw = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ relw[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **p0 = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ p0[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pmutm1 = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pmutm1[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pmutm = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pmutm[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **psel = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ psel[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pgtypexp = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pgtypexp[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }





/* Open the output file. */

remove("dataout.txt ");


efpopn[21] = 205636404;
efpopn[20] = 157180983;
efpopn[19] = 117462706;
efpopn[18] = 85822271;
efpopn[17] = 61305579;
efpopn[16] = 42815399;
efpopn[15] = 29234791;
efpopn[14] = 19516413;
efpopn[13] = 12737963;
efpopn[12] = 8128305;
efpopn[11] = 5071075;
efpopn[10] = 3093000;
efpopn[9] = 1844591;
efpopn[8] = 1075474;
efpopn[7] = 613056;
efpopn[6] = 341665;
efpopn[5] = 186166;
efpopn[4] = 99174;
efpopn[3] = 51654;
efpopn[2] = 26300;
efpopn[1] = 13095;

delmutrate[21] = 0.0000000005290;
delmutrate[20] = 0.0000000006488;
delmutrate[19] = 0.0000000008096;
delmutrate[18] = 0.000000001028;
delmutrate[17] = 0.000000001327;
delmutrate[16] = 0.000000001743;
delmutrate[15] = 0.000000002330;
delmutrate[14] = 0.000000003167;
delmutrate[13] = 0.000000004380;
delmutrate[12] = 0.000000006163;
delmutrate[11] = 0.000000008821;
delmutrate[10] = 0.00000001284;
delmutrate[9] = 0.00000001900;
delmutrate[8] = 0.00000002870;
delmutrate[7] = 0.00000004390;
delmutrate[6] = 0.00000006850;
delmutrate[5] = 0.0000001090;
delmutrate[4] = 0.0000001750;
delmutrate[3] = 0.0000002880;
delmutrate[2] = 0.0000004810;
delmutrate[1] = 0.0000008170; 




scalef[21] = 10000;
scalef[20] = 10000;
scalef[19] = 10000;
scalef[18] = 10000;
scalef[17] = 10000;
scalef[16] = 1000;
scalef[15] = 1000;
scalef[14] = 1000;
scalef[13] = 1000;
scalef[12] = 1000;
scalef[11] = 1000;
scalef[10] = 1000;
scalef[9] = 1000;
scalef[8] = 100;
scalef[7] = 100;
scalef[6] = 100;
scalef[5] = 100;
scalef[4] = 10;
scalef[3] = 10;
scalef[2] = 10;
scalef[1] = 10;

gmaxrate = (((double)ellmin) * scomin) + (((double)ellmaj) * scomaj);



for (itera = 1; itera <= 21; ++itera) {							/* Start iterations over the set of population sizes and mutation rates. */
		
	stream=fopen("dataout.txt", "w");
	

	ne = efpopn[itera];											/* effective population size */

	u10 = delmutrate[itera];									/* deleterious and beneficial mutation rates */
	u01 = u10 * mutbeta;

	kfac = scalef[itera];										/* scaling factor for changing ne, mutation rates, and selection coefficients */

	ne = ne / kfac;
	u10 = ((long double)kfac) * u10;
	u01 = ((long double)kfac) * u01;




	/* Set the fitnesses for the allelic classes. */
	/* The first (0,0) and final (ellmin,ellmaj) indices are the worst and best classes. */

	selcomaj = ((long double)kfac) * scomaj;
	selcomin = ((long double)kfac) * scomin;


	exterm = 2.0 * ((double)ne) * selcomaj;
	exterm = mutbeta * exp(exterm);

	pseqmaj = exterm / (1.0 + exterm);
	
	exterm = 2.0 * ((double)ne) * selcomin;
	exterm = mutbeta * exp(exterm);

	pseqmin = exterm / (1.0 + exterm);




	for (igmin = 0; igmin <= ellmin; ++igmin) {							/* genotypic fitnesses */

		sumfreqmin[igmin] = 0.0;

		for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {

			sumfreqmaj[igmaj] = 0.0;

			difmin = ((long double)ellmin) - ((long double)igmin);		/* number of mismatches at loci with minor effects */
			difmaj = ((long double)ellmaj) - ((long double)igmaj);		/* number of mismatches at loci with major effects */

			wfit[igmin][igmaj] = exp(-selcomin * powl(difmin, powerexp)) * exp(-selcomaj * powl(difmaj, powerexp));
			relw[igmin][igmaj] = wfit[igmin][igmaj];

		}
	}





	/* Set the mutation constants. */

	phi = u01 / u10;											/* ratio of beneficial and deleterious mutation rates per site */

	mut0min = 1.0 - (((long double)ellmin) * u01);				/* fraction remaining in minor class 0 */
	mutdown0min = u10;											/* fraction of minor class 1 degrading to minor class 0 */
	mutemin = 1.0 - (((long double)ellmin) * u10);				/* fraction remaining in minor class ellmin */
	mutupemin = u01;											/* fraction of minor class (ellmin-1) moving to minor class ellmin */

	mut0maj = 1.0 - (((long double)ellmaj) * u01);				/* fraction remaining in major class 0 */
	mutdown0maj = u10;											/* fraction of major class 1 degrading to major class 0 */
	mutemaj = 1.0 - (((long double)ellmaj) * u10);				/* fraction remaining in major class ellmaj */
	mutupemaj = u01;											/* fraction of major class (ellmaj-1) moving to major class ellmaj */

	for (igmin = 1; igmin <= (ellmin - 1); ++igmin) {
		mutransupmin[igmin] = u01 * (((long double)ellmin) - ((long double)(igmin - 1))); 											/* gain of a + from any of the - in next lowest class */
		mutransdownmin[igmin] = u10 * ((long double)(igmin + 1));																	/* loss of a + from next highest class */
		mutrans0min[igmin] = 1.0 - (u10 * ((long double)igmin)) - (u01 * (((long double)ellmin) - ((long double)igmin))); 	}		/* stays unchanged */

	for (igmaj = 1; igmaj <= (ellmaj - 1); ++igmaj) {
		mutransupmaj[igmaj] = u01 * (((long double)ellmaj) - ((long double)(igmaj - 1))); 											/* gain of a + from next lowest class */
		mutransdownmaj[igmaj] = u10 * ((long double)(igmaj + 1));																	/* loss of a + from next highest class */
		mutrans0maj[igmaj] = 1.0 - (u10 * ((long double)igmaj)) - (u01 * (((long double)ellmaj) - ((long double)igmaj))); 	}		/* stays unchanged */




	/* Set the initial genotype frequencies to be fixation at midpoint. */

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
			p0[igmin][igmaj] = 0.0; } }


	/*
	p0[ellmin / 2][ellmaj / 2] = 1.0;*/
	
	p0[0][ellmaj] = 1.0;





	/* Initiate the allele frequencies, counters, and test statistics. */

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
			pmutm[igmin][igmaj] = 0.0;						/* zero the various allele-frequency counters */
			psel[igmin][igmaj] = 0.0;
			pgtypexp[igmin][igmaj] = 0.0;
			pmutm1[igmin][igmaj] = 0.0;  } }
	
	igen = 0;
	tcount = 0;
	tint = 0;
	counter = 0;
	totw = 0.0;
	totwsq = 0.0;

	meanminhet = 0.0; 
	
	totigmin = 0.0;
	totigsqmin = 0.0;
	totwvarigmin = 0.0;

	totigmaj = 0.0;
	totigsqmaj = 0.0;
	totwvarigmaj = 0.0;
	
	newhighmin = ellmin;
	newlowmin = 0;
	newhighmaj = ellmaj;
	newlowmaj = 0;

	increment = ne / xinc;						/* increment in generations between statistic calculations (set as a fraction of Ne). */



	/* ******************************************************************************************************************************************* */


	/* Iterate the recursion equations to obtain the equilibrium expectations. */

	while (tint < ngens)  										/* iterate until the stopping criterion has been met. */
	{
		igen = igen + 1;



		/* Set the running upper and lower boundaries to the allelic count classes. */

		
		lowigmin = newlowmin - 1;								
		highigmin = newhighmin + 1;

		if (lowigmin < 0) { lowigmin = 0; }
		if (highigmin > ellmin) { highigmin = ellmin; }
		if (newlowmin < 0) { newlowmin = 0; }
		if (newhighmin > ellmin) { newhighmin = ellmin; }

		lowigmaj = newlowmaj - 1;
		highigmaj = newhighmaj + 1;

		if (lowigmaj < 0) { lowigmaj = 0; }
		if (highigmaj > ellmaj) { highigmaj = ellmaj; }
		if (newlowmaj < 0) { newlowmaj = 0; }
		if (newhighmaj > ellmaj) { newhighmaj = ellmaj; }

		

		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {		/* zero the frequencies of the classes */
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				psel[igmin][igmaj] = 0.0; } }



		/* Impose selection. */

		meanfit = 0.0;

		for (igmin = newlowmin; igmin <= newhighmin; ++igmin) {									/* calculate mean fitness */
			for (igmaj = newlowmaj; igmaj <= newhighmaj; ++igmaj) {
				meanfit = meanfit + (p0[igmin][igmaj] * relw[igmin][igmaj]); } }

		for (igmin = newlowmin; igmin <= newhighmin; ++igmin) {									/* weight the prior genotype frequencies by relative fitness */
			for (igmaj = newlowmaj; igmaj <= newhighmaj; ++igmaj) {
				psel[igmin][igmaj] = p0[igmin][igmaj] * relw[igmin][igmaj] / meanfit; } }





		/* Impose mutation on the genotypic classes. */

		sump = 0.0;

		if (lowigmin == 0) { lowmutmin = 1; }
		else { lowmutmin = lowigmin; }
		if (highigmin == ellmin) { highmutmin = ellmin - 1; }
		else { highmutmin = highigmin; }

		if (lowigmaj == 0) { lowmutmaj = 1; }
		else { lowmutmaj = lowigmaj; }
		if (highigmaj == ellmaj) { highmutmaj = ellmaj - 1; }
		else { highmutmaj = highigmaj; }


		for (igmin = lowmutmin; igmin <= highmutmin; ++igmin) {									/* first, do the minor-allele classes */
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[igmin][igmaj] = (mutransupmin[igmin] * psel[igmin - 1][igmaj]) + (mutransdownmin[igmin] * psel[igmin + 1][igmaj]) + (mutrans0min[igmin] * psel[igmin][igmaj]); }}

		if (lowigmin == 0) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[0][igmaj] = (mut0min * psel[0][igmaj]) + (mutdown0min * psel[1][igmaj]);	}}

		if (highigmin == ellmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[ellmin][igmaj] = (mutemin * psel[ellmin][igmaj]) + (mutupemin * psel[ellmin - 1][igmaj]); }}




		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {									/* next, do the major-allele classes */
			for (igmaj = lowmutmaj; igmaj <= highmutmaj; ++igmaj) {
				pmutm[igmin][igmaj] = (mutransupmaj[igmaj] * pmutm1[igmin][igmaj - 1]) + (mutransdownmaj[igmaj] * pmutm1[igmin][igmaj + 1]) + (mutrans0maj[igmaj] * pmutm1[igmin][igmaj]);
				sump = sump + pmutm[igmin][igmaj]; 	}}

		if (lowigmaj == 0) {
			for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
				pmutm[igmin][0] = (mut0maj * pmutm1[igmin][0]) + (mutdown0maj * pmutm1[igmin][1]);
				sump = sump + pmutm[igmin][0]; }}

		if (highigmaj == ellmaj) {
			for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
				pmutm[igmin][ellmaj] = (mutemaj * pmutm1[igmin][ellmaj]) + (mutupemaj * pmutm1[igmin][ellmaj - 1]);
				sump = sump + pmutm[igmin][ellmaj]; }}



		/* Reset the next generation's expected genotype frequencies, and ensure that they sum to 1.0. */

		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pgtypexp[igmin][igmaj] = pmutm[igmin][igmaj] / sump; }}

		for (igmin = 0; igmin <= ellmin; ++igmin) {
			for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
				p0[igmin][igmaj] = 0.0; }}

		


		/* Sample the population for new genotype frequencies. */
		/* Uses a Poisson approximation when expected frequencies are extreme. */

		ntot = ne;
		sump = 0.0;
		
		newlowmin = ellmin;
		newhighmin = 0;
		newlowmaj = ellmaj;
		newhighmaj = 0;
		
		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {

			if ((pgtypexp[igmin][igmaj] > 0.0) && (ntot > 0))  {
				pp = pgtypexp[igmin][igmaj] / (1.0 - sump);											/* this is the remaining frequency to sample */

				if (pp >= 1.0000000000000) {														/* if remaining frequency = 1.0, then numerator is equal to remaining sample */
					draw = ntot;
					p0[igmin][igmaj] = ((long double)draw) / ((long double)ne); 	}

				else if (pp < 0.0000001) {															/* if expected frequency is very small, just draw from a Poisson */
					draw = -1;
					epoi = exp(-pp*((long double)ntot));
					rnum = 1.0;
					while (rnum >= epoi) {
						rnum = rnum * ranmar();
						draw = draw + 1; }
					p0[igmin][igmaj] = ((long double)draw) / ((long double)ne); }

				else if (pp > 0.9999999) {															/* if expected frequency is very high, just draw the minor allele from a Poisson */
					draw = -1;
					epoi = exp(-(1.0 - pp)*((long double)ntot));
					rnum = 1.0;
					while (rnum >= epoi) {
						rnum = rnum * ranmar();
						draw = draw + 1; }
					draw = ntot - draw;
					p0[igmin][igmaj] = ((long double)draw) / ((long double)ne); }

				else {																				/* for all other frequencies, draw a binomial based on the remaining frequencies */
					draw = ignbin(ntot, pp);
					p0[igmin][igmaj] = ((long double)draw) / ((long double)ne); }

				ntot = ntot - draw;
				sump = sump + pgtypexp[igmin][igmaj];

				if (p0[igmin][igmaj] > 0.0) {
					if (igmin < newlowmin) { newlowmin = igmin; }
					if (igmin > newhighmin) { newhighmin = igmin; }
					if (igmaj < newlowmaj) { newlowmaj = igmaj; }
					if (igmaj > newhighmaj) { newhighmaj = igmaj; } 	}
								
		}}}


		

		/* Calculate the summary statistics if the sampling interval is completed. */

		if (igen == increment) {
			igen = 0;
			counter = counter + 1;


			if (counter > burnin) {
				meanw = 0.0;
				meanigmin = 0.0;
				ssqigmin = 0.0;
				meanigmaj = 0.0;
				ssqigmaj = 0.0;


				for (igmin = 0; igmin <= ellmin; ++igmin) {
					for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
						
						sumfreqmin[igmin] = sumfreqmin[igmin] + p0[igmin][igmaj];
						sumfreqmaj[igmaj] = sumfreqmaj[igmaj] + p0[igmin][igmaj];

						meanw = meanw + (p0[igmin][igmaj] * wfit[igmin][igmaj]);

						meanigmin = meanigmin + (p0[igmin][igmaj] * ((long double)igmin));
						ssqigmin = ssqigmin + (p0[igmin][igmaj] * powl(((long double)igmin), 2.0));

						meanigmaj = meanigmaj + (p0[igmin][igmaj] * ((long double)igmaj));
						ssqigmaj = ssqigmaj + (p0[igmin][igmaj] * powl(((long double)igmaj), 2.0)); } }
				
				wvarigmin = (ssqigmin - (meanigmin * meanigmin)) ;
				wvarigmaj = (ssqigmaj - (meanigmaj * meanigmaj)) ;

				totw = totw + meanw;
				totwsq = totwsq + powl(meanw, 2.0);
				
				totigmin = totigmin + meanigmin;
				totigsqmin = totigsqmin + powl(meanigmin, 2.0);
				totwvarigmin = totwvarigmin + wvarigmin;

				totigmaj = totigmaj + meanigmaj;
				totigsqmaj = totigsqmaj + powl(meanigmaj, 2.0);
				totwvarigmaj = totwvarigmaj + wvarigmaj;

				tint = tint + 1;
				tcount = tcount + 1;

				grandmeanigmin = totigmin / ((long double)tint);
				varigmin = totwvarigmin / ((long double)tint);

				grandmeanigmaj = totigmaj / ((long double)tint);
				varigmaj = totwvarigmaj / ((long double)tint);



				totmin0 = 0.0;											/* Calculates the heterozygosity at the minor site. */

				for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {				/* Is only relevant if there is one such site, and it is set to be neutral. */
					totmin0 = totmin0 + p0[0][igmaj]; 	}

				meanminhet = meanminhet + (2.0 * totmin0 * (1.0 - totmin0));


				
				if (tcount > tintprint) {
					  printf("%9d, %9d, %7.6f, %10.6Lf, %7.6Lf, %6d, %4d, %4d, %6d, %4d, %4d, %10.4Lf, %7.4Lf\n", (ne*kfac), tint, (totw / ((long double)tint)), (grandmeanigmin / ((double)ellmin)), (grandmeanigmaj / ((double)ellmaj)),
						newlowmin, newhighmin, (newhighmin - newlowmin), newlowmaj, newhighmaj, (newhighmaj - newlowmaj), powl(varigmin, 0.5), powl(varigmaj, 0.5));
					
					tcount = 0; }

			}
		}								/* ends the summary statistic analysis for this point */

	}									/* ends the loop for generations. */



	/* Calculate the final statistics. */

	for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
		meanpfreqmaj[ig] = sumfreqmaj[ig] / ((long double)tint); }

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		meanpfreqmin[ig] = sumfreqmin[ig] / ((long double)tint); }

	grandmeanw = totw / ((long double)tint);

	grandmeanigmin = totigmin / ((long double)tint);
	grandmeanigmaj = totigmaj / ((long double)tint);

	varw = (totwsq / ((long double)tint)) - powl(grandmeanw, 2.0);
	varigmin = totwvarigmin / ((long double)tint);
	varigmaj = totwvarigmaj / ((long double)tint);


	meanperform = ((((double)ellmin) * scomin * grandmeanigmin / ((long double)ellmin)) + (((double)ellmaj) * scomaj * grandmeanigmaj / ((long double)ellmaj))) / gmaxrate;

	meanminhet = meanminhet / ((long double)tint);

	neuthet = 4.0 * ne * u01 * u10 / (u01 + u10);

	neestimate = scalef[itera] * meanminhet / (4.0 * u01 * u10 / (u01 + u10));


	fprintf(stream, " %11d, %11d, %11d,, %12.11Lf, %12.11Lf, %3.2Lf,, %12.11Lf, %12.11Lf ,,  %9d ,, %11d, %11d, %17.0Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf,, %12.11Lf, %12.11Lf ,, %12d\n  ",
		(ne*kfac), ellmin, ellmaj,
		(selcomin / ((double)kfac)), (selcomaj / ((double)kfac)), powerexp,
		(u10 / ((double)kfac)), mutbeta, kfac, 
		ngens, burnin, (((double)tint)*((double)increment)), grandmeanw, varw, 
		(grandmeanigmin / ((long double)ellmin)), (grandmeanigmaj / ((long double)ellmaj)), pseqmin, pseqmaj, powl(varigmin, 0.5), powl(varigmaj, 0.5), meanperform, meanminhet, neuthet, ((long)neestimate));

	printf("\n");

	fclose(stream);

}									/* End the set of iterations over all population sizes and mutation rates. */


exit(0);

}





