/* 

This program estimates the quasi-steady-sate features of a two-allele / two-allele system in a finite population, allowing for arbitrary selection and reversible mutation.
The goal is to estimate the steady-state distribution of the four genotypes over time. The equilibrium long-term result is a balance between the forces of mutation, selection, and random genetic drift. 

Results are obtained for an array of actual population sizes (N), which are constant in time.

The population experiences sequential episodes of mutation, recombination, selection, and random genetic drit.

Haploidy is assumed.

Haplotype designations: AB = [1][1], Ab = [1][0], aB = [0][1], and ab = [0][0], with fitness schemes 1, 1 - sa, 1 - sb, and 1 + sab.

All mutation and selection processes are treated deterministically, prior to random sampling.

Mutation:
	There are just two types of mutations at each type of site: - to +, ux01 (beneficial); and + to -, ux10 (deleterious), where x = a or b; the rates are allowed to differ among loci.
	
Drift:
	Separate effective population sizes are allowed for the two loci, but if complete linkage is assumed, these should be set equal to each other. 

	The two loci are allowed to have different Ne, which requires application of two drift samplings, the first at Ne = Nea over the full genotype; and the second applying Nex only to be B locus, 
	keeping the previous A frequencies constant. 

	Letting 1 - (1/Neb) = [1 - (1/Nea)] * [1 - (1/Nex)], Nea >= Neb, leads to 

	Nex = Neb * (Nea - 1) / (Nea - Neb) to be used in the 2nd interval of sampling.
	


Running averages of the population features are kept track of: equilibrium frequencies of the four types; average ne at the two loci.

The run starts with an allele frequency distribution fixed for AB, but this can be modified internally.

After a burnin, statistics are then taken at intervals of N/xinc generations, with a total of ngen sampling intervals.

An array of population sizes to run is set internally.

********** NOTE AN ARRAY OF SCALING FACTORS (scalef) CONTAINED INTERNALLY CAN BE UTILIZED TO SCALE UP THE SPEED OF RUNS. IF THE N IS SCALED DOWN BY THE SAME FACTOR AND THE MUTATION RATES UP BY THE SAME FACTOR,
THIS ENABLES THE OVERALL PROGRAM TO RUN FASTER, AS ALL IS A FUNCTION OF THE PRODUCTS N*S AND N*U. ***** NOTE THOUGH THAT KFAC CAN'T BE SO LARGE AS TO MAKE THE SELECTION COEFFICIENTS OR MUTATIONS APPROACH 1.0. ******

THE LOOP IMMEDIATELY BELOW THESE ARRAYS (itera) NEEDS TO BE SET TO DETERMINE THE RANGE OF POPULATION SIZES TO RUN. 

*/



/* ********************************************************************************************************************** */

#define ua		0.00000001					/* mutation rate of A to a */

#define ub		0.00000001					/* mutation rate of B to b */

#define muta		1.0						/* ratio of beneficial to deleterious mutation rates, locus A */

#define mutb		1.0						/* ratio of beneficial to deleterious mutation rates, locus B */

#define scoa		0.000001				/* selection coefficient against aB (positive if deleterious) */

#define scob		0.000001				/* selection coefficient against Ab (positive if deleterious) */

#define scoab		0.000					/* selection coefficient in favor of AB (positive if advantageous) */

#define neratio		1.00					/* ratio of Ne at locus B to that at locus A (array in program is for locus A, so this should be <= 1.0, with B having the smaller Ne, as in an organelle genome) */

#define recr		0.0						/* recombination rate: recr = 0.5 is free recombination. */


#define xinc		10						/* statistics to be recorded every ne/xinc generations */

#define burnin		100000					/* number of initial burn-in sampling increments, ignored in the statistics */

#define tintprint	1000000					/* number of sampling increments between screen printing */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include    <string.h> 


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
char filename[100];

 void main(int argc, char *argv[]) 
{                                                                               
 int f0, f1;                                                                     
                                                                                 
 if (argc > 1)                                                                 
 {                                                                           
     f0 = atoi(argv[1]);                                                             
     f1=f0;                                                                      
     sprintf(filename, "dataout_%d.txt", f0);                                        
 }                                                                               
 else                                                                          
 {                                                                           
     f0 = 1;                                                                     
     f1 = 19;                                                                    
     sprintf(filename, "dataout.txt");                                               
 }  


	/* NOTE: The seed variables can have values between:    0 <= IJ <= 31328 */
	/*                                                      0 <= KL <= 30081 */
	/* Default random seeds for the ranmar() random number generator:        */

	int ij = 1802;
	int kl = 9373;

	/* Uncomment the next line for system clock seeding of ranmar(): */
	ij = time((long)0) % 31329;

	/* Initialization of ranmar() */
	/* i.e., must always call rmarin(ij,kl) before first usage of ranmar */

	rmarin(ij, kl);





	/* ***************************************************************************************** */

	/* MAIN BODY OF PROGRAM. */


	int iga, igb;											/* counters for the classes */

	long igen;												/* generation counter */

	long double tint;												/* interval number for statistics */

	long nea, neb, nex;										/* effective population sizes */

	long double u10a, u01a, u10b, u01b;						/* mutation rates; 01 = beneficial, 10 = deleterious */

	long double selcoa, selcob, selcoab;					/* selection coefficients */

	long kfac;												/* scaling factor for speeding up runs, from scalef[] */

	long efpopn[40];										/* effective population sizes at locus A to run -- NOTE THESE ARE SCALED DOWN ACCORDING TO SCALEF TO INCREASE RUN SPEED */
	long scalef[40];											/* SCALING FACTORS FOR SPEEDING UP RUNS */
    long double rlng[40];

    long double ngens;                                             /* time iterations in run */

	int itera;												/* counter for population-size iterations */

	long increment;											/* increment between surveys */
	long tcount;											/* counter for printing to screen to monitor simulation progress */
	long counter;											/* initial number of surveys to be skipped */

	long double meanfit;									/* mean fitness */

	long double ldiseq;										/* linkage disequilibrium */

	long double wfit[2][2];									/* genotypic fitnesses */

	long double p0[2][2];									/* genotypic frequencies at start of generation */

	long double psel[2][2];									/* after selection */

	long double pmutm[2][2];								/* after mutation */

	long double prec[2][2];									/* after recombination */

	long double pgtypexp[2][2];								/* expected frequencies prior to first epsisode of random drift */

	long double pnew[2][2];									/* frequencies after 2nd sampling episode */

	long double pa1tot, pa0tot;								/* conditional frequencies of A/a alleles */

	long double sumfreqab[2][2];							/* summations of genotype frequencies*/

	long double totw;										/* summations for grand means and variances */

	long double sump;										/* sum of frequencies */

	long double pp;											/* probability associated with the binomial for drift */
	long ntot;												/* integer associated with the binomial for drift */
	long draw;												/* drift drawn from the binomial */
	long double epoi, rnum;									/* terms for Poisson draws */

	long double meanw, grandmeanw;							/* generational mean for fitness */

	long double mean[3][3];									/* mean genotypic frequencies */

	int oldfix, newfix;										/* indicators for fixation states */

	long double fixgens, numfix[5][5], genfix[5][5];		/* counters for numbers and times of fixations of different types */
	
	long double totgens, totratea, totrateb;

	long double meanfixtime[5][5];							/* mean fixation times */
	long double fixrate[5][5];								/* mean fixation rates */

	long double ratea, rateb, neutratea, neutrateb;			/* overall average observed substitution rates and expected netural rates */



	/* Open the output file. */

	remove("dataout.txt ");


	/* Population sizes to run. */

	efpopn[19] = 1000000000;
	efpopn[18] = 465000000;
	efpopn[17] = 216000000;
	efpopn[16] = 100000000;
	efpopn[15] = 46500000;
	efpopn[14] = 21600000;
	efpopn[13] = 10000000;
	efpopn[12] = 4650000;
	efpopn[11] = 2160000;
	efpopn[10] = 1000000;
	efpopn[9] = 465000;
	efpopn[8] = 216000;
	efpopn[7] = 100000;
	efpopn[6] = 46500;
	efpopn[5] = 21600;
	efpopn[4] = 10000;
	efpopn[3] = 4650;
	efpopn[2] = 2160;
	efpopn[1] = 1000;


	/* Scaling factors for N, u, and s to reduce run times. */

	scalef[19] = 10000;
	scalef[18] = 10000;
	scalef[17] = 10000;
	scalef[16] = 10000;
	scalef[15] = 10000;
	scalef[14] = 1000;
	scalef[13] = 1000;
	scalef[12] = 1000;
	scalef[11] = 1000;
	scalef[10] = 100;
	scalef[9] = 100;
	scalef[8] = 100;
	scalef[7] = 10;
	scalef[6] = 10;
	scalef[5] = 10;
	scalef[4] = 1;
	scalef[3] = 1;
	scalef[2] = 1;
	scalef[1] = 1;


	/* Number of sampling increments in run; each increment is (ne/10) generations */

	rlng[19] = 50000000.0;
	rlng[18] = 50000000.0;
	rlng[17] = 100000000.0;
	rlng[16] = 100000000.0;
	rlng[15] = 100000000.0;
	rlng[14] = 200000000.0;
	rlng[13] = 200000000.0;
	rlng[12] = 300000000.0;
	rlng[11] = 300000000.0;
	rlng[10] = 300000000.0;
	rlng[9] = 400000000.0;
	rlng[8] = 400000000.0;
	rlng[7] = 400000000.0;
	rlng[6] = 500000000.0;
	rlng[5] = 500000000.0;
	rlng[4] = 500000000.0;
	rlng[3] = 600000000.0;
	rlng[2] = 600000000.0;
	rlng[1] = 600000000.0;




for (itera = f0; itera <= f1; ++itera){							/* Start iterations over the set of population sizes and mutation rates. */
        
        stream=fopen(filename, "a");

        /* Set the run length. */
        
        ngens = rlng[itera];



		/* Set the initial population-genetic parameters. */

		nea = efpopn[itera];										/* effective population size for locus A */
		neb = int(((double) nea) * neratio);						/* effective population size for locus B */

		nex = 0;
		if (neb != nea) {
			nex = int(((double)neb) * (((double)nea) - 1.0) / (((double)nea) - ((double)neb))); }		/* effective population size to be used for second sampling episode to allow for lower Ne at locus B */
																										/* with the double sampling scheme used for locus B, this renders the desired effective size = neb */

		kfac = scalef[itera];										/* scaling factor for changing ne, mutation rates, and selection coefficients */

		nea = nea / kfac;
		neb = neb / kfac;
		nex = nex / kfac;

		u10a = ((double) kfac) * ua;							/* deleterious and beneficial mutation rates at locus A */
		u01a = muta * u10a;

		u10b = ((double) kfac) * ub;							/* deleterious and beneficial mutation rates at locus B */
		u01b = mutb * u10b;

		selcoa = ((double) kfac) * scoa;						/* selection coefficients */
		selcob = ((double) kfac) * scob;
		selcoab = ((double) kfac) * scoab;

		wfit[1][1] = 1.0;										/* genotypic fitnesses */
		wfit[0][1] = 1.0 - selcoa;
		wfit[1][0] = 1.0 - selcob;
		wfit[0][0] = 1.0 + selcoab;

		p0[1][1] = 0.25;										/* set the initial genotype frequencies to be random. THIS IS ARBITRARY, AND THE HISTORICAL EFFECT ON RUNS IS MINIMAL WITH A LONG BURN-IN PERIOD. */
		p0[1][0] = 0.25;
		p0[0][1] = 0.25;
		p0[0][0] = 0.25;


		/* Initiate the counters for estimating transition rates between states. */

		oldfix = 1;												
		fixgens = 0.0;
		totgens = 0.0;
		totratea = 0.0;
		totrateb = 0.0;

		for (iga = 1; iga <= 4; ++iga) {
			for (igb = 1; igb <= 4; ++igb) {
				numfix[iga][igb] = 0.0;
				genfix[iga][igb] = 0.0;
				meanfixtime[iga][igb] = 0.0;
				fixrate[iga][igb] = 0.0; } }


		/* Initiate the genotype frequencies and counters. */

		for (iga = 0; iga <= 1; ++iga) {
			for (igb = 0; igb <= 1; ++igb) {
				pmutm[iga][igb] = 0.0;									/* zero the various allele-frequency counters */
				psel[iga][igb] = 0.0;
				pgtypexp[iga][igb] = 0.0;
				pnew[iga][igb] = 0.0;
				sumfreqab[iga][igb] = 0.0; 	} }

		igen = 0;
		tcount = 0;
		tint = 0.0;
		counter = 0;
		totw = 0.0;

		increment = nea / xinc;											/* increment in generations between statistic calculations (set as a fraction of nea). */




		/* ******************************************************************************************************************************************* */


		/* Iterate the recursion equations to obtain the equilibrium expectations. */

		while (tint < ngens)  										/* iterate until the stopping criterion has been met. */
		{
			igen = igen + 1;


			/* Impose selection on the genotypic classes. */

			meanfit = 0.0;

			for (iga = 0; iga <= 1; ++iga) {						/* calculate mean relative fitness */
				for (igb = 0; igb <= 1; ++igb) {
					meanfit = meanfit + (p0[iga][igb] * wfit[iga][igb]); } 	}

			for (iga = 0; iga <= 1; ++iga) {						/* genotypic frequencies after selection */
				for (igb = 0; igb <= 1; ++igb) {
					psel[iga][igb] = p0[iga][igb] * wfit[iga][igb] / meanfit; } }




			/* Impose mutation on the post-selection genotypic classes. */

			pmutm[1][1] = (u01a * u01b * psel[0][0]) + (u01a * (1.0 - u10b) * psel[0][1]) + ((1.0 - u10a) * u01b * psel[1][0]) + ((1.0 - u10a) * (1.0 - u10b) * psel[1][1]);

			pmutm[1][0] = (u01a * (1.0 - u01b) * psel[0][0]) + (u01a * u10b * psel[0][1]) + ((1.0 - u10a) * (1.0 - u01b) * psel[1][0]) + ((1.0 - u10a) * u10b * psel[1][1]);

			pmutm[0][1] = ((1.0 - u01a) * u01b * psel[0][0]) + ((1.0 - u01a) * (1.0 - u10b) * psel[0][1]) + (u10a * u01b * psel[1][0]) + (u10a * (1.0 - u10b) * psel[1][1]);

			pmutm[0][0] = ((1.0 - u01a) * (1.0 - u01b) * psel[0][0]) + ((1.0 - u01a) * u10b * psel[0][1]) + (u10a * (1.0 - u01b) * psel[1][0]) + (u10a * u10b * psel[1][1]);




			/* Impose recombination on the post-mutation classes. */

			ldiseq = (pmutm[0][1] * pmutm[1][0]) - (pmutm[1][1] * pmutm[0][0]);

			prec[1][1] = pmutm[1][1] + (recr * ldiseq);
			prec[0][0] = pmutm[0][0] + (recr * ldiseq);
			prec[0][1] = pmutm[0][1] - (recr * ldiseq);
			prec[1][0] = pmutm[1][0] - (recr * ldiseq);



			/* Reset the next generation's expected genotype frequencies, and ensure that they sum to 1.0. */

			sump = prec[1][1] + prec[0][1] + prec[1][0] + prec[0][0];

			for (iga = 0; iga <= 1; ++iga) {
				for (igb = 0; igb <= 1; ++igb) {
					pgtypexp[iga][igb] = prec[iga][igb] / sump; 
					p0[iga][igb] = 0.0;	}}
				
				


			/* Sample the population for new genotype frequencies after the first two-locus sampling episode. */
					/* This allows for sampling at the first level with nea, and is all that is required if nea = neb. */
			/* Uses a Poisson approximation when expected frequencies are extreme. */

			ntot = nea;
			sump = 0.0;

			for (iga = 0; iga <= 1; ++iga) {
				for (igb = 0; igb <= 1; ++igb) {

					if ((pgtypexp[iga][igb] > 0.0) && (ntot > 0))  {
						pp = pgtypexp[iga][igb] / (1.0 - sump);												/* this is the remaining frequency to sample */

						if (pp >= 1.0000000000000) {														/* if remaining frequency = 1.0, then numerator is equal to remaining sample */
							draw = ntot;
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						else if (pp < 0.0000001) {															/* if expected frequency is very small, just draw from a Poisson */
							draw = -1;
							epoi = exp(-pp*((double)ntot));
							rnum = 1.0;
							while (rnum >= epoi) {
								rnum = rnum * ranmar();
								draw = draw + 1;
							}
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						else if (pp > 0.9999999) {															/* if expected frequency is very high, just draw the minor allele from a Poisson */
							draw = -1;
							epoi = exp(-(1.0 - pp)*((double)ntot));
							rnum = 1.0;
							while (rnum >= epoi) {
								rnum = rnum * ranmar();
								draw = draw + 1;
							}
							draw = ntot - draw;
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						else {																				/* for all other frequencies, draw a binomial based on the remaining frequencies */
							draw = ignbin(ntot, pp);
							p0[iga][igb] = ((double)draw) / ((double)nea);
						}

						ntot = ntot - draw;
						sump = sump + pgtypexp[iga][igb];

					}
				}
			}



			/* Sample the population for new genotype frequencies after the second B-locus sampling episode, if neb <> nea, using nex as the effective population size. */

			if (nea != neb) {

				for (iga = 0; iga <= 1; ++iga) {
					for (igb = 0; igb <= 1; ++igb) {
						pnew[iga][igb] = 0.0; } }

				pa1tot = p0[1][1] + p0[1][0];									/* frequency of allele A at first locus */
				pa0tot = p0[0][1] + p0[0][0];									/* frequency of allele a at first locus */

				if (pa1tot > 0.0) {

					pp = p0[1][1] / pa1tot;
					ntot = int(pa1tot * ((double)nex));
					if (ntot < 1){
						ntot = 1; }

					if (pp == 1.0) {
						pnew[1][1] = p0[1][1]; }
					else if (pp == 0.0) {
						pnew[1][0] = p0[1][0]; 	}
					else {
						draw = ignbin(ntot, pp);
						pnew[1][1] = pa1tot * ((double)draw) / ((double)ntot);
						pnew[1][0] = pa1tot - pnew[1][1];	}
				}

				if (pa0tot > 0.0) {

					pp = p0[0][0] / pa0tot;
					ntot = int(pa0tot * ((double)nex));
					if (ntot < 1){
						ntot = 1; }

					if (pp == 1.0) {
						pnew[0][0] = p0[0][0]; }
					else if (pp == 0.0) {
						pnew[0][1] = p0[0][1]; 	}
					else {
						draw = ignbin(ntot, pp);
						pnew[0][0] = pa0tot * ((double)draw) / ((double)ntot);
						pnew[0][1] = pa0tot - pnew[0][0]; 	}
				}

				for (iga = 0; iga <= 1; ++iga) {
					for (igb = 0; igb <= 1; ++igb) {
						p0[iga][igb] = pnew[iga][igb];
						if (p0[iga][igb] < 0.0){
							p0[iga][igb] = 0.0;	}}}
			}



			/* Check for new fixation, and if it occurs, record statistics; cutoff frequency is arbitrarily set as 0.999. */
					/* Legend: AB = 1; Ab = 2; aB = 3; ab = 4. */

			totgens = totgens + 1.0;

			if (p0[1][1] > 0.999){ newfix = 1; }
			else if (p0[1][0] > 0.999){ newfix = 2; }
			else if (p0[0][1] > 0.999){ newfix = 3; }
			else if (p0[0][0] > 0.999){ newfix = 4; }
			else { newfix = 0; }

			fixgens = fixgens + 1.0;

			if ((newfix != oldfix) && (newfix != 0)) {
				numfix[oldfix][newfix] = numfix[oldfix][newfix] + 1.0;
				genfix[oldfix][newfix] = genfix[oldfix][newfix] + fixgens;

				if ((oldfix==1 && newfix==3) || (oldfix==1 && newfix==4) || (oldfix==2 && newfix==3) || (oldfix==2 && newfix==4) || (oldfix==3 && newfix==1) || (oldfix==3 && newfix==2) || (oldfix==4 && newfix==1) || (oldfix==4 && newfix==2)) {
					totratea = totratea + 1.0; }
				else { totrateb = totrateb + 1.0; }
				
				oldfix = newfix;
				fixgens = 0.0; 	
			}



			/* Calculate the summary statistics if the sampling interval is completed. */

			if (igen == increment) {
				igen = 0;
				counter = counter + 1;

				if (counter > burnin) {

					meanw = 0.0;

					for (iga = 0; iga <= 1; ++iga) {
						for (igb = 0; igb <= 1; ++igb) {
							meanw = meanw + (p0[iga][igb] * wfit[iga][igb]);

							sumfreqab[iga][igb] = sumfreqab[iga][igb] + p0[iga][igb];} 	}

					totw = totw + meanw;
					tint = tint + 1.0;
					tcount = tcount + 1;


					if (tcount > tintprint) {

						for (iga = 0; iga <= 1; ++iga) {									/* cumulative mean genotypic frequencies */
							for (igb = 0; igb <= 1; ++igb) {
								mean[iga][igb] = sumfreqab[iga][igb] / tint; } }

						printf("%9d, %9d, %9d, %10.0Lf, %9.5Lf, %9.5Lf, %6.5Lf, %6.5Lf, %6.5Lf, %9.5Lf, %6.5Lf \n", (nea*kfac), (neb*kfac), (nex*kfac), tint, (totw / tint),
							mean[1][1], mean[1][0], mean[0][1], mean[0][0], (mean[1][1] + mean[1][0]), (mean[1][1] + mean[0][1]));

						tcount = 0;
					}
				}
			}								/* ends the summary statistic analysis for this point */
		}									/* ends the loop for generations for this population size. */



		/* Calculate the final statistics. */
		
		grandmeanw = totw / tint;								/* mean fitness */

		for (iga = 1; iga <= 4; ++iga) {										/* mean transtion times from one fixed state to another */
			for (igb = 1; igb <= 4; ++igb) {
				if (numfix[iga][igb] >= 1.0) {
					meanfixtime[iga][igb] = genfix[iga][igb] / numfix[iga][igb]; } }	}

		for (iga = 1; iga <= 4; ++iga) {										/* transition rates */
			for (igb = 1; igb <= 4; ++igb) {
				if (iga != igb) {
					if (meanfixtime[iga][igb] > 0.0) {
						fixrate[iga][igb] = 1.0 / meanfixtime[iga][igb];	} }	} }


		/* Rates of substitution at the two loci. */

		ratea = totratea / totgens;
		rateb = totrateb / totgens;

		neutratea = 2.0 * u10a * u01a / (u10a + u01a);							/* neutral rates of substitution */
		neutrateb = 2.0 * u10b * u01b / (u10b + u01b);



		/* Print the output. */

		fprintf(stream, " %11d, %11d, %11d ,, %12.11f, %12.11f, %4.3f, %4.3f ,, %4.3f,, %12.11f, %12.11f, %12.11f ,,  %6d ,, %17.0Lf, %13d, %17.0Lf ,, %12.11Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf, %12.11Lf, %12.11Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf,, %12.11Lf, %9.1Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf ,, %12.11Lf, %12.11Lf \n  ",
			(nea*kfac), (neb*kfac), (nex*kfac),
			ua, ub, muta, mutb, recr, scoa, scob, scoab,
			kfac, ngens, burnin, (tint*((long double)increment)),
			grandmeanw,
			(mean[1][1] + mean[1][0]), (mean[1][1] + mean[0][1]), 
			mean[1][1], mean[1][0], mean[0][1], mean[0][0],
			fixrate[1][2], numfix[1][2],
			fixrate[1][3], numfix[1][3],
			fixrate[1][4], numfix[1][4],
			fixrate[2][1], numfix[2][1],
			fixrate[2][3], numfix[2][3],
			fixrate[2][4], numfix[2][4],
			fixrate[3][1], numfix[3][1],
			fixrate[3][2], numfix[3][2],
			fixrate[3][4], numfix[3][4],
			fixrate[4][1], numfix[4][1],
			fixrate[4][2], numfix[4][2],
			fixrate[4][3], numfix[4][3],
			ratea, rateb,
			neutratea, neutrateb,
			(ratea / neutratea), (rateb / neutrateb));


		printf("\n");

		fclose(stream);


	}									/* End of the population-size loop. */


exit(0);

}





