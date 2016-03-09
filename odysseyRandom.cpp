#include "odysseyRandom.h"
#include <math.h>

#define LOW 0.02425
#define HIGH 0.97575

bool odysseyRandom::instanceFlag = false;
odysseyRandom* odysseyRandom::singleRnd = NULL;

odysseyRandom::odysseyRandom()
{
	odysseyRandom::setDefaultSeed();

	randseqcurstat = 0;
	randseqmaxx    = 0;
	randseqmaxy    = 0;
	randseqmaxz    = 0;
	randnumcurstat = 0;
	randnummaxx    = 0;
	randnummaxz    = 0;
	 
}

odysseyRandom::odysseyRandom(int seed)
{
	odysseyRandom::setSeed(seed);
	 
}

odysseyRandom::~odysseyRandom()
{

}

odysseyRandom* odysseyRandom::getInstance() {
	if(!instanceFlag) {
		singleRnd = new odysseyRandom();
		instanceFlag = true;
		return singleRnd;
	}
	else {
		return singleRnd;
	}
}

odysseyRandom* odysseyRandom::getInstance(int seed) {
	if(!instanceFlag) {
		singleRnd = new odysseyRandom(seed);
		instanceFlag = true;
		return singleRnd;
	}
	else {
		return singleRnd;
	}
}

void odysseyRandom::removeInstance() {
    if(instanceFlag) {
        delete singleRnd;
		instanceFlag = false;		 
	}
}


double odysseyRandom::unirand(int seed)
{
	double rand0;
	int k = seed / IQ;

	seed = RAND_OPT * (seed - k * IQ) - IR * k;

	if(seed < 0) seed += RAND_PRIME;
	seedbase = seed;
	rand0 = (double)seed / (double)RAND_PRIME;

	return rand0;
}

void odysseyRandom::setDefaultSeed()
{
	odysseyRandom::setSeed(DEF_SEEDBASE);
	instanceFlag = false;
}

void odysseyRandom::setSeed(int seed)
{
	seedbase = seed;
	instanceFlag = false;
}

// default seed값을 기초로 random number 생성
double odysseyRandom::unirand()
{
    return unirand(seedbase);
}
 
double odysseyRandom::unirand(double from, double to){
    return ((to - from) * unirand() + from);
}


/*
    Random Number가 들어있는  배열위치를 Random하게 생성하고 가져온다
    X    : Number of Iterations    
    Y    : Length(Time Step) of Random Sequence배열의 크기(보통 만기일수)
    SEQS : Random Number Sequence가 들어갈 배열
    return : 0 (Success), 1 (Failure)
   */

int odysseyRandom::randseq(int method,long x,long y,long *seqs)
{
   /*
    Random Number Sequence
    X    : Number of Iterations    
    Y    : Length(Time Step) of Random Seq 
    Z    : Number of Dimension
    SEQS : Random Number Sequence
    return : 0 (Success), 1 (Failure)
   */
   
   //static int curstat=0;
   //static long maxx=0, maxy=0, maxz=0;
   //static long **randomxy;
   long *randomy;
   long randomstmp;
   long i, j,  pseq;

   /* Method : 0 : Create New Data Sets, Size (x*y*z)     */
   /* Method : 1 : Read Data from Generated Data Sets     */
   /* Method : 2 : Destroy Data Sets                      */
 
   if(method==1) {

      if(randseqcurstat!=1) return(-100);
	  
      /* Method 1 : Effective Parameters
         X : Iteration Number
         Y : Vector Size
         Z : Dimension Number
         seqs : Returning Vector
      */
      /* Checking Validity */
      if((x>=randseqmaxx) || (x<0)) return(-2);
	  
      if((y>randseqmaxy)||(y<1)) return(-3);
	  
      /* for(i=0;i<y;i++) seqs[i]=randomxy[(randomxz[x][z])][i]; */
      memcpy(seqs,randomxy[x],((size_t)((y)*sizeof(long)))); 
	  
      return(0);
   }
   if(method==0) {
      /* Method 0 : Ignore seqs, seqsize */
      if((x<1)||(y<1)) return(-8);
      randomxy=(long **)malloc((size_t)((x)*sizeof(long *)));
      randomy=(long *)malloc((size_t)((x)*sizeof(long)));
      if((!randomxy)||(!randomy)) return(-9);
      for(i=0;i<x;i++) {
         randomxy[i]=(long *)malloc((size_t)((y)*sizeof(long)));
         if((!randomxy[i])) return(-10);
         randomy[i]=i;
      }
      setDefaultSeed();
    for(i=0;i<x;i++) randomxy[i][0]=randomy[i];
      for(j=1;j<y;j++) {
         for(i=(x-1);i>=0;i--) {
            pseq=(long)floor(i*unirand(0,1));
            randomstmp=randomy[i]; randomy[i]=randomy[pseq]; randomy[pseq]=randomstmp;
            randomxy[i][j]=randomy[i];
         }
      }
      free(randomy);
      randseqmaxx=x; randseqmaxy=y; 
      randseqcurstat=1;
   }
   if(method==2) {
      if(randseqcurstat!=1) return(-11);
      for(i=0;i<randseqmaxx;i++) { free(randomxy[i]); };
      free(randomxy); 
      randseqmaxx=randseqmaxy=0;
      randseqcurstat=0;
   }
   return(0);
}



double odysseyRandom::ltqnorm(double p)
{
/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */
   double q, r;
   errno = 0;

   if(p<0||p>1) {
      errno = EDOM;
      return 0.0;
   } else if(p==0) {
      errno = ERANGE;
      return -HUGE_VAL /* minus "infinity" */;
   } else if(p==1) {
      errno = ERANGE;
      return HUGE_VAL /* "infinity" */;
   } else if(p<LOW) {
      /* Rational approximation for lower region */
      q = sqrt(-2*log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
   } else if(p>HIGH) {
      /* Rational approximation for upper region */
      q  = sqrt(-2*log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
   } else {
      /* Rational approximation for central region */
      q = p - 0.5;
      r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
   }
}



/** 
   sobol 난수 발생을 위한 초기화
*/
int odysseyRandom::sobol_init(void * state, unsigned int dimension)
{
  sobol_state_t * s_state = (sobol_state_t *) state;
  unsigned int i_dim;
  int j, k;
  int ell;

  if(dimension < 1 || dimension > SOBOL_MAX_DIMENSION) {
    return(1);
  }

  /* Initialize direction table in dimension 0. */
  for(k=0; k<SOBOL_BIT_COUNT; k++) s_state->v_direction[k][0] = 1;

  /* Initialize in remaining dimensions. */
  for(i_dim=1; i_dim<dimension; i_dim++) {

    const int poly_index = i_dim;
    const int degree_i = degree_table[poly_index];
    int includ[8];

    /* Expand the polynomial bit pattern to separate
     * components of the logical array includ[].
     */
    int p_i = primitive_polynomials[poly_index];
    for(k = degree_i-1; k >= 0; k--) {
      includ[k] = ((p_i % 2) == 1);
      p_i /= 2;
    }

    /* Leading elements for dimension i come from v_init[][]. */
    for(j=0; j<degree_i; j++) s_state->v_direction[j][i_dim] = v_init[j][i_dim];

    /* Calculate remaining elements for this dimension,
     * as explained in Bratley+Fox, section 2.
     */
    for(j=degree_i; j<SOBOL_BIT_COUNT; j++) {
      int newv = s_state->v_direction[j-degree_i][i_dim];
      ell = 1;
      for(k=0; k<degree_i; k++) {
        ell *= 2;
        if(includ[k]) newv ^= (ell * s_state->v_direction[j-k-1][i_dim]);
      }
      s_state->v_direction[j][i_dim] = newv;
    }
  }

  /* Multiply columns of v by appropriate power of 2. */
  ell = 1;
  for(j=SOBOL_BIT_COUNT-1-1; j>=0; j--) {
    ell *= 2;
    for(i_dim=0; i_dim<dimension; i_dim++) {
      s_state->v_direction[j][i_dim] *= ell;
    }
  }

  /* 1/(common denominator of the elements in v_direction) */
  s_state->last_denominator_inv = 1.0 /(2.0 * ell);

  /* final setup */
  s_state->sequence_count = 0;
  for(i_dim=0; i_dim<dimension; i_dim++) s_state->last_numerator_vec[i_dim] = 0;

  return(0);
}

//static int optRandom::sobol_get(void * state, unsigned int dimension, double * v)

/**
  내용 : sobol 난수 발생 함수
*/
int odysseyRandom::sobol_get(void * state, unsigned int dimension, double * v)
{
  sobol_state_t * s_state = (sobol_state_t *) state;

  unsigned int i_dimension;

  /* Find the position of the least-significant zero in count. */
  int ell = 0;
  int c = s_state->sequence_count;
  while(1) {
    ++ell;
    if((c % 2) == 1) c /= 2;
    else break;
  }

  /* Check for exhaustion. */
  if(ell > SOBOL_BIT_COUNT) return(1); /* FIXME: good return code here
*/

  for(i_dimension=0; i_dimension<dimension; i_dimension++) {
    const int direction_i     = s_state->v_direction[ell-1][i_dimension];
    const int old_numerator_i = s_state->last_numerator_vec[i_dimension];
    const int new_numerator_i = old_numerator_i ^ direction_i;
    s_state->last_numerator_vec[i_dimension] = new_numerator_i;
    v[i_dimension] = new_numerator_i * s_state->last_denominator_inv;
  }

  s_state->sequence_count++;

  return(0);
}


/**
   난수를 생성 및 return  함수
*/

int odysseyRandom::randnum(int method, long x, long z, double **nums, int randtype)
{
/*
    Random Number Sequence
    X    : Number of Iterations
    Z    : Number of Dimension
    SEQS : Random Number Sequence
    return : 0 (Success), 1 (Failure)
   */
   sobol_state_t s_state;
   
   //static int curstat=0;
   //static long maxx=0, maxz=0;
   //static double **randomzx;
   double *randomz, *randomz1;
   long i, j;
   int seednum;

   /* Method : 0 : Create New Data Sets, Size (z*x)       */
   /* Method : 1 : Read Data from Generated Data Sets     */
   /* Method : 2 : Destroy Data Sets                      */

   /* RandType 
      0 : Switch to "7"
      
      7 : (ND) Sobol Sequence + Inverse
       
   */
   if(method==1) {
	  	  	

      if(randnumcurstat!=1) return(1);
      /* Method 1 : Effective Parameters
         X    : Iteration Number
         Z    : Dimension Number
         nums : Returning Matrix
      */
      /* Checking Validity */

      if((x>randnummaxx)||(z>randnummaxz)||(x<0)||(z<0)) return(1);

      for(i=0;i<z;i++) {
         memcpy(nums[i],randomzx[i],((size_t)(x*sizeof(double))));
      }
	  


      return(0);
   }
   if(method==0) {
      setDefaultSeed();

	  
      if((x<1)||(z<1)) return(1);

      randomzx=(double **)malloc((size_t)((z)*sizeof(double *)));
      if((!randomzx)) return(1);

      for(i=0;i<z;i++) {
         randomzx[i]=(double *)malloc((size_t)((x)*sizeof(double)));
         if((!randomzx[i])) return(1);
      }

      randomz=(double *)malloc((size_t)((z)*sizeof(double)));
      randomz1=(double *)malloc((size_t)((z)*sizeof(double)));
      seednum=-1;

      if(randtype==0) randtype=7;

      /* Method 0 : Ignore seqs, seqsize */
   
      if(randtype==7) {
         sobol_init(&s_state,z);
         for(i=0;i<x;i++) {
            sobol_get(&s_state,z,randomz);
            for(j=0;j<z;j++){
				randomzx[j][i]=ltqnorm(randomz[j]);
			//	printf("randomzx[%d][%d]=%lf\n",j,i,randomzx[j][i]);
			}
         }
      }
  
      free(randomz); free(randomz1);
      randnummaxx=x; randnummaxz=z;
      randnumcurstat=1;
  
	  
   }
   if(method==2) {
      if(randnumcurstat!=1) return(1);
      for(i=0;i<randnummaxz;i++) free(randomzx[i]);
      free(randomzx);
      randnummaxx=randnummaxz=0;
      randnumcurstat=0;
   }
   return(0);
}
