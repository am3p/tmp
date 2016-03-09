#include "Basic.h"

#define MAXAST 3 // 최대기초자산개수
//#define MAXDAYS 1261 //최대 randomnumber 생성 날짜 개수 252 * 5 + 1
#define MAXDAYS 2000

#define SOBOL_MAX_DIMENSION 40
#define SOBOL_BIT_COUNT 30




//랜덤넘버 생성에 필요한 상수값들
const long RAND_PRIME = 2147483647;
const long RAND_OPT = 16807;
const long IQ = 127773;
const int IR = 2836;
const int DEF_SEEDBASE = 895553;
//static int seedbase;



static const double a[] =
{
        -3.969683028665376e+01,
         2.209460984245205e+02,
        -2.759285104469687e+02,
         1.383577518672690e+02,
        -3.066479806614716e+01,
         2.506628277459239e+00
};

static const double b[] =
{
        -5.447609879822406e+01,
         1.615858368580409e+02,
        -1.556989798598866e+02,
         6.680131188771972e+01,
        -1.328068155288572e+01
};

static const double c[] =
{
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
         4.374664141464968e+00,
         2.938163982698783e+00
};

static const double d[] =
{
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00
};

/* primitive polynomials in binary encoding for sobol */
static const int primitive_polynomials[SOBOL_MAX_DIMENSION] =
{
  1,     3,   7,  11,  13,  19,  25,  37,  59,  47,
  61,   55,  41,  67,  97,  91, 109, 103, 115, 131,
  193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
  213, 191, 253, 203, 211, 239, 247, 285, 369, 299
};

/* degrees of the primitive polynomials for sobol */
static const int degree_table[SOBOL_MAX_DIMENSION] =
{
  0, 1, 2, 3, 3, 4, 4, 5, 5, 5,
  5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 8, 8, 8
};

/*
  initial values for direction tables, following
  Bratley+Fox, taken from [Sobol+Levitan, preprint 1976]
*/   // int -> short
static const int v_init[8][SOBOL_MAX_DIMENSION] =
{
  {
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  },
  {
    0, 0, 1, 3, 1, 3, 1, 3, 3, 1,
    3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
    1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
    3, 1, 1, 3, 1, 3, 1, 3, 1, 3
  },
  {
    0, 0, 0, 7, 5, 1, 3, 3, 7, 5,
    5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
    5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
    5, 1, 1, 5, 7, 7, 5, 1, 3, 3
  },
  {
    0,  0,  0,  0,  0,  1,  7,  9, 13, 11,
    1,  3,  7,  9,  5, 13, 13, 11,  3, 15,
    5,  3, 15,  7,  9, 13,  9,  1, 11,  7,
    5, 15,  1, 15, 11,  5,  3,  1,  7,  9
  },
  {
     0,  0,  0,  0,  0,  0,  0,  9,  3, 27,
    15, 29, 21, 23, 19, 11, 25,  7, 13, 17,
     1, 25, 29,  3, 31, 11,  5, 23, 27, 19,
    21,  5,  1, 17, 13,  7, 15,  9, 31,  9
  },
  {
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0, 37, 33,  7,  5, 11, 39, 63,
    27, 17, 15, 23, 29,  3, 21, 13, 31, 25,
     9, 49, 33, 19, 29, 11, 19, 27, 15, 25
  },
  {
     0,   0,  0,  0,  0,  0,    0,  0,  0,   0,
     0,   0,  0,  0,  0,  0,    0,  0,  0,  13,
    33, 115, 41, 79, 17,  29, 119, 75, 73, 105,
     7,  59, 65, 21,  3, 113,  61, 89, 45, 107
  },
  {
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
    0, 0, 0, 0, 0, 0, 0, 7, 23, 39
  }
};

typedef struct
{
  unsigned int  sequence_count;
  double        last_denominator_inv;
  int           last_numerator_vec[SOBOL_MAX_DIMENSION];
  int           v_direction[SOBOL_BIT_COUNT][SOBOL_MAX_DIMENSION];
} sobol_state_t;




#pragma once
using namespace std;

	 
class odysseyRandom : public Basic
{
private:
	static bool instanceFlag;
	static odysseyRandom *singleRnd;

	int randseqcurstat;
	long randseqmaxx;
	long randseqmaxy;
	long randseqmaxz;
	long **randomxy;

	int randnumcurstat;
	long randnummaxx;
	long randnummaxz;
	double **randomzx;

	int seedbase;
public:	 
	odysseyRandom();
	odysseyRandom(int seed);
	
	long SimulNo; 
	static odysseyRandom* getInstance();
	static odysseyRandom* getInstance(int seed);
	static void removeInstance();
		
	double unirand(double from, double to);
	double unirand(int seed);
	double unirand();
	void setSeed(int seed);
	void setDefaultSeed();
	int randseq(int method,long x,long y,long *seqs);
	double ltqnorm(double p);
	int sobol_get(void * state, unsigned int dimension, double * v);
	int sobol_init(void * state, unsigned int dimension);
	int randnum(int method,long x, long z,double **nums,int randtype);

	virtual ~odysseyRandom();
 
};


