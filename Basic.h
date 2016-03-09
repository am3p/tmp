#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <string>

#pragma once


#define MAX(a,b) (((a)>(b))? (a):(b))
#define MIN(a,b) (((a)>(b))? (b):(a))
#define Heviside(a,b) (((a)>=(b))? (1.0):(0.0))
 
#define Admesh 0.02
#define AdmeshR 0.001
#define myeps 0.000001  

#define kCodeN 100
#define minVol 0.1
#define FXminVol 0.01
#define PI 4.0*atan(1.0)
#define YFactor 365.0 

#define AMOUNT (1)
#define MAXGREEKS1 20
#define MAXGREEKS2 30
#define MAXGREEKS3 30 // 16
#define MAXVEGAS 153
#define MAXRHOS 10
#define MAXCORRDELTA 10

#define GreeksN1 12
#define BasicGRsN1 7 //0:price, 1:udc, 2:ddc, 3:ugc, 4:dgc, 5:t, 6:charmc
#define VegaGRsN1 3 // 7:vega, 8:vanna, 9:zomma
#define RhoGRsN1 2 // 10:rf rho, 11: dc rho

#define GreeksN2 22
#define BasicGRsN2 13  
#define VegaGRsN2 6 
#define RhoGRsN2 2  
#define CorrDeltaGRsN2 1 // corrlation delta

#define Priceidx 0

#define UpDeltaCidx 1
#define DownDeltaCidx 2

#define UpGammaCidx 3
#define DownGammaCidx 4

#define Thetaidx 5
#define CharmCidx 6

#define Vegaidx 7
#define VannaCidx 8
#define ZommaCidx 9

#define Rhorfidx 10
#define Rhodcidx 11

#define SoutRange 20
#define RhoPertubation 0.001
#define VegaPertubation 0.01
#define CorrDeltaPertubation 0.1


#define twoPriceidx 0

#define twoUpDeltaCidx1 1
#define twoDownDeltaCidx1 2
#define twoUpDeltaCidx2 3
#define twoDownDeltaCidx2 4

#define twoUpGammaCidx1 5
#define twoDownGammaCidx1 6
#define twoUpGammaCidx2 7
#define twoDownGammaCidx2 8

#define twoCrossGammaCidx 9

#define twoThetaidx 10
#define twoCharmCidx1 11
#define twoCharmCidx2 12

#define twoVegaidx1 13 
#define twoVegaidx2 14
#define twoVannaCidx1 15
#define twoVannaCidx2 16
#define twoZommaCidx1 17
#define twoZommaCidx2 18

#define twoRhorfidx 19
#define twoRhodcidx 20

#define twoCorrDeltaidx 21
 

// 3S

#define GreeksN3 16
#define BasicGRsN3 8 // p, d1, d2, d3, g1, g2, g3, t 
#define VegaGRsN3 3 // v1, v2, v3
#define RhoGRsN3 5 // r1, r2, r3, rdc, rcd

#define S3Priceidx 0

#define S3Delta1idx 1
#define S3Delta2idx 2
#define S3Delta3idx 3
#define S3Gamma1idx 4
#define S3Gamma2idx 5
#define S3Gamma3idx 6

#define S3Thetaidx 7

#define S3Vega1idx 8
#define S3Vega2idx 9
#define S3Vega3idx 10

#define S3Rhorf1idx 11
#define S3Rhorf2idx 12
#define S3Rhorf3idx 13
#define S3Rhodcidx 14
#define S3RhoCDidx 15



#define maxhgdelta 0.05
#define minhgdelta -0.05





class Basic
{
public:
	Basic();
	virtual ~Basic();
 

protected:
	double cnorm(double x);
	double interp1(double x, double *Dx, double *Dy, int N, int meshType, int interpType);
	double interp1(int x, int *Dx, double *Dy, int N, int meshType, int interpType);
	//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
	int IsinArr(int x, int N, int *Arr); // Arr array 안에 x가 있는지? 있을 때 index + 1 return, 없을 땐 0 return!!
	double interp2(double x, double y, double *Dx, double *Dy, double **Dz, int Nx, int Ny);
	double interp2(int x, int y, int *Dx, int *Dy, double **Dz, int Nx, int Ny);

};
