#include "Basic.h"
#include "odysseyRandom.h"
 
#define AsNUM1 1

#pragma once
using namespace std;
 
class mcPlain : public Basic
{
public:
	mcPlain(int cpFlag, double sval,  double xval, 		  
	      int matuN, int *matu,   
		  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  long SimulNo);
/*
(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr);
*/
	virtual ~mcPlain();

	 
	double getValue();
	

	 
private:	

	odysseyRandom *rnd;
	long m_noiters;

	int m_cpFlag;
	double m_sval[AsNUM1]; // 현재가
 
	double m_xval[AsNUM1];  //m_xval[AsNUM1]

		 
	
	int m_matuN;
	int *m_matu;    // 잔존만기
	 
 
	

	int m_irateN;   // 금리 term 개수
	int *m_irateBDay; 
 
	double *m_iRFrate; // m_iRFrate[AsNUM1][m_irateN];

	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate[AsNUM1]; // 연속 배당
	double m_divbval[AsNUM1];
	int m_divN[AsNUM1];
	int m_MAXdivN;
	int **m_divBDay;  //m_divBDay[AsNUM1][m_MAXdivN]
	double **m_divCash; // m_divCash[AsNUM1][m_MAXdivN] 

	double m_divApy;

	int m_voltype; // 볼 type 0:상수, 1;atm term vol, 2: local 
	double m_volbval[AsNUM1];	

	int m_volTN;
	int m_volSN;
	int *m_volBDay;
	double *m_volSmness;
	double ***m_vol;  //m_vol[AsNUM1][m_volSN][m_volTN]; 
	
	
	 
	double mc_price();

	void curttimeVol(double *vol, int day, int Asidx);  // vol[AsNUM1][m_volSN];
	 
	 
	
	int findinteridx(double x, double *Dx, int DxN);
};


