#include "Basic.h"
#include "odysseyRandom.h"
 
#define AsNUM1 1
 

#pragma once
using namespace std;
 
class TomA : public Basic
{
public:
	TomA(double sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,
	      int matuN, int *matu, double *prate, double ctime,		   
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
	virtual ~TomA();

	 
	double getValue();
	double getBatch(double tsval);
	double getDeltaGammaTheta();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);
	//double getCorrDelta(int corrdeltaN,int *corrdeltaidx); 

	double m_Ogreeks[MAXGREEKS1];	
	double m_Otvegas[AsNUM1*MAXVEGAS];	 
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	 

	int m_modeltype;
	  
	 
	
	//int m_iRateSeq[AsNUM1+1];    // <<-QT
 
	//int m_iRateSeqCD;
	 
private:	

	odysseyRandom *rnd;
	long m_noiters;
	 
 
	 
	 
	 
	double m_sval[AsNUM1]; // 현재가
	double m_stom;
	double m_btom;

	int m_STDayN;
	int *m_STStartDay;
	double *m_IXT;
	double *m_IndexT;
	int *m_STEndDay;
	double *m_IXTheta;
		 
	


	int m_matuN;
	int *m_matu;    // 잔존만기
	double *m_prate;

	 
 
	

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
	
	 
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	 
	double mc_price();

	void curttimeVol(double *vol, int day, int Asidx);  // vol[AsNUM1][m_volSN];
	 
	 
	
	int findinteridx(double x, double *Dx, int DxN);
};


