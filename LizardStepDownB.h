#include "Basic.h"
#include "odysseyRandom.h"
 
#define AsNUM 2

#pragma once
using namespace std;
 
class LizardStepDownB : public Basic
{
public:
	LizardStepDownB(int kiFlag, int koFlag, double *sval, double *bval, double *xval, double *kihval, double *kohval, double *Cpn,
	      int matuN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int corrtype, int corrN, int *corrBDay, double *corr, long SimulNo);

	virtual ~LizardStepDownB();

	 
	double getValue();
	double getDeltaGammaTheta();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);
	//double getCorrDelta(int corrdeltaN,int *corrdeltaidx); 

	double m_Ogreeks[MAXGREEKS2];	
	double m_Otvegas1[MAXVEGAS];
 
	double m_Otvegas2[MAXVEGAS];
 
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	double m_Otcorrdelta[MAXCORRDELTA];

	int m_modeltype;
	 
private:	
	
	odysseyRandom *rnd;
	long m_noiters;


	double m_sval[AsNUM];
	double m_bval[AsNUM];
	double **m_xval;

	double m_kihval[AsNUM];
	

	int m_kiFlag;

	  
	double *m_Cpn; // Cpn & dummy
	 
	
	int m_matuN;	

	int *m_matu;    // 잔존만기
	int m_lizardmatu;

	double m_ctime; // 계산 호출시간 종가전 10, 종가후 16 정도?

	int m_maxevalN;
	int *m_evalN;

	double **m_psval;
	
	int m_irateN;   // 금리 term 개수
	int *m_irateBDay; 
 
	double *m_iRFrate;
	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate[AsNUM]; // 연속 배당
	double m_divbval[AsNUM];
	int m_divN[AsNUM];
	int m_MAXdivN;
	int **m_divBDay;  //m_divBDay[AsNUM][m_MAXdivN]
	double **m_divCash; // m_divCash[AsNUM][m_MAXdivN] 

	double m_divApy;

	int m_voltype; // 볼 type 0:상수, 1;atm term vol, 2: local 

	double m_volbval[AsNUM];	

	int m_volTN;
	int m_volSN;
	int *m_volBDay;
	double *m_volSmness;
	double ***m_vol;  //m_vol[AsNUM][m_volSN][m_volTN]; 

	int m_corrtype;
	int m_corrN;
	int *m_corrBDay;
	double *m_corr;	 
	 
	 
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double mc_price();
	 
	void curttimeVol(double *vol, int day, int Asidx);  // vol[AsNUM][m_volSN];

	int curttimeCorridx(int day);

	int findinteridx(double x, double *Dx, int DxN);
};


