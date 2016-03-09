#include "Basic.h"
#include "odysseyRandom.h"
 
#define AsNUM3 3

#pragma once
using namespace std;
 
class StepDownCPCDC : public Basic
{
public:
	StepDownCPCDC(int kiFlag, int koFlag, double *sval, double *bval, double *xval, double *kihval, double *kohval, double *Cpn,
	      int *matuN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  double *CDlegInfo,
		  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int *corrtype, int *corrN, int *corrBDay, double *corr, long SimulNo);
/*
(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr);
*/
	virtual ~StepDownCPCDC();

	 
	double getValue();
	double getDeltaGammaTheta();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);
	//double getCorrDelta(int corrdeltaN,int *corrdeltaidx); 

	double m_Ogreeks[MAXGREEKS3];	
	double m_Otvegas[AsNUM3*MAXVEGAS];	 
	double m_Otrfrhos[(AsNUM3+1)*MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	 

	int m_modeltype;

	int m_iRateSeq[AsNUM3+1];    // <<-QT
 
	//int m_iRateSeqCD;
	 
private:	

	odysseyRandom *rnd;
	long m_noiters;

	//pay
	int m_paymatuN;
	double **m_payxval;
	double *m_payCpn;
	int *m_paymatu;

 
	 //CD
	double m_CDadjval;
	int m_CDN;
	int *m_ObBDay;
	int *m_PayBDay;
	double *m_ObCDrate;
	double *m_FlegFactor;


	 
	double m_sval[AsNUM3]; // ���簡
 
	double m_bval[AsNUM3]; // ���ذ�
	 
	double **m_xval;  //m_xval[AsNUM3][m_matuN]
	
	double m_kihval[AsNUM3]; //ki barrier
	 
	int m_kiFlag;
	double *m_Cpn; // Cpn & dummy

	 
	
	int m_matuN;
	int *m_matu;    // ��������
	 

	int m_maxevalN;
	int *m_evalN;

	double **m_psval;  //m_psval[AsNUM3][m_maxevalN];
	

	int m_irateN;   // �ݸ� term ����
	int *m_irateBDay; 
 
	double **m_iRFrate; // m_iRFrate[AsNUM3][m_irateN];

	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate[AsNUM3]; // ���� ���
	double m_divbval[AsNUM3];
	int m_divN[AsNUM3];
	int m_MAXdivN;
	int **m_divBDay;  //m_divBDay[AsNUM3][m_MAXdivN]
	double **m_divCash; // m_divCash[AsNUM3][m_MAXdivN] 

	double m_divApy;

	int m_voltype; // �� type 0:���, 1;atm term vol, 2: local 
	double m_volbval[AsNUM3];	

	int m_volTN;
	int m_volSN;
	int *m_volBDay;
	double *m_volSmness;
	double ***m_vol;  //m_vol[AsNUM3][m_volSN][m_volTN]; 
	
	
	int m_FXvoltype; // �� type 0:���, 1;atm term vol, 2: local 
	double m_FXvolbval[AsNUM3];
	
	int m_FXvolTN;
	int m_FXvolSN;
	int *m_FXvolBDay;
	double *m_FXvolSmness;
	double ***m_FXvol;  //m_FXvol[AsNUM3][m_FXvolSN][m_FXvolTN];	 

	int m_corrtype;
	int m_corrN;
	int *m_corrBDay;
	double ***m_corr;	 //m_corr[m_corrN][AsNUM3][AsNUM3];
	
	int m_FXcorrtype;
	int m_FXcorrN;
	int *m_FXcorrBDay;
	double **m_FXcorr;	  //m_FXcorr[AsNUM3][m_FXcorrN];
		 
	 
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	 
	double mc_price();

	void curttimeVol(double *vol, int day, int Asidx);  // vol[AsNUM3][m_volSN];
	 
	int curttimeCorridx(int day); // return [AsNUM3][AsNUM3];
	
	void curttimeFXVol(double *vol, int day, int Asidx); // vol[AsNUM3][m_FXvolSN];
 
	double curttimeFXCorr(int day, int Asidx); // return [AsNUM3]
	
	int findinteridx(double x, double *Dx, int DxN);
};


