#include "Basic.h"
 
 

#pragma once
using namespace std;

	 
class TwoAssetDigital : public Basic
{
public:
	TwoAssetDigital(int cpFlag, double *sval, double *bval, double *xval, double Cpn,
	      int matu, double ctime,
		  int evalN, double *psval, 
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int corrtype, int corrN, int *corrBDay, double *corr,
		  int batchFlag, char *kCode,
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
		  double shiftRate, double spreadRate);

	virtual ~TwoAssetDigital();

	double getValue();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);
	double getCorrDelta(int corrdeltaN,int *corrdeltaidx); 

	double m_Ogreeks[MAXGREEKS2];	
	double m_Otvegas1[MAXVEGAS];
	double m_Otvannas1[MAXVEGAS];
	double m_Otzommas1[MAXVEGAS];		
	double m_Otvegas2[MAXVEGAS];
	double m_Otvannas2[MAXVEGAS];
	double m_Otzommas2[MAXVEGAS];
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	double m_Otcorrdelta[MAXCORRDELTA];
	 
	 
private:
	 
	char *m_kCode;
	int m_cpFlag;  // call=1 or put=-1 관련 Flag
	double m_sval1; // 현재가
	double m_sval2;
	double m_bval1; // 기준가
	double m_bval2;
	double m_xval1; // 행사가
	double m_xval2;
	double m_Cpn;

	int m_matu;    // 잔존만기
	double m_ctime; // 계산 호출시간 종가전 10, 종가후 16 정도?
	
	 
	int m_evalN;
	double *m_psval1;
	double *m_psval2;
	
	int m_irateN;   // 금리 term 개수
	int *m_irateBDay; 
	int *m_irateCDay;
	double *m_iRFrate;
	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate1; // 연속 배당
	double m_divbval1;
	int m_divN1;
	int *m_divBDay1;
	double *m_divCash1;

	double m_divrate2; // 연속 배당
	double m_divbval2;
	int m_divN2;
	int *m_divBDay2;
	double *m_divCash2;

	double m_divApy;

	int m_voltype; // 볼 type 0:상수, 1;atm term vol, 2: local 
	double m_volbval1;
	double m_volbval2;
	int m_volTN;
	int m_volSN;
	int *m_volBDay;
	double *m_volSmness;
	double **m_vol1;
	double **m_vol2;

	int m_corrtype;
	int m_corrN;
	int *m_corrBDay;
	double *m_corr;	 

	int m_batchFlag;
	double m_Sminlevel1;
	double m_Smaxlevel1;
	double m_Sminlevel2;
	double m_Smaxlevel2;

	int m_meshSM1;
	int m_meshSM2;
	int m_meshSM3;
	int m_meshTN1;
	int m_meshTN2;

	double m_shiftRate;
	double m_spreadRate;
	double m_spreadS1;
	double m_spreadS2;
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks); 
	void curttimeVol1(double *vol, int day);
	void curttimeVol2(double *vol, int day);
	double curttimeCorr(int day);

};


