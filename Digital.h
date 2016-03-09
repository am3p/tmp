#include "Basic.h"
 

#pragma once
using namespace std;

class Digital : public Basic
{
public:
	Digital(int cpFlag, double sval, double bval, double xval, double Cpn,
	      int matu, double ctime,
		  int evalN, double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int batchFlag, char *kCode,
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
		  double shiftRate, double spreadRate);

	virtual ~Digital();

	double getValue();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);

	double m_Ogreeks[MAXGREEKS1];	
	double m_Otvegas[MAXVEGAS];
	double m_Otvannas[MAXVEGAS];
	double m_Otzommas[MAXVEGAS];
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	 
	 
private:
	 
	char *m_kCode;
	int m_cpFlag;  // call=1 or put=-1 관련 Flag
	double m_sval; // 현재가
	double m_bval; // 기준가
	double m_xval; // 행사가
	double m_Cpn;

	int m_matu;    // 잔존만기
	double m_ctime; // 계산 호출시간 종가전 10, 종가후 16 정도?
	
	 
	int m_evalN;
	double *m_psval;
	
	int m_irateN;   // 금리 term 개수
	int *m_irateBDay; 
	int *m_irateCDay;
	double *m_iRFrate;
	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate; // 연속 배당
	double m_divbval;
	int m_divN;
	int *m_divBDay;
	double *m_divCash;
	double m_divApy;

	int m_voltype; // 볼 type 0:상수, 1;atm term vol, 2: local 
	double m_volbval;
	int m_volTN;
	int m_volSN;
	int *m_volBDay;
	double *m_volSmness;
	double **m_vol;
	
	 

	int m_batchFlag;
	double m_Sminlevel;
	double m_Smaxlevel;
	int m_meshSM1;
	int m_meshSM2;
	int m_meshSM3;
	int m_meshTN1;
	int m_meshTN2;

	double m_shiftRate;
	double m_spreadRate;
	double m_spreadS;
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks);
	void curttimeVol(double *vol, int day);

};


