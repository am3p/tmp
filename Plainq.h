#include "Basic.h"
 

#pragma once
using namespace std;

class Plainq : public Basic
{
public:
	Plainq(int cpFlag,  double sval,  double xval, int matu,  		   
		  int irateN, int *irateCDay, double *iRFrate,  
		  int divrateN, int *divrateCDay, double *divrate, 
		  double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		   
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN);

	virtual ~Plainq();

	double getValue();
	double getcfValue();
	  
	double m_Ogreeks[MAXGREEKS1];	
	 
	 
	 
private:
	 
	 
	int m_cpFlag;  // call=1 or put=-1 관련 Flag
	 
	double m_sval; // 현재가
	 
	double m_xval; // 행사가
	int m_matu;    // 잔존만기
	  
	 
	int m_irateN;   // 금리 term 개수
	 
	int *m_irateCDay;
	double *m_iRFrate;
	 
	int m_divrateN;
	int *m_divrateCDay;
	double *m_divrate; // 연속 배당

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
	
	 
	 
	double m_Sminlevel;
	double m_Smaxlevel;
	int m_meshSM1;
	int m_meshSM2;
	int m_meshSM3;
	int m_meshTN1;
	int m_meshTN2;
	
	///////////////////////////
	 
	 				
	 
	double fd_price(double *greeks);
	 
	double volST(double S, int T);
	double volinter2d(double S, int T, int SN, int TN, double *SArr, int *TArr, double **VArr); 
	void curttimeVol(double *vol, int day);
};


