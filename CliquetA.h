#include "Basic.h"
 

#pragma once
using namespace std;

class CliquetA : public Basic
{
public:
	CliquetA(double *sval, double *bval, double prate, double GFloor, double LCap, double LFloor, double NAQ,
	      int matuN, int *matu, double ctime,		  
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int batchFlag, char *kCode,
		  double Xmaxlevel, double Zmesh, int XmeshM, int *TmeshN);

	virtual ~CliquetA();

	double getValue();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);

	double m_Ogreeks[MAXGREEKS1];	
	double m_Otvegas[MAXVEGAS];
	double m_Otvannas[MAXVEGAS];
	double m_Otzommas[MAXVEGAS];
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	
	//int m_modeltype; 현재는 payoff = NA * prate * max(sum(TRi) , GF) 한가지
	 
	 int m_modeltype;

private:
	 
	double m_xalpha; // 기준지수 대비 행사가 비율 보통은.. 1 행사가 102% 면 1.02로...
 
	double m_sval; // 현재가
	double *m_bval; // 기준가

	double m_prate;
	double m_GFloor;
	double m_LCap;
	double m_LFloor;
	double m_NAQ;

	 
	int m_matuN;
	int *m_matu;    // 잔존만기
	double m_ctime; // 계산 호출시간 종가전 10, 종가후 16 정도?

	//int m_maxevalN;
	//int *m_evalN;
	//double *m_psval;
	
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
	char *m_kCode;

	double m_Xmaxlevel;
	 
	int m_ZFactor;
	double m_Zmesh;

	int m_meshZL;
	int m_meshXM;
	int m_meshTN1;
	int m_meshTN2;
 
	   
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks);
	void curttimeVol(double *vol, int day);
	double realS(double x, double z, double s0, double z0);

};


