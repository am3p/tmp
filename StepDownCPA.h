#include "Basic.h"
 

#pragma once
using namespace std;

class StepDownCPA : public Basic
{
public:
	StepDownCPA(int kiFlag, int koFlag, double sval, double bval, double *xval, double kihval, double kohval, double *Cpn,
	      int *matuN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int batchFlag, char *kCode,
		  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
		  double shiftHRate, double spreadHRate, double *spreadXRate);

	virtual ~StepDownCPA();

	double getValue();
	double getVega(int vegaN, int *vegaidx);
	double getRho(int rhoN,int *rhoidx);

	double m_Ogreeks[MAXGREEKS1];	
	double m_Otvegas[MAXVEGAS];
	double m_Otvannas[MAXVEGAS];
	double m_Otzommas[MAXVEGAS];
	double m_Otrfrhos[MAXRHOS];
	double m_Otdcrhos[MAXRHOS];
	
	int m_modeltype;
	 
	 
private:
	 
 
	double *m_payxval;
	double *m_payCpn;
	int m_paymatuN;
	int *m_paymatu;


	double m_sval; // 현재가
	double m_bval; // 기준가
	double *m_xval; // 행사가
	double m_kihval; //ki barrier
	int m_kiFlag;
	double *m_Cpn; // Cpn & dummy

	double m_kohval;
	int m_koFlag;
	
	int m_matuN;
	int *m_matu;    // 잔존만기
	double m_ctime; // 계산 호출시간 종가전 10, 종가후 16 정도?

	int m_maxevalN;
	int *m_evalN;
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
	char *m_kCode;

	double m_Sminlevel;
	double m_Smaxlevel;
	int m_meshSM1;
	int m_meshSM2;
	int m_meshSM3;
	int m_meshTN1;
	int m_meshTN2;

	double m_startS;
	double m_endS;

	 
	double m_shiftHRate;
	double m_spreadHRate;
	double *m_spreadXRate;
	double m_spreadHS;
	//double *m_spreadXS;
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks);
	void curttimeVol(double *vol, int day);

};


