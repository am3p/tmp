#include "Basic.h"
 

#pragma once
using namespace std;

class FRN : public Basic
{
public:
	//int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	FRN(double ctime,	
		  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
		  int irateN, int *irateBDay,   double *iRFrate, double *iDCrate, double *iIRSrate);

	virtual ~FRN();

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
	 
 
	double m_sval; // ���簡
	double m_bval; // ���ذ�
	double *m_xval; // ��簡
	double m_kihval; //ki barrier
	double m_kiFlag;
	double *m_Cpn; // Cpn & dummy
	
	int m_matuN;
	int *m_matu;    // ��������
	double m_ctime; // ��� ȣ��ð� ������ 10, ������ 16 ����?

	int m_maxevalN;
	int *m_evalN;
	double *m_psval;

	//cd
	int m_CDN;
	double m_CDadjval;
	int *m_ObBDay;
	int *m_PayBDay;
	double *m_FlegFactor;
	double *m_ObCDrate;


	
	int m_irateN;   // �ݸ� term ����
	int *m_irateBDay; 
	 
	double *m_iRFrate;
	double *m_iDCrate;
	double  *m_iIRSrate;

	double m_divrate; // ���� ���
	double m_divbval;
	int m_divN;
	int *m_divBDay;
	double *m_divCash;
	double m_divApy;

	int m_voltype; // �� type 0:���, 1;atm term vol, 2: local 
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
	double *m_spreadXS;
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks);
	void curttimeVol(double *vol, int day);

};

