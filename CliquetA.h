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
	
	//int m_modeltype; ����� payoff = NA * prate * max(sum(TRi) , GF) �Ѱ���
	 
	 int m_modeltype;

private:
	 
	double m_xalpha; // �������� ��� ��簡 ���� ������.. 1 ��簡 102% �� 1.02��...
 
	double m_sval; // ���簡
	double *m_bval; // ���ذ�

	double m_prate;
	double m_GFloor;
	double m_LCap;
	double m_LFloor;
	double m_NAQ;

	 
	int m_matuN;
	int *m_matu;    // ��������
	double m_ctime; // ��� ȣ��ð� ������ 10, ������ 16 ����?

	//int m_maxevalN;
	//int *m_evalN;
	//double *m_psval;
	
	int m_irateN;   // �ݸ� term ����
	int *m_irateBDay; 
	int *m_irateCDay;
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


