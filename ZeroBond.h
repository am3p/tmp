#include "Basic.h"
 

#pragma once
using namespace std;

// ���� binarybarrier�� �̿��Ͽ� in flat ���·� ����
class ZeroBond : public Basic
{
public:
	ZeroBond(int matu, 		  		  
	      int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate);

	virtual ~ZeroBond();

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
	int m_duFlag; // down = 1 or up = -1 ���� Flag
	int m_ioFlag; // in = 1 or out = -1 ���� Flag     ----->>> ���� 1, �ڰ� -1 ����!!
	int m_cpFlag;  // call=1 or flat(bond) = 0 or put=-1 ���� Flag      //��ü 12
	int m_ZeroBondCode; // m_ZeroBondCode = 6*max(m_duFlag,0) + 3*max(m_ioFlag,0) + 1*(m_cpFlag+1) ����
	  /*    11	13	    'DIDC'
			10	5	    'DIF'
			9	17	    'DIDP'
			8	21	    'DODC'
			7	9	    'DOF'
			6	25	    'DODP'
			5	14	    'UIDC'
			4	6	    'UIF'
			3	18	    'UIDP'
			2	22	    'UODC'
			1	10	    'UOF'
			0	26	    'UODP'   */


	                   
	int m_knockFlag; // KI or KO flag 
	double m_sval; // ���簡
	double m_bval; // ���ذ�
	double m_xval; // ��簡
	double m_hval; //barrier
	double m_Cpn; // Cash
	 

	int m_matu;    // ��������
	double m_ctime; // ��� ȣ��ð� ������ 10, ������ 16 ����?
	
	 
	int m_evalN;
	double *m_psval;
	
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
	double m_shiftXRate;
	double m_spreadHRate;
	double m_spreadXRate;
	double m_spreadHS;
	double m_spreadXS;
	
	///////////////////////////
	int m_curtgreek;
	int m_curtvegapertubation;
	int	m_curtvegaidx;
	 				
	
	double cf_price();
	double fd_price(double *greeks);
	void curttimeVol(double *vol, int day);

};


