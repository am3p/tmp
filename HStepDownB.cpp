#include "HStepDownB.h"
#include <iostream>
#include <fstream>

 

HStepDownB::HStepDownB(int kiFlag, int koFlag, int ki1Flag, int ki2Flag, double *sval, double *bval, double *xval, double *kihval, double *kohval, double *Cpn,
	      int matuN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int corrtype, int corrN, int *corrBDay, double *corr,
		  int batchFlag, char *kCode,
		  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
		  double shiftHRate, double spreadHRate, double *spreadXRate)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	int i, j;	 
	 
	//입력 데이터 내부 전역 변수로 저장   
	  
	m_matuN = matuN;
	m_matu = new int[m_matuN];
	m_sval1 = sval[0];
	m_sval2 = sval[1];
	m_bval1 = bval[0];	
	m_bval2 = bval[1];
	m_xval1 = new double[m_matuN];
	m_xval2 = new double[m_matuN];
	
	m_kihval1 = kihval[0];
	m_kihval2 = kihval[1];
	//m_kiFlag = kiFlag;

	m_kohval1 = kohval[0];  // koadd
	m_kohval2 = kohval[1];  // koadd
	//m_koFlag = koFlag;  // koadd

	m_ki1Flag = ki1Flag; //SH add
	m_ki2Flag = ki2Flag;

	m_Cpn = new double[m_matuN+2];

	m_shiftHRate = shiftHRate;
	m_spreadHRate = abs(spreadHRate);
	m_spreadXRate = new double[m_matuN];
	 
	//m_spreadXS1 = new double[m_matuN];
	//m_spreadXS2 = new double[m_matuN];

	for(i=0; i<m_matuN; i++) {
		m_xval1[i] = xval[i];
		m_xval2[i] = xval[m_matuN+i];
		m_Cpn[i] = Cpn[i];
		m_matu[i] = matu[i];
		m_spreadXRate[i] = abs(spreadXRate[i]);				
	}
	m_Cpn[m_matuN] = Cpn[m_matuN]; // dummy cpn
	m_Cpn[m_matuN+1] = Cpn[m_matuN+1]; // 원금보장 조정 cpn
		 
	m_ctime = ctime;

	m_maxevalN = maxevalN;
	m_evalN = new int[m_matuN];
	for(i=0; i<m_matuN; i++) 
		m_evalN[i] = evalN[i];
	m_psval1 = new double[m_maxevalN];
	m_psval2 = new double[m_maxevalN];
	for(i=0; i<m_maxevalN; i++) {
		m_psval1[i] = psval[i];	 
		m_psval2[i] = psval[m_maxevalN+i];
	}
	
	m_irateN = irateN;
	m_irateBDay = new int[m_irateN];
	m_irateCDay = new int[m_irateN];
	m_iRFrate = new double[m_irateN];
	m_iDCrate = new double[m_irateN];
	m_iIRSrate = new double[m_irateN];
	for(i=0; i<m_irateN; i++) {
		m_irateBDay[i] = irateBDay[i];
		m_irateCDay[i] = irateCDay[i];
		m_iRFrate[i] = iRFrate[i];
		m_iDCrate[i] = iDCrate[i];
		m_iIRSrate[i] = iIRSrate[i];
	}

	m_divrate1 = divrate[0];
	m_divrate2 = divrate[1];

	m_divbval1 = divbval[0];
	m_divbval2 = divbval[1];

	m_divN1 = divN[0];
	m_divN2 = divN[1];

	m_divBDay1 = new int[m_divN1];
	m_divBDay2 = new int[m_divN2];
	m_divCash1 = new double[m_divN1];
	m_divCash2 = new double[m_divN2];

	for(i=0; i<m_divN1; i++) {
		m_divBDay1[i] = divBDay[i];		
		m_divCash1[i] = divCash[i];
	}	
	for(i=0; i<m_divN2; i++) {
		m_divBDay2[i] = divBDay[m_divN1+i];		
		m_divCash2[i] = divCash[m_divN1+i];
	}

	m_divApy = divApy;

	m_voltype = voltype;
	m_volbval1 = volbval[0];
	m_volbval2 = volbval[1];
	m_volTN = volTN;
	m_volSN = volSN;
	m_volBDay = new int[m_volTN];
	for(i=0; i<m_volTN; i++)
		m_volBDay[i] = volBDay[i];
	m_volSmness = new double[m_volSN];
	for(i=0; i<m_volSN; i++)
		m_volSmness[i] = volSmness[i];

	m_vol1 = new double*[m_volSN]; // m_volSN은.. 1이상
	m_vol2 = new double*[m_volSN];
	for(i=0; i<m_volSN; i++) {
		m_vol1[i] = new double[m_volTN];
		m_vol2[i] = new double[m_volTN];
	}
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++) {
			m_vol1[i][j] = vol[i*m_volTN+j];
			m_vol2[i][j] = vol[m_volSN*m_volTN + i*m_volTN+j];
		}

	m_corrtype = corrtype;
	m_corrN = corrN;
	m_corrBDay = new int[m_corrN];
	m_corr = new double[m_corrN];
	for(i=0; i<m_corrN; i++) {
		m_corrBDay[i] = corrBDay[i];
		m_corr[i] = corr[i];
	}
    
	m_batchFlag = batchFlag;
	m_kCode = new char[kCodeN+1];
	strcpy_s(m_kCode,kCodeN,kCode);


	m_Sminlevel1 = Sminlevel;
	m_Sminlevel2 = Sminlevel;
	m_Smaxlevel1 = Smaxlevel;
	m_Smaxlevel2 = Smaxlevel;	
	
	m_meshSM1 = SmeshM[0];
	m_meshSM2 = SmeshM[1];
	m_meshSM3 = SmeshM[2];

	m_meshTN1 = TmeshN[0];
	m_meshTN2 = TmeshN[1];


	 
	m_kihval1 = m_kihval1 + m_bval1*m_shiftHRate;
	m_kihval2 = m_kihval2 + m_bval2*m_shiftHRate;

	m_spreadHS1 = m_kihval1 - m_bval1*m_spreadHRate;
	m_spreadHS2 = m_kihval2 - m_bval2*m_spreadHRate;

	/*
	for(i=0; i<m_matuN; i++) { 	  
		m_spreadXS1[i] = m_xval1[i] - m_bval1*m_spreadXRate[i];
		m_spreadXS2[i] = m_xval2[i] - m_bval2*m_spreadXRate[i];
	}
	*/
	
	m_startS1 = m_spreadHS1;
	m_startS2 = m_spreadHS2;
	//m_endS = m_xval[m_matuN-2]; // 마지막 조기상환일에 필요 조기상환일에 세팅!! 



	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/kzhere1.txt",ios::out);	 
	 
			 
			logsgnuplot << "m_ki1Flag = " << m_ki1Flag << endl;
			
			logsgnuplot << "m_ki2Flag = " << m_ki2Flag << endl;

		 
			logsgnuplot.close();

 
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	 

	// ki 체크
	if(m_sval1 <= m_kihval1)			
		m_ki1Flag = 1;

	if(m_sval2 <= m_kihval2)			
		m_ki2Flag = 1;
	 
	  
}

HStepDownB::~HStepDownB()
{	
	int i;	 
	delete m_evalN;
		 
	delete m_psval1;
	delete m_psval2;

	delete m_irateBDay;
	delete m_irateCDay;  
	delete m_iRFrate; 
	delete m_iDCrate;  
	delete m_iIRSrate; 	 
			
	delete m_divBDay1;
	delete m_divCash1;
	delete m_divBDay2;
	delete m_divCash2;
			
	delete m_volBDay;	 
	delete m_volSmness;	 

	for(i=0; i<m_volSN; i++) {
		delete m_vol1[i];
		delete m_vol2[i];
	}	
	delete m_vol1;
	delete m_vol2;
	 
	delete m_corrBDay;
	delete m_corr;
   

	delete m_matu;
	delete m_xval1;
	delete m_xval2;
	delete m_Cpn;
	delete m_spreadXRate;
	//delete m_spreadXS1;
	//delete m_spreadXS2;

	delete m_kCode;

}

double HStepDownB::getValue() {
 
	double greeks[MAXGREEKS2];
	double value;
	int i;

	 
	m_curtgreek = 1;
	value = fd_price(greeks);
	 

	for(i=0;i<BasicGRsN2; i++) {
		m_Ogreeks[i] = greeks[i];
	}
 

	return value;

}


double HStepDownB::getVega(int vegaN, int *vegaidx) {

	double greeks[MAXGREEKS2];
	double value;
	int i,j,k,udi;

	double **tempv;
	
	int pertubationX[3];
	double pertubationY[3];		
	double curtpertubation;
	double Price[3][2];
	double VegaPertubationR;
	
	m_curtgreek = 2;
		
	tempv = new double*[m_volSN]; // m_volSN은.. 1이상
	for(i=0; i<m_volSN; i++)
		tempv[i] = new double[m_volTN];
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol1[i][j]; // vol 저장
	

	m_Ogreeks[twoVegaidx1] = 0;
	m_Ogreeks[twoVannaCidx1] = 0;
	m_Ogreeks[twoZommaCidx1] = 0;

	//vol + 계산
	if(m_voltype == 0)
		vegaN = 1;
	for(i=0; i<vegaN; i++) {
		m_curtvegaidx = i; ////batchFile을 만들 때.. 화일이름에 idx 넣기 위해 

		for(udi=0; udi<=1; udi++) { //udi = 0은 vol+ , udi = 1은 vol-
			m_curtvegapertubation =udi;  //batchFile을 만들 때.. 0은 화일이름 + 넣기위해
	 
			VegaPertubationR = pow(-1.0,udi)*VegaPertubation;

			if(m_voltype == 0 || vegaN == 1) { // constant vol vegaN 1 은 평행이동 vega
				for(j=0; j<m_volTN; j++) {
					curtpertubation = VegaPertubationR;
					for(k=0; k<m_volSN; k++)
						m_vol1[k][j] = tempv[k][j] + curtpertubation;
				}				 
			}
			else { // term vol을 갖는 imp vol 또는 local vol //else 1
				if(i == 0) {
					pertubationX[0] = m_volBDay[vegaidx[i]];
					pertubationX[1] = m_volBDay[vegaidx[i+1]];

					pertubationY[0] = VegaPertubationR;
					pertubationY[1] = 0.0;

					for(j=0; j<m_volTN; j++) {
						curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,2,1,0); //  // 2:Data 2개, 1:비균등, 0:양쪽 flat
						for(k=0; k<m_volSN; k++) 
							m_vol1[k][j] = tempv[k][j] + curtpertubation;
					}
				} //if(i == 0)
				else { // else 2
					if(i == vegaN-1) {
						pertubationX[0] = m_volBDay[vegaidx[i-1]];
						pertubationX[1] = m_volBDay[vegaidx[i]];

						pertubationY[0] = 0.0;
						pertubationY[1] = VegaPertubationR;

						for(j=0; j<m_volTN; j++) {
							curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,2,1,0); //  // 2:Data 2개, 1:비균등, 0:양쪽 flat
							for(k=0; k<m_volSN; k++) 
								m_vol1[k][j] = tempv[k][j] + curtpertubation;
						}
					} // if(i == vegaN-1)
					else { //else 3					
						pertubationX[0] = m_volBDay[vegaidx[i-1]];
						pertubationX[1] = m_volBDay[vegaidx[i]];
						pertubationX[2] = m_volBDay[vegaidx[i+1]];

						pertubationY[0] = 0.0;
						pertubationY[1] = VegaPertubationR;
						pertubationY[2] = 0.0;

						for(j=0; j<m_volTN; j++) {
							curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,3,1,0); //  // 3:Data 3개, 1:비균등, 0:양쪽 flat
							for(k=0; k<m_volSN; k++) 
								m_vol1[k][j] = tempv[k][j] + curtpertubation;
						}
					} //else 3 if(i == vegaN-1)
				} //else 2 if(i == 0)
			}// else 1 if(m_voltype == 0)

			value = fd_price(greeks);
			for(j=0; j<3; j++) 			
				Price[j][udi] = greeks[j]; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+
		} //for(vi)
		//
		/*
		m_Otvegas1[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas1[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas1[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/
		
		//2012-05-23 1% vanna cash / 1% zomma cash
		m_Otvegas1[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas1[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/4; // term 1% vanna cash
		m_Otzommas1[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/2; // term 1% zomma cash

		m_Ogreeks[twoVegaidx1] += m_Otvegas1[i];
		m_Ogreeks[twoVannaCidx1] += m_Otvannas1[i];
		m_Ogreeks[twoZommaCidx1] += m_Otzommas1[i];
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol1[i][j] = tempv[i][j]; // vol 복원

	/////////////////////v1 end

		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol2[i][j]; // vol 저장
	

	m_Ogreeks[twoVegaidx2] = 0;
	m_Ogreeks[twoVannaCidx2] = 0;
	m_Ogreeks[twoZommaCidx2] = 0;

	//vol + 계산
	if(m_voltype == 0)
		vegaN = 1;
	for(i=0; i<vegaN; i++) {
		m_curtvegaidx = i; ////batchFile을 만들 때.. 화일이름에 idx 넣기 위해 

		for(udi=0; udi<=1; udi++) { //udi = 0은 vol+ , udi = 1은 vol-
			m_curtvegapertubation =udi + 2;  //batchFile을 만들 때.. 0은 화일이름 + 넣기위해
	 
			VegaPertubationR = pow(-1.0,udi)*VegaPertubation;

			if(m_voltype == 0 || vegaN == 1) { // constant vol vegaN 1 은 평행이동 vega
				for(j=0; j<m_volTN; j++) {
					curtpertubation = VegaPertubationR;
					for(k=0; k<m_volSN; k++)
						m_vol2[k][j] = tempv[k][j] + curtpertubation;
				}				 
			}
			else { // term vol을 갖는 imp vol 또는 local vol //else 1
				if(i == 0) {
					pertubationX[0] = m_volBDay[vegaidx[i]];
					pertubationX[1] = m_volBDay[vegaidx[i+1]];

					pertubationY[0] = VegaPertubationR;
					pertubationY[1] = 0.0;

					for(j=0; j<m_volTN; j++) {
						curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,2,1,0); //  // 2:Data 2개, 1:비균등, 0:양쪽 flat
						for(k=0; k<m_volSN; k++) 
							m_vol2[k][j] = tempv[k][j] + curtpertubation;
					}
				} //if(i == 0)
				else { // else 2
					if(i == vegaN-1) {
						pertubationX[0] = m_volBDay[vegaidx[i-1]];
						pertubationX[1] = m_volBDay[vegaidx[i]];

						pertubationY[0] = 0.0;
						pertubationY[1] = VegaPertubationR;

						for(j=0; j<m_volTN; j++) {
							curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,2,1,0); //  // 2:Data 2개, 1:비균등, 0:양쪽 flat
							for(k=0; k<m_volSN; k++) 
								m_vol2[k][j] = tempv[k][j] + curtpertubation;
						}
					} // if(i == vegaN-1)
					else { //else 3					
						pertubationX[0] = m_volBDay[vegaidx[i-1]];
						pertubationX[1] = m_volBDay[vegaidx[i]];
						pertubationX[2] = m_volBDay[vegaidx[i+1]];

						pertubationY[0] = 0.0;
						pertubationY[1] = VegaPertubationR;
						pertubationY[2] = 0.0;

						for(j=0; j<m_volTN; j++) {
							curtpertubation = interp1(m_volBDay[j],pertubationX,pertubationY,3,1,0); //  // 3:Data 3개, 1:비균등, 0:양쪽 flat
							for(k=0; k<m_volSN; k++) 
								m_vol2[k][j] = tempv[k][j] + curtpertubation;
						}
					} //else 3 if(i == vegaN-1)
				} //else 2 if(i == 0)
			}// else 1 if(m_voltype == 0)

			value = fd_price(greeks);
			for(j=0; j<3; j++) 			
				Price[j][udi] = greeks[j+3]; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+
		} //for(vi)
		//
		/*
		m_Otvegas2[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas2[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas2[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/

		//2012-05-23 1% vanna cash / 1% zomma cash
		m_Otvegas2[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas2[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/4; // term 1% vanna cash
		m_Otzommas2[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/2; // term 1% zomma cash

		m_Ogreeks[twoVegaidx2] += m_Otvegas2[i];
		m_Ogreeks[twoVannaCidx2] += m_Otvannas2[i];
		m_Ogreeks[twoZommaCidx2] += m_Otzommas2[i];
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol2[i][j] = tempv[i][j]; // vol 복원
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	 
	return m_Ogreeks[twoVegaidx1];

}


double HStepDownB::getRho(int rhoN,int *rhoidx) {
	
	double greeks[MAXGREEKS2];
	double value;

	double *tempirate;
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		 

	m_curtgreek = 3;
 


	tempirate = new double[m_irateN];
	//rf rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iRFrate[i]; // RF rate 저장
	}
	 
	m_Ogreeks[twoRhorfidx] = 0; //rf rho sum
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // 평행이동 rho
			for(j=0; j<m_irateN; j++) 
				m_iRFrate[j] = tempirate[j] + RhoPertubation;		
		}
		else {
			if(i == 0) {		
				pertubationX[0] = m_irateBDay[rhoidx[i]];
				pertubationX[1] = m_irateBDay[rhoidx[i+1]];

				pertubationY[0] = RhoPertubation;
				pertubationY[1] = 0.0;

				for(j=0; j<m_irateN; j++) {			
					m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
				}
			}
			else {
				if(i == rhoN-1) {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];			

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;			
			
					for(j=0; j<m_irateN; j++) {							
						m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat				
					}		
				}
				else {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];
					pertubationX[2] = m_irateBDay[rhoidx[i+1]];

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;	
					pertubationY[2] = 0.0;
			
					for(j=0; j<m_irateN; j++) {							
						m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat				
					}	

				}
			}
		}// if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otrfrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[twoRhorfidx] += m_Otrfrhos[i];
	}
 

	for(i=0; i<m_irateN; i++) {
		m_iRFrate[i] = tempirate[i]; // RF rate 복원
	}
	

//////////////////////////////
	//dc rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iDCrate[i]; // DC rate 저장
	}
	 
	m_Ogreeks[twoRhodcidx] = 0; // dc rho sum
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // 평행이동 rho
			for(j=0; j<m_irateN; j++) 
				m_iDCrate[j] = tempirate[j] + RhoPertubation;
		}
		else {
			if(i == 0) {		
				pertubationX[0] = m_irateBDay[rhoidx[i]];
				pertubationX[1] = m_irateBDay[rhoidx[i+1]];

				pertubationY[0] = RhoPertubation;
				pertubationY[1] = 0.0;

				for(j=0; j<m_irateN; j++) {			
					m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
				}
			}
			else {
				if(i == rhoN-1) {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];			

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;			
			
					for(j=0; j<m_irateN; j++) {							
						m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat				
					}		
				}
				else {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];
					pertubationX[2] = m_irateBDay[rhoidx[i+1]];

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;	
					pertubationY[2] = 0.0;
			
					for(j=0; j<m_irateN; j++) {							
						m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat				
					}	

				}
			}
		} // if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otdcrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[twoRhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate 복원
	}
 	


	delete tempirate;

		
	return m_Ogreeks[twoRhorfidx] + m_Ogreeks[twoRhodcidx];

}



double HStepDownB::getCorrDelta(int corrdeltaN,int *corrdeltaidx) {
	
	double greeks[MAXGREEKS2];
	double value;

	double *tempcorr;
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		 

	m_curtgreek = 4;
 


	tempcorr = new double[m_corrN];
	 
	for(i=0; i<m_corrN; i++) {
		tempcorr[i] = m_corr[i];  
	}
	 
	m_Ogreeks[twoCorrDeltaidx] = 0;  
	for(i=0; i<corrdeltaN; i++) {
		if(corrdeltaN == 1) { // 평행이동 rho
			for(j=0; j<m_corrN; j++) 
				m_corr[j] = tempcorr[j] + CorrDeltaPertubation;		
		}
		else {
			if(i == 0) {		
				pertubationX[0] = m_corrBDay[corrdeltaidx[i]];
				pertubationX[1] = m_corrBDay[corrdeltaidx[i+1]];

				pertubationY[0] = CorrDeltaPertubation;
				pertubationY[1] = 0.0;

				for(j=0; j<m_corrN; j++) {			
					m_corr[j] = tempcorr[j] + interp1(m_corrBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
				}
			}
			else {
				if(i == corrdeltaN-1) {					
					pertubationX[0] = m_corrBDay[corrdeltaidx[i-1]];			
					pertubationX[1] = m_corrBDay[corrdeltaidx[i]];			

					pertubationY[0] = 0.0;			
					pertubationY[1] = CorrDeltaPertubation;			
			
					for(j=0; j<m_corrN; j++) {							
						m_corr[j] = tempcorr[j] + interp1(m_corrBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat				
					}		
				}
				else {					
					pertubationX[0] = m_corrBDay[corrdeltaidx[i-1]];			
					pertubationX[1] = m_corrBDay[corrdeltaidx[i]];
					pertubationX[2] = m_corrBDay[corrdeltaidx[i+1]];

					pertubationY[0] = 0.0;			
					pertubationY[1] = CorrDeltaPertubation;	
					pertubationY[2] = 0.0;
			
					for(j=0; j<m_corrN; j++) {							
						m_corr[j] = tempcorr[j] + interp1(m_corrBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat				
					}	

				}
			}
		}// if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otcorrdelta[i] = value - m_Ogreeks[0]; // 0.1 term correlationdelta 
		m_Ogreeks[twoCorrDeltaidx] += m_Otcorrdelta[i];
	}
 

	for(i=0; i<m_corrN; i++) {
		m_corr[i] = tempcorr[i]; // RF rate 복원
	}
	
 

	delete tempcorr;

		
	return m_Ogreeks[twoCorrDeltaidx];

}

double HStepDownB::cf_price()
{	
	return 7777;
}

double HStepDownB::fd_price(double *greeks)
{
	 
		double DCValue12;
		double DCValue1;
		double DCValue2;
		
		int i, j, matuLoop, intraLoop;
		int SM1, SM2, SmeshType1, SmeshType2, intraTN;

		double dt;
		double dx, dx3, dy;

		double slevel1, slevel2, xlevel1, xlevel2, hlevel1, hlevel2;
		double spreadHSlevel1, spreadHSlevel2;//, spreadXSlevel1, spreadXSlevel2;
		double *S1, *h1, *S2, *h2;
		double **oldV, **newV, **tU;
		double **knockoldV, **knocknewV;	

		//SH
		double **knock1oldV, **knock1newV;
		double **knock2oldV, **knock2newV;

		double **XMAu, **XMAd, **XMAl;
		double **YMAu, **YMAd, **YMAl;
		double **knockXMAu, **knockXMAd, **knockXMAl;
		double **knockYMAu, **knockYMAd, **knockYMAl;
			
		double *vecV1, *vecV2;
		double **coeV1, **coeV2;
	 
		double startSlevel1, startSlevel2, endSlevel1, endSlevel2;

		int xidx1, xidx2;
		int startSidx1, startSidx2, endSidx1, endSidx2;
		int spreadHSidx1, spreadHSidx2;//, spreadXSidx1, spreadXSidx2;
		int hidx1, hidx2;
		int evali;
		int evalq, expiryq;


		double value;

		double *SrTemp1, *SrTemp2, **SrVTemp, **SrV1Temp; // greek 계산을 위한 변수 당일값, 다음날 값 저장
		double sleveltemp1, sleveltemp2;		
		int outN;	
		int Baseidx;
		int curtN;
	 
		double disc_div1, disc_div2;
		double Srv1, Srv2;
		double *volskew1, *volskew2;

		int IsMidCompu, mid_n;
 
		curtN = m_matuN;

		evalq = m_evalN[curtN-1] - m_matu[curtN-1];
	
		if(evalq <=1 ) {	
			slevel1 = m_sval1/m_bval1;
			slevel2 = m_sval2/m_bval2;
		}
		else {
			slevel1 = m_sval1/m_bval1 * (m_evalN[curtN-1] - (evalq-1));
			for(evali=1; evali<evalq; evali++)
				slevel1 = slevel1 + m_psval1[evali]/m_bval1;
			slevel1 = slevel1/m_evalN[curtN-1];

			slevel2 = m_sval2/m_bval2 * (m_evalN[curtN-1] - (evalq-1));
			for(evali=1; evali<evalq; evali++)
				slevel2 = slevel2 + m_psval2[evali]/m_bval2;
			slevel2 = slevel2/m_evalN[curtN-1];
		}	 
	 


		xlevel1 = m_xval1[curtN-1]/m_bval1;
		hlevel1 = m_kihval1/m_bval1;
		startSlevel1 = m_startS1/m_bval1;	 
		spreadHSlevel1 = m_spreadHS1/m_bval1;
		//spreadXSlevel1 = m_spreadXS1[curtN-1]/m_bval1;
	
		xlevel2 = m_xval2[curtN-1]/m_bval2;
		hlevel2 = m_kihval2/m_bval2;
		startSlevel2 = m_startS2/m_bval2;	
		spreadHSlevel2 = m_spreadHS2/m_bval2;
		//spreadXSlevel2 = m_spreadXS2[curtN-1]/m_bval2;
	 
		double meshxlevel1, meshxlevel2;
		for(i=0; i<m_matuN; i++) {
			if(m_matu[i]>=0) {
				meshxlevel1 = m_xval1[i]/m_bval1;
				meshxlevel2 = m_xval2[i]/m_bval2;
				break;
			}
		}
			
		dx = (m_Smaxlevel1 - m_Sminlevel1)/m_meshSM1;
		xidx1 = (int)floor((meshxlevel1-m_Sminlevel1)/dx);
	 
		
		// Adaptive mesh 조건 만족 못할 때는 균등 분할 적용
		if((m_meshSM2 <=0) || (m_meshSM3 <=0) || dx <= AdmeshR ) { // || (xidx1 - m_meshSM2 < 0) || (m_meshSM1 < xlevel1 + m_meshSM2)) { 
			SmeshType1 = 0; // 균등
			SM1 = m_meshSM1;
		 
			S1 = (double *)malloc((size_t)(SM1+1)*sizeof(double));
			//S mesh generate
			for(i=0; i<=SM1; i++) 
				S1[i] = m_Sminlevel1 + dx*i;
			h1 = (double *)malloc((size_t)(1)*sizeof(double)); // free(h)를 위한 더미 설정
			h1[0] = dx;
		
		}
		else {		
			SmeshType1 = 1; // Adaptive mesh
			if(m_meshSM3 == 1) {
				double *tempgenS;
				tempgenS = new double[m_meshSM1 + int(4*Admesh/AdmeshR)+4];

				int cri;
				double temph, temps;
			
				cri = 0;
				tempgenS[cri] = m_Sminlevel1;
				temph = dx;
				temps = tempgenS[cri] + temph;
				while (temps <= m_Smaxlevel1 + 2*myeps) {
					if(temps <= hlevel1 - Admesh + myeps)
						temph = dx;
					if(temps > hlevel1 - Admesh + myeps )
						temph = AdmeshR;
					if(temps > hlevel1 + Admesh + myeps)
						temph = dx;
					if(temps > meshxlevel1 - Admesh + myeps)
						temph = AdmeshR;
					if(temps > meshxlevel1 + Admesh + myeps)
						temph = dx;

					temps = tempgenS[cri] + temph;
					cri = cri + 1;
					tempgenS[cri] = temps;
					temps = tempgenS[cri] + temph;		
					if(tempgenS[cri] > m_Smaxlevel1)
						tempgenS[cri] = m_Smaxlevel1;
				}	

				SM1 = cri;
				S1 = (double *)malloc((size_t)(SM1+1)*sizeof(double));
				for(i=0; i<=SM1; i++)
					S1[i] = tempgenS[i];

				delete tempgenS;			
			}
			else {
				SM1 = m_meshSM1 + 2*m_meshSM2*(m_meshSM3-1);
				dx3 = dx/m_meshSM3;

				S1 = (double *)malloc((size_t)(SM1+1)*sizeof(double));
				//S mesh generate
				for(i=0; i<=m_meshSM1; i++) {
					if(i<=xidx1-m_meshSM2) {
						S1[i] = m_Sminlevel1 + dx*i;
					}
					else {
						if(xidx1-m_meshSM2 < i && i<= xidx1 + m_meshSM2) {
							for(j=1; j<=m_meshSM3; j++) {					
								S1[xidx1-m_meshSM2 + m_meshSM3*(i-(xidx1-m_meshSM2+1)) +j] = S1[xidx1-m_meshSM2] + dx3*(m_meshSM3*(i-(xidx1-m_meshSM2+1))+j);
							}
						}
						else {
							S1[i+2*m_meshSM2*(m_meshSM3-1)] = S1[xidx1+m_meshSM2*(2*m_meshSM3-1)] + dx*(i-(xidx1+m_meshSM2));
						}
					}
				}
			
			}

			h1 = (double *)malloc((size_t)(SM1)*sizeof(double));
			for(i=0; i<SM1; i++) 
				h1[i] = S1[i+1] - S1[i];

		}

		dy = (m_Smaxlevel2 - m_Sminlevel2)/m_meshSM1;
		xidx2 = (int)floor((meshxlevel2-m_Sminlevel2)/dy);
 
		/*
		SmeshType2 = SmeshType1;
		SM2 = SM1;
		S2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
		memcpy(S2,S1,(size_t)((SM2+1)*sizeof(double)));
		if(SmeshType2 == 0) {
			h2 = (double *)malloc((size_t)(1)*sizeof(double)); // free(h)를 위한 더미 설정
			h2[0] = dy;
		}
		else {
			h2 = (double *)malloc((size_t)(SM2)*sizeof(double));
			memcpy(h2,h1,(size_t)((SM2)*sizeof(double)));
		}
 
         */  
			
		// Adaptive mesh 조건 만족 못할 때는 균등 분할 적용
		if((m_meshSM2 <=0) || (m_meshSM3 <=0) || dy <= AdmeshR ) { // || (xidx1 - m_meshSM2 < 0) || (m_meshSM1 < xlevel1 + m_meshSM2)) { 
			SmeshType2 = 0; // 균등
			SM2 = m_meshSM1;
		 
			S2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			//S mesh generate
			for(i=0; i<=SM2; i++) 
				S2[i] = m_Sminlevel2 + dy*i;
			h2 = (double *)malloc((size_t)(1)*sizeof(double)); // free(h)를 위한 더미 설정
			h2[0] = dy;
		
		}
		else {		
			SmeshType2 = 1; // Adaptive mesh
			if(m_meshSM3 == 1) {
				double *tempgenS;
				tempgenS = new double[m_meshSM1 + int(4*Admesh/AdmeshR)+4]; 

				int cri;
				double temph, temps;
			
				cri = 0;
				tempgenS[cri] = m_Sminlevel2;
				temph = dy;
				temps = tempgenS[cri] + temph;
				while (temps <= m_Smaxlevel2 + 2*myeps) {
					if(temps <= hlevel2 - Admesh + myeps)
						temph = dy;
					if(temps > hlevel2 - Admesh + myeps )
						temph = AdmeshR;
					if(temps > hlevel2 + Admesh + myeps)
						temph = dy;
					if(temps > meshxlevel2 - Admesh + myeps)
						temph = AdmeshR;
					if(temps > meshxlevel2 + Admesh + myeps)
						temph = dy;
		 


					temps = tempgenS[cri] + temph;
					cri = cri + 1;
					tempgenS[cri] = temps;
					temps = tempgenS[cri] + temph;		
					if(tempgenS[cri] > m_Smaxlevel2)
						tempgenS[cri] = m_Smaxlevel2;
				}	

				SM2 = cri;
				S2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
				for(i=0; i<=SM2; i++)
					S2[i] = tempgenS[i];

				delete tempgenS;			
			}
			else {
				SM2 = m_meshSM1 + 2*m_meshSM2*(m_meshSM3-1);
				dx3 = dy/m_meshSM3;

				S2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
				//S mesh generate
				for(i=0; i<=m_meshSM1; i++) {
					if(i<=xidx2-m_meshSM2) {
						S2[i] = m_Sminlevel2 + dy*i;
					}
					else {
						if(xidx2-m_meshSM2 < i && i<= xidx2 + m_meshSM2) {
							for(j=1; j<=m_meshSM3; j++) {					
								S2[xidx2-m_meshSM2 + m_meshSM3*(i-(xidx2-m_meshSM2+1)) +j] = S2[xidx2-m_meshSM2] + dx3*(m_meshSM3*(i-(xidx2-m_meshSM2+1))+j);
							}
						}
						else {
							S2[i+2*m_meshSM2*(m_meshSM3-1)] = S2[xidx2+m_meshSM2*(2*m_meshSM3-1)] + dy*(i-(xidx2+m_meshSM2));
						}
					}
				}
			
			}

			h2 = (double *)malloc((size_t)(SM2)*sizeof(double));
			for(i=0; i<SM2; i++) 
				h2[i] = S2[i+1] - S2[i];

		}
 
	 

		startSidx1 = findinteridx(startSlevel1, S1, SM1);
		spreadHSidx1 = findinteridx(spreadHSlevel1, S1, SM1);
		//spreadXSidx1 = findinteridx(spreadXSlevel1, S1, SM1);
		hidx1 = findinteridx(hlevel1, S1, SM1);
		xidx1 = findinteridx(xlevel1, S1, SM1);
	
		startSidx2 = findinteridx(startSlevel2, S2, SM2);
		spreadHSidx2 = findinteridx(spreadHSlevel2, S2, SM2);
		//spreadXSidx2 = findinteridx(spreadXSlevel2, S2, SM2);
		hidx2 = findinteridx(hlevel2, S2, SM2);
		xidx2 = findinteridx(xlevel2, S2, SM2);
	 	

		// greek 계산을 위한 작업
		if(m_curtgreek == 1 || m_curtgreek == 2) { // 기본 그릭		 		 

			outN = 2*SoutRange + 1; // +- 20%출력
			SrTemp1 = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // up down gamma를 위해 위 아래 2개씩 더 
			SrTemp2 = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // up down gamma를 위해 위 아래 2개씩 더 
			SrVTemp = (double **)malloc((size_t)(outN + 4)*sizeof(double *)); // 당일 저장
			SrV1Temp = (double **)malloc((size_t)(outN + 4)*sizeof(double *)); // 다음날 저장
			for(i=0; i<outN+4; i++) {
				SrVTemp[i] = (double *)malloc((size_t)(outN + 4)*sizeof(double));
				SrV1Temp[i] = (double *)malloc((size_t)(outN + 4)*sizeof(double));
			}

			for(i=0; i<outN + 4; i++) {
				SrTemp1[i] = m_sval1*(1-(SoutRange+2)*0.01 + 0.01*i); //m_sval의 78%~122%까지 1%씩 저장 0~44 45개 중심 22	
				SrTemp2[i] = m_sval2*(1-(SoutRange+2)*0.01 + 0.01*i); //m_sval의 78%~122%까지 1%씩 저장 0~44 45개 중심 22	
			}		 
		}


		
		oldV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		newV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		tU = (double **)malloc((size_t)(SM1+1)*sizeof(double *));

		XMAl = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		XMAd = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		XMAu = (double **)malloc((size_t)(SM1+1)*sizeof(double *));

		YMAl = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		YMAd = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		YMAu = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
	
		
		knockoldV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knocknewV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));	 

		//SH
		knock1oldV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knock1newV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));	
		knock2oldV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knock2newV = (double **)malloc((size_t)(SM1+1)*sizeof(double *));	

		knockXMAl = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knockXMAd = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knockXMAu = (double **)malloc((size_t)(SM1+1)*sizeof(double *));

		knockYMAl = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knockYMAd = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		knockYMAu = (double **)malloc((size_t)(SM1+1)*sizeof(double *));

		coeV1 = (double **)malloc((size_t)(SM1+1)*sizeof(double *));
		coeV2 = (double **)malloc((size_t)(SM1+1)*sizeof(double *));


		for(i=0; i<SM1+1; i++) {
			oldV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			newV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			tU[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));

			XMAl[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			XMAd[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			XMAu[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));

			YMAl[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			YMAd[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			YMAu[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
				
			knockoldV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			knocknewV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	

			//SH
			knock1oldV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			knock1newV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			knock2oldV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
			knock2newV[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));	
	 
			knockXMAl[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			knockXMAd[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			knockXMAu[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));

			knockYMAl[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			knockYMAd[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			knockYMAu[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));

			coeV1[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
			coeV2[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
		}
	 
		vecV1 = (double *)malloc((size_t)(SM1+1)*sizeof(double));	 
		vecV2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
 
		volskew1 = (double*)malloc((size_t)(m_volSN)*sizeof(double));
		volskew2 = (double*)malloc((size_t)(m_volSN)*sizeof(double));
	 
	
 
		// 만기 payoff setting
	
		expiryq = 1;

	 
		if(m_modeltype == 10) {
			
	 
			m_ki1Flag = 1;
			m_ki2Flag = 1;
	 
	 

						//KI12
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {						
					if(i>=xidx1 && j>=xidx2) 
						knocknewV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
					else
						if(i>hidx1 && j>hidx2)	// no KI				
							knocknewV[i][j] = (1+m_Cpn[curtN]); //dummy
						else {			
							if(i<=hidx1 && j<=hidx2) // KI12						
								knocknewV[i][j] = AMOUNT * (1-1.0*(1-MIN(S1[i],S2[j])));
							else //KI1 or KI2
								knocknewV[i][j] = AMOUNT * (1-0.5*(1-MIN(S1[i],S2[j])));
						}				
				}
			}

		 
		 
 
		}
		else {

			//KI12
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(i>=xidx1 && j>=xidx2) 
						knocknewV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
					else
						knocknewV[i][j] = AMOUNT * (1-1.0*(1-MIN(S1[i],S2[j])));				
				}
			}

		
			//KI1
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(i>=xidx1 && j>=xidx2) 
						knock1newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
					else
						if(j>hidx2)	// KI1만				
							knock1newV[i][j] = AMOUNT * (1-0.5*(1-MIN(S1[i],S2[j])));				
						else
							knock1newV[i][j] = AMOUNT * (1-1.0*(1-MIN(S1[i],S2[j])));
				}
			}
		
			//KI2
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(i>=xidx1 && j>=xidx2) 
						knock2newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
					else
						if(i>hidx1)	// KI2만				
							knock2newV[i][j] = AMOUNT * (1-0.5*(1-MIN(S1[i],S2[j])));				
						else
							knock2newV[i][j] = AMOUNT * (1-1.0*(1-MIN(S1[i],S2[j])));
				}
			}

			//KI0
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(i>=xidx1 && j>=xidx2) 
						newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
					else
						if(i>hidx1 && j>hidx2)	// no KI				
							newV[i][j] = (1+m_Cpn[curtN]); //dummy
						else {			
							if(i<=hidx1 && j<=hidx2) // KI12						
								newV[i][j] = AMOUNT * (1-1.0*(1-MIN(S1[i],S2[j])));
							else //KI1 or KI2
								newV[i][j] = AMOUNT * (1-0.5*(1-MIN(S1[i],S2[j])));
						}
				}
			}
		}


		/*
		if(m_spreadHRate > 0) { // KI 쪽 spread 만기에만 
			double spreadvalue;
			double spmin, spmax;
			double jumpvalue;
		
			//KI1
			spmin = S2[spreadHSidx2];				
			spmax = S2[hidx2];			
			if(spmin == spmax)
				spmin = S2[hidx2-1];

			for(i=0; i<=SM1; i++) {
				jumpvalue = AMOUNT * (1-0.5*(1-MIN(S1[i],S2[hidx2])));
				for(j=spreadHSidx2; j<=hidx2; j++) {				
					spreadvalue = jumpvalue/(spmax-spmin) * (S2[j]-spmin);						
					knock1newV[i][j] = MAX(knock1newV[i][j],spreadvalue);
				 
				}
			}
			
			
			//KI2
			spmin = S1[spreadHSidx1];				
			spmax = S1[hidx1];			
			if(spmin == spmax)
				spmin = S1[hidx1-1];
			for(j=0; j<=SM2; j++) {
				jumpvalue = AMOUNT * (1-0.5*(1-MIN(S1[hidx1],S2[j])));
				for(i=spreadHSidx1; i<=hidx1; i++) {				
					spreadvalue = jumpvalue/(spmax-spmin) * (S1[i]-spmin);						
					knock2newV[i][j] = MAX(knock2newV[i][j],spreadvalue);
				 
				}
			}
 
			 
			//KI0
			for(i=spreadHSidx1; i<=SM1; i++) {
				for(j=spreadHSidx2; j<=SM2; j++) {
					if(i<=hidx1 || j<=hidx2) {
						if(S1[i]<=S2[j]) {
							spmin = S1[spreadHSidx1];
							spmax = S1[hidx1];									
							if(spmin == spmax)				
								spmin = S1[hidx1-1];
						}
						else {
							spmin = S2[spreadHSidx2];
							spmax = S2[hidx2];													
							if(spmin == spmax)				
								spmin = S2[hidx2-1];
						}

						spreadvalue = AMOUNT * (1+m_Cpn[curtN])/(spmax-spmin) * (MIN(S1[i],S2[j])-spmin);
						newV[i][j] = MAX(newV[i][j],spreadvalue);
					}
				}
			}
		} // H spread
		*/
	  
		//if(m_spreadXRate[curtN-1] > 0 && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
		if(m_spreadXRate[curtN-1] > 0) { // X spread는 만기와 조기상환일에 가능
			double spreadvalue;
			double spmax;
		


			//KI0
			//for(i=spreadXSidx1; i<=SM1; i++) {
				//for(j=spreadXSidx2; j<=SM2; j++) {
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(i<=xidx1 || j<=xidx2) {
					 
						if(S1[i]<=S2[j]) {
							 
							spmax = S1[xidx1];													
							 
						}
						else {
						 
							spmax = S2[xidx2];													
						 
						}

						//spreadvalue = AMOUNT * (1+m_Cpn[curtN-1])/(spmax-spmin) * (MIN(S1[i],S2[j])-spmin);
						//spreadvalue = AMOUNT * ( (1+m_Cpn[curtN-1]-1)/(spmax-spmin) * (MIN(S1[i],S2[j])-spmin) + 1);
						spreadvalue = AMOUNT * ( m_spreadXRate[curtN-1] * (MIN(S1[i],S2[j])-spmax) + 1+m_Cpn[curtN-1]);
						knocknewV[i][j] = MAX(knocknewV[i][j], spreadvalue);
						knock1newV[i][j] = MAX(knock1newV[i][j], spreadvalue); //SH
						knock2newV[i][j] = MAX(knock2newV[i][j], spreadvalue);
						newV[i][j] = MAX(newV[i][j],spreadvalue);
					}
				}
			}
		} // X spread

	 
		// KI 종류에 따라..
		if(m_ki1Flag == 1 && m_ki2Flag == 1) { 
			for(i=0; i<=SM1; i++)
				memcpy(newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
			for(i=0; i<=SM1; i++)
				memcpy(knock1newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
			for(i=0; i<=SM1; i++)
				memcpy(knock2newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
		}
		
		if(m_ki1Flag == 1 && m_ki2Flag == 0) { 
			for(i=0; i<=SM1; i++)
				memcpy(newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
			for(i=0; i<=SM1; i++)
				memcpy(knock2newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
		}
		
		if(m_ki1Flag == 0 && m_ki2Flag == 1) { 
			for(i=0; i<=SM1; i++)
				memcpy(newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
			for(i=0; i<=SM1; i++)
				memcpy(knock1newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
		}
	
	 	
	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/kzhere2.txt",ios::out);	 
	 
			 
			logsgnuplot << "m_ki1Flag = " << m_ki1Flag << endl;
			
			logsgnuplot << "m_ki2Flag = " << m_ki2Flag << endl;

			logsgnuplot << "newV[49][100] = " << newV[49][100] << endl;
			logsgnuplot << "knocknewV[49][100] = " << knocknewV[49][100] << endl;


		 
			logsgnuplot.close();

 
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	  
  
	
		int sim_knock1Flag, sim_knock2Flag;
		if(m_curtgreek == 1 && (m_matu[curtN-1] == 1 || m_matu[curtN-1] == 0)) { // 기본 그릭	Theta를 위한 +1일 저장
			for(i=0; i<outN + 4; i++) {
				for(j=0; j<outN + 4; j++) {
					sim_knock1Flag = 0;
					sim_knock2Flag = 0;
					if(SrTemp1[i] <= m_kihval1)
						sim_knock1Flag = 1;
					if(SrTemp2[j] <= m_kihval2)
						sim_knock2Flag = 1;

					if(m_matu[curtN-1] == 1) {
						if(m_ctime>15) { // 장 종료 후 는.. 오늘 종가가 내일 psval의 [0]에 들어 가는 경우
							sleveltemp1 = SrTemp1[i]/m_bval1; // *(m_evalN[curtN-1] - evalq) = 1 = m_matu[curtN-1] 해당
							for(evali=0; evali<evalq; evali++)
								sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
							sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
												
							sleveltemp2 = SrTemp2[j]/m_bval2; // *(m_evalN[curtN-1] - evalq) = 1 = m_matu[curtN-1] 해당
							for(evali=0; evali<evalq; evali++)
								sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
							sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];
						}
						else { // 장중
							sleveltemp1 = (1+MIN(evalq,1))*SrTemp1[i]/m_bval1;					
							for(evali=1; evali<evalq; evali++)
								sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
							sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
												
							sleveltemp2 = (1+MIN(evalq,1))*SrTemp2[j]/m_bval2;					
							for(evali=1; evali<evalq; evali++)
								sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
							sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];
						}
					}
					else {
						sleveltemp1 = SrTemp1[i]/m_bval1;
						for(evali=1; evali<evalq; evali++)
							sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
						sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
					
						sleveltemp2 = SrTemp2[j]/m_bval2;
						for(evali=1; evali<evalq; evali++)
							sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
						sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];
					}

					if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
						SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knocknewV, SM1, SM2);  
					if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
						SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock1newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
						SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock2newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
						SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2);
				}		 
			}
		}
 	 

		double r, q1, q2, v1, v2, rd, cor;
		double *dailyIFrf, *dailyIFdc;	
		// intra Day에서는 입력 파라미터 고정
		// 하루마다 해당 파라미터 계산
		q1 = m_divrate1;
		q2 = m_divrate2;
	
		if(m_matu[curtN-1]>0) { // 현재 curtN 는 m_matuN <-- 만기해당..
			dailyIFrf = (double *)malloc((size_t)(m_matu[curtN-1])*sizeof(double));
			dailyIFdc = (double *)malloc((size_t)(m_matu[curtN-1])*sizeof(double));
		
			for(i=0; i<m_matu[curtN-1]; i++) {
				dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
				dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
			}
		} //if(m_matu[curtN-1]>0)
	
	 


		// time.. 계산 시작
		DCValue12 = knocknewV[0][0];
		DCValue1 = knock1newV[0][SM2];
		DCValue2 = knock2newV[SM1][0];


		for(matuLoop=m_matu[curtN-1]-1; matuLoop>=0; matuLoop--) {
			//time Adaptive mesh
			mid_n = IsinArr(matuLoop, m_matuN, m_matu); // 조기상환일 중에 하나면 1 ~ m_matuN-1 값 리턴
			if(mid_n > 0) { //
				curtN = mid_n;
				xlevel1 = m_xval1[curtN-1]/m_bval1;
				xlevel2 = m_xval2[curtN-1]/m_bval2;
 
				//n 영업일 평균 관련 처리
				evalq = m_evalN[curtN-1] - m_matu[curtN-1];
	
				if(evalq <=1 ) {	
					slevel1 = m_sval1/m_bval1;	
					slevel2 = m_sval2/m_bval2;
				}
				else {
					slevel1 = m_sval1/m_bval1 * (m_evalN[curtN-1] - (evalq-1));
					for(evali=1; evali<evalq; evali++)
						slevel1 = slevel1 + m_psval1[evali]/m_bval1;
					slevel1 = slevel1/m_evalN[curtN-1];

					slevel2 = m_sval2/m_bval2 * (m_evalN[curtN-1] - (evalq-1));
					for(evali=1; evali<evalq; evali++)
						slevel2 = slevel2 + m_psval2[evali]/m_bval2;
					slevel2 = slevel2/m_evalN[curtN-1];
				}	 

				m_endS1 = m_xval1[curtN-1]; 
				endSlevel1 = m_endS1/m_bval1;

				m_endS2 = m_xval2[curtN-1]; 
				endSlevel2 = m_endS2/m_bval2;

				//spreadXSlevel1 = m_spreadXS1[curtN-1]/m_bval1;
				//spreadXSlevel2 = m_spreadXS2[curtN-1]/m_bval2;

				endSidx1 = findinteridx(endSlevel1, S1, SM1);
				//spreadXSidx1 = findinteridx(spreadXSlevel1, S1, SM1);
				xidx1 = findinteridx(xlevel1, S1, SM1);		
			
				endSidx2 = findinteridx(endSlevel2, S2, SM2);
				//spreadXSidx2 = findinteridx(spreadXSlevel2, S2, SM2);
				xidx2 = findinteridx(xlevel2, S2, SM2);		
		 
	 
			} // 조기상환일 관련 setting


			if(matuLoop == m_matu[curtN-1]-1) 
				intraTN = m_meshTN1;		
			else  
				intraTN = m_meshTN2;		
			dt = (1.0/YFactor)/intraTN;

			r = dailyIFrf[matuLoop]; // 금리: One Day implied forward rate 를 이용
			rd = dailyIFdc[matuLoop];

			if(m_voltype == 0) {
				volskew1[0] = m_vol1[0][0];
				volskew2[0] = m_vol2[0][0];
			}
			else {
				curttimeVol1(volskew1, matuLoop);
				curttimeVol2(volskew2, matuLoop);
			}

			cor = curttimeCorr(matuLoop);
		 
			// X 생성		 
			for(j=1; j<SM2; j++) {  //y축 고정시켜서.. 
				if(m_voltype == 2) {				
					Srv2 = S2[j]*m_bval2 / m_volbval2;				
					v2 = interp1(Srv2, m_volSmness,volskew2,m_volSN,1,0);
				}
				else {
					v2 = volskew2[0];
				}

				for(i=1; i<SM1; i++) { //x축 풀기 위하여 행렬 만들기
					if(m_voltype == 2) {
						Srv1 = S1[i]*m_bval1 / m_volbval1;
						v1 = interp1(Srv1, m_volSmness, volskew1, m_volSN,1,0);
					}
					else {
						v1 = volskew1[0];
					}

					if(SmeshType1 == 0) { //균등 mesh
						knockXMAl[i][j] = -0.5*(v1*v1)*(i*i)*dt + (r-q1)*i*dt/2.0;// + 0.5*cor*v1*v2*i*j*dt/2.0;
						knockXMAd[i][j] = v1*v1*i*i*dt + 0.5*rd*dt + 1;
						knockXMAu[i][j] = -0.5*(v1*v1)*(i*i)*dt - (r-q1)*i*dt/2.0;// - 0.5*cor*v1*v2*i*j*dt/2.0;

						coeV1[i][j] = 0.5*cor*v1*v2*i*j*dt/4.0; //0.5*cor*v1*v2*i*j*dt/2.0;
					}
					else {						
						knockXMAl[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*(h1[i]+h1[i-1])) + (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]);// + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
						knockXMAd[i][j] = (v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*h1[i]) + 0.5*rd*dt + 1;
						knockXMAu[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i]*(h1[i]+h1[i-1])) - (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]);// - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);								
					
						coeV1[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*(h2[j]+h2[j-1]));	//0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
					}				
			
					if(m_ki1Flag == 0) {							
						XMAl[i][j] = knockXMAl[i][j];				
						XMAd[i][j] = knockXMAd[i][j];					
						XMAu[i][j] = knockXMAu[i][j];
					}

				} // 행렬 만들기					
				knockXMAl[0][j] = 0.0;				
				knockXMAd[0][j] = 1.0;				
				knockXMAu[0][j] = 0.0;		

				//경계 선형 조건
				if(SmeshType1 == 0) {
					knockXMAl[SM1-1][j] = knockXMAl[SM1-1][j] - knockXMAu[SM1-1][j];
					knockXMAd[SM1-1][j] = knockXMAd[SM1-1][j] + 2*knockXMAu[SM1-1][j];
				}
				else {			 
					knockXMAl[SM1-1][j] = knockXMAl[SM1-1][j] - knockXMAu[SM1-1][j]*h1[SM1-1]/h1[SM1-2];
					knockXMAd[SM1-1][j] = knockXMAd[SM1-1][j] + knockXMAu[SM1-1][j]*(h1[SM1-1]+h1[SM1-2])/h1[SM1-2];			
				}

				if(m_ki1Flag == 0) {
					XMAl[SM1-1][j] = knockXMAl[SM1-1][j];			
					XMAd[SM1-1][j] = knockXMAd[SM1-1][j];

					XMAl[startSidx1][j] = 0;		
					XMAd[startSidx1][j] = 1;			
					XMAu[startSidx1][j] = 0;	
						
					//LU decomposition			
					for(i=startSidx1+1; i<=SM1-1; i++) {				
						XMAl[i][j] = XMAl[i][j]/XMAd[i-1][j];				
						XMAd[i][j] = XMAd[i][j] - XMAl[i][j]*XMAu[i-1][j];			
					}			
				}

				//LU decomposition
				for(i=1; i<=SM1-1; i++) {
					knockXMAl[i][j] = knockXMAl[i][j]/knockXMAd[i-1][j];
					knockXMAd[i][j] = knockXMAd[i][j] - knockXMAl[i][j]*knockXMAu[i-1][j];
				}				
			} //X 행렬 끝
	 
			// Y 생성		 
			for(i=1; i<SM1; i++) {  //y축 고정시켜서.. 
				if(m_voltype == 2) {				
					Srv1 = S1[i]*m_bval1 / m_volbval1;				
					v1 = interp1(Srv1, m_volSmness,volskew1,m_volSN,1,0);
				}
				else {
					v1 = volskew1[0];
				}

				for(j=1; j<SM2; j++) { //x축 풀기 위하여 행렬 만들기
					if(m_voltype == 2) {
						Srv2 = S2[j]*m_bval2 / m_volbval2;
						v2 = interp1(Srv2, m_volSmness, volskew2, m_volSN,1,0);
					}
					else {
						v2 = volskew2[0];
					}

					if(SmeshType2 == 0) { //균등 mesh
						knockYMAl[i][j] = -0.5*(v2*v2)*(j*j)*dt + (r-q2)*j*dt/2.0;// + 0.5*cor*v1*v2*i*j*dt/2.0;
						knockYMAd[i][j] = v2*v2*j*j*dt + 0.5*rd*dt + 1;
						knockYMAu[i][j] = -0.5*(v2*v2)*(j*j)*dt - (r-q2)*j*dt/2.0;// - 0.5*cor*v1*v2*i*j*dt/2.0;

						coeV2[i][j] = 0.5*cor*v1*v2*i*j*dt/4.0;  //0.5*cor*v1*v2*i*j*dt/2.0;
					}
					else {						
						knockYMAl[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*(h2[j]+h2[j-1])) + (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]);// + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);
						knockYMAd[i][j] = (v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*h2[j]) + 0.5*rd*dt + 1;
						knockYMAu[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j]*(h2[j]+h2[j-1])) - (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]);// - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);							
					
						coeV2[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*(h1[i]+h1[i-1])); //0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);	
					}				
			
					if(m_ki2Flag == 0) {							
						YMAl[i][j] = knockYMAl[i][j];				
						YMAd[i][j] = knockYMAd[i][j];					
						YMAu[i][j] = knockYMAu[i][j];
					}

				} // 행렬 만들기					
				knockYMAl[i][0] = 0.0;				
				knockYMAd[i][0] = 1.0;				
				knockYMAu[i][0] = 0.0;	

				//경계 선형 조건
				if(SmeshType2 == 0) {
					knockYMAl[i][SM2-1] = knockYMAl[i][SM2-1] - knockYMAu[i][SM2-1];
					knockYMAd[i][SM2-1] = knockYMAd[i][SM2-1] + 2*knockYMAu[i][SM2-1];
				}
				else {			 
					knockYMAl[i][SM2-1] = knockYMAl[i][SM2-1] - knockYMAu[i][SM2-1]*h2[SM2-1]/h2[SM2-2];
					knockYMAd[i][SM2-1] = knockYMAd[i][SM2-1] + knockYMAu[i][SM2-1]*(h2[SM2-1]+h2[SM2-2])/h2[SM2-2];			
				}

				if(m_ki2Flag == 0) {
					YMAl[i][SM2-1] = knockYMAl[i][SM2-1];			
					YMAd[i][SM2-1] = knockYMAd[i][SM2-1];

					YMAl[i][startSidx2] = 0;		
					YMAd[i][startSidx2] = 1;			
					YMAu[i][startSidx2] = 0;	
						
					//LU decomposition			
					for(j=startSidx2+1; j<=SM2-1; j++) {				
						YMAl[i][j] = YMAl[i][j]/YMAd[i][j-1];				
						YMAd[i][j] = YMAd[i][j] - YMAl[i][j]*YMAu[i][j-1];			
					}			
				}

				//LU decomposition
				for(j=1; j<=SM2-1; j++) {
					knockYMAl[i][j] = knockYMAl[i][j]/knockYMAd[i][j-1];									
					knockYMAd[i][j] = knockYMAd[i][j] - knockYMAl[i][j]*knockYMAu[i][j-1];
				}				
			} //Y 행렬 끝		 
	   
		
			for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {
				DCValue12 = DCValue12 * 1/(1+rd/(YFactor*intraTN));
				DCValue1 = DCValue1 * 1/(1+rd/(YFactor*intraTN));
				DCValue2 = DCValue2 * 1/(1+rd/(YFactor*intraTN));

				//memcpy(knockoldV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
				for(i=0; i<=SM1; i++)			
					memcpy(knockoldV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
			
				if(m_ki2Flag == 0)
					for(i=0; i<=SM1; i++)				
						memcpy(knock1oldV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));

				if(m_ki1Flag == 0)
					for(i=0; i<=SM1; i++)				
						memcpy(knock2oldV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));

				if(m_ki1Flag == 0 && m_ki2Flag == 0)
					for(i=0; i<=SM1; i++)				
						memcpy(oldV[i],newV[i],(size_t)((SM2+1)*sizeof(double)));

				/*
				for(i=0; i<=SM1; i++)
					for(j=0; j<=SM2; j++) {
						knockoldV[i][j] = knocknewV[i][j];
						if(m_kiFlag != 1) 
							oldV[i][j] = newV[i][j];
					}
				*/
			
				IsMidCompu = 0;

				if(intraLoop == 1 && mid_n > 0) { // 조기 상환일 중에 하나
					if(matuLoop == 0 && m_ctime > 15 && (slevel1 < xlevel1 || slevel2 < xlevel2))   // 오늘이 조기상환일 이면서 종가 이후 상환 안된 때..
						IsMidCompu = 0; // 조기상환 안된 next 헤지를 위해서..														
					else
						IsMidCompu = 1;
				}

				if(IsMidCompu == 0) { //조기 상환일이 아닌 평시..형태로 풀기 (조기상환일에 상환 안된 경우도 해당)				
					for(i=0; i<=SM1; i++) 
						tU[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;

					for(j=1; j<SM2; j++) { //y축 고정 풀기

						for(i=1; i<SM1; i++) 
							vecV1[i] = knockoldV[i][j] + coeV1[i][j]*(-knockoldV[i+1][j-1]+knockoldV[i-1][j-1] +knockoldV[i+1][j+1]-knockoldV[i-1][j+1]); // herexxx

						tU[0][j] = knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
						for(i=1; i<SM1; i++)
							tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];
						tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
						for(i=SM1-2; i>=0; i--) 				
							tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
										
						if(SmeshType1 == 0)				
							tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
						else
							tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				

					}// for(j=1

					for(i=0; i<=SM1; i++) {				
						if(SmeshType2 == 0) 					
							tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
						else					
							tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
					} // 일차 tU 완결 

					//////////////////////////////////////////////////////// newV 구하기!!
					for(j=0; j<=SM2; j++)
						knocknewV[0][j]=knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
		 
					for(i=1; i<SM1; i++) {  //x축 고정시켜서.. 				 

						for(j=1; j<SM2; j++) 				
							vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

						knocknewV[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
						for(j=1; j<SM2; j++) 
							knocknewV[i][j] = vecV2[j] - knockYMAl[i][j]*knocknewV[i][j-1];					
						knocknewV[i][SM2-1] = knocknewV[i][SM2-1]/knockYMAd[i][SM2-1];
						for(j=SM2-2; j>=0; j--) 
							knocknewV[i][j] = (knocknewV[i][j] - knockYMAu[i][j]*knocknewV[i][j+1])/knockYMAd[i][j];					
				
						if(SmeshType2 == 0)				
							knocknewV[i][SM2] = 2*knocknewV[i][SM2-1] - knocknewV[i][SM2-2];
						else
							knocknewV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knocknewV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knocknewV[i][SM2-2];
				
					} // for(i x축.. 고정 다 품

					for(j=0; j<=SM2; j++) {
						if(SmeshType1 == 0) 
							knocknewV[SM1][j] = 2*knocknewV[SM1-1][j] - knocknewV[SM1-2][j];
						else
							knocknewV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knocknewV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knocknewV[SM1-2][j];
					} // knocknewV 완결

					//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knocknewV[0][j] = 2*knocknewV[1][j] - knocknewV[2][j];
								 
							}
							else {
								knocknewV[0][j] = (h1[0]+h1[1])/h1[1] * knocknewV[1][j] - h1[0]/h1[1] * knocknewV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knocknewV[i][0] = 2*knocknewV[i][1] - knocknewV[i][2];
								 
							}
							else {
								knocknewV[i][0] = (h1[0]+h1[1])/h1[1] * knocknewV[i][1] - h1[0]/h1[1] * knocknewV[i][2];
								 
							}
						} // knocknewV 완결



	//////////////////////////// KI1				//herework!!					
					if(m_ki2Flag == 0) {						 
						 
						for(i=0; i<=SM1; i++)
							tU[i][startSidx2] = knocknewV[i][startSidx2];

						for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기						
							for(i=1; i<SM1; i++) 
								vecV1[i] = knock1oldV[i][j] + coeV1[i][j]*(-knock1oldV[i+1][j-1]+knock1oldV[i-1][j-1] +knock1oldV[i+1][j+1]-knock1oldV[i-1][j+1]); //herexxx
                      
							tU[0][j] = knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1;
							for(i=1; i<SM1; i++)
								tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];

							tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
							for(i=SM1-2; i>=0; i--) 											
								tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
													
							if(SmeshType1 == 0) {		
								//tU[0][j] = 2*tU[1][j] - tU[2][j];
								tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
							}
							else {
								//tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
								tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
							}
						 			

						}// for(j=1

						for(i=0; i<=SM1; i++) {				
							if(SmeshType2 == 0) 			
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							else					
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
						} // 일차 tU 완결  
						
						//호
						for(j=startSidx2+1; j<=SM2; j++) {
							if(SmeshType1 == 0) 			
								tU[0][j] = 2*tU[1][j] - tU[2][j];				
							else					
								tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
						}

					//////////////////////////////////////////////////////// newV 구하기!!
						for(j=startSidx2+1; j<=SM2; j++)
							knock1newV[0][j] = tU[0][j]; //knock1oldV[0][j]* 1/(1+rd/(252*intraTN)); //DCValue1; //0.0;
		 
						for(i=1; i<SM1; i++) {  //x축 고정시켜서.. 				 

							for(j=startSidx2+1; j<SM2; j++) 				
								vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

							knock1newV[i][startSidx2] =  knocknewV[i][startSidx2]; 
							for(j=startSidx2+1; j<SM2; j++) 								
								knock1newV[i][j] = vecV2[j] - YMAl[i][j]*knock1newV[i][j-1];					

							knock1newV[i][SM2-1] = knock1newV[i][SM2-1]/YMAd[i][SM2-1];
							for(j=SM2-2; j>=startSidx2; j--) 						
								knock1newV[i][j] = (knock1newV[i][j] - YMAu[i][j]*knock1newV[i][j+1])/YMAd[i][j];			
				
							if(SmeshType2 == 0)				
								knock1newV[i][SM2] = 2*knock1newV[i][SM2-1] - knock1newV[i][SM2-2];
							else
								knock1newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock1newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock1newV[i][SM2-2];
				
						} // for(i x축.. 고정 다 품

						//호
						for(j=startSidx2+1; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock1newV[0][j] = 2*knock1newV[1][j] - knock1newV[2][j];
								knock1newV[SM1][j] = 2*knock1newV[SM1-1][j] - knock1newV[SM1-2][j];
							}
							else {
								knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
								knock1newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock1newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock1newV[SM1-2][j];
							}
						} // knocknewV 완결

						
					//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock1newV[0][j] = 2*knock1newV[1][j] - knock1newV[2][j];
								 
							}
							else {
								knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock1newV[i][0] = 2*knock1newV[i][1] - knock1newV[i][2];
								 
							}
							else {
								knock1newV[i][0] = (h1[0]+h1[1])/h1[1] * knock1newV[i][1] - h1[0]/h1[1] * knock1newV[i][2];
								 
							}
						} // knocknewV 완결



						for(i=0; i<=SM1; i++)
							for(j=0; j<=startSidx2; j++)  									
								knock1newV[i][j] = knocknewV[i][j];
 					 
					} 
	//////////////////////////////// KI1 end
	 
 	//////////////////////////// KI2				//herework!!					
					if(m_ki1Flag == 0) {						 
						 
						for(i=startSidx1+1; i<=SM1; i++)
							tU[i][0] = knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); //DCValue2;

						for(j=1; j<SM2; j++) { //y축 고정 풀기						
							for(i=startSidx1+1; i<SM1; i++) 
								vecV1[i] = knock2oldV[i][j] + coeV1[i][j]*(-knock2oldV[i+1][j-1]+knock2oldV[i-1][j-1] +knock2oldV[i+1][j+1]-knock2oldV[i-1][j+1]); //herexxx
                      
							tU[startSidx1][j] = knocknewV[startSidx1][j];
							for(i=startSidx1+1; i<SM1; i++)
								tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];

							tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
							for(i=SM1-2; i>=startSidx1; i--) 											
								tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
													
							if(SmeshType1 == 0)									
								tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];							
							else								
								tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
							
						 			

						}// for(j=1

						for(i=startSidx1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								//tU[i][0] = 2*tU[i][1] - tU[i][2];
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							}
							else  {
								//tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
							}
						} // 일차 tU 완결  

						//호
						for(i=startSidx1+1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								tU[i][0] = 2*tU[i][1] - tU[i][2];												
							}
							else  {
								tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];											
							}
						} // 일차 tU 완결  
						
					//////////////////////////////////////////////////////// newV 구하기!!
						for(j=0; j<=SM2; j++)
							knock2newV[startSidx1][j]=knocknewV[startSidx1][j];
		 
						for(i=startSidx1+1; i<SM1; i++) {  //x축 고정시켜서.. 				 

							for(j=1; j<SM2; j++) 				
								vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

							knock2newV[i][0] = knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); // DCValue2; 
							for(j=1; j<SM2; j++) 								
								knock2newV[i][j] = vecV2[j] - knockYMAl[i][j]*knock2newV[i][j-1];					

							knock2newV[i][SM2-1] = knock2newV[i][SM2-1]/knockYMAd[i][SM2-1];
							for(j=SM2-2; j>=0; j--) 						
								knock2newV[i][j] = (knock2newV[i][j] - knockYMAu[i][j]*knock2newV[i][j+1])/knockYMAd[i][j];			
				
							if(SmeshType2 == 0) {
								//knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];
								knock2newV[i][SM2] = 2*knock2newV[i][SM2-1] - knock2newV[i][SM2-2];
							}
							else {
								//knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];
								knock2newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock2newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock2newV[i][SM2-2];
							}
				
						} // for(i x축.. 고정 다 품

						for(j=1; j<=SM2; j++) {
							if(SmeshType1 == 0) 
								knock2newV[SM1][j] = 2*knock2newV[SM1-1][j] - knock2newV[SM1-2][j];
							else
								knock2newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock2newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock2newV[SM1-2][j];
						} // knocknewV 완결

						//호
						for(i=startSidx1+1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];												
							}
							else  {
								knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];											
							}
						} // 일차 tU 완결  

											
					//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock2newV[0][j] = 2*knock2newV[1][j] - knock2newV[2][j];
								 
							}
							else {
								knock2newV[0][j] = (h1[0]+h1[1])/h1[1] * knock2newV[1][j] - h1[0]/h1[1] * knock2newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];
								 
							}
							else {
								knock2newV[i][0] = (h1[0]+h1[1])/h1[1] * knock2newV[i][1] - h1[0]/h1[1] * knock2newV[i][2];
								 
							}
						} // knocknewV 완결




						for(i=0; i<=startSidx1; i++)
							for(j=0; j<=SM2; j++)  									
								knock2newV[i][j] = knocknewV[i][j];
 					 
					} 
	//////////////////////////////// KI2 end
	 
					// no KI 영역 풀기 
					if(m_ki1Flag == 0 && m_ki2Flag == 0) {				
					 
						for(i=startSidx1; i<=SM1; i++)
							tU[i][startSidx2] = knock2newV[i][startSidx2];

						for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기						
							for(i=startSidx1+1; i<SM1; i++) 
								vecV1[i] = oldV[i][j] + coeV1[i][j]*(-oldV[i+1][j-1]+oldV[i-1][j-1] +oldV[i+1][j+1]-oldV[i-1][j+1]); //herexxx
                        
							tU[startSidx1][j] = knock1newV[startSidx1][j];
							for(i=startSidx1+1; i<SM1; i++)
								tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
							tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
							for(i=SM1-2; i>=startSidx1; i--) 				
								tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
										
							if(SmeshType1 == 0)				
								tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
							else
								tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				

						}// for(j=1

						for(i=startSidx1; i<=SM1; i++) {				
							if(SmeshType2 == 0) 					
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							else					
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
						} // 일차 tU 완결 

						//////////////////////////////////////////////////////// newV 구하기!!
						for(j=startSidx2; j<=SM2; j++)
							newV[startSidx1][j]=knock1newV[startSidx1][j]; // D.C 적용
		 
						for(i=startSidx1+1; i<SM1; i++) {  //y축 고정시켜서.. 				 

							for(j=startSidx2+1; j<SM2; j++) 				
								vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);	//herexxx				   				 		 		

							newV[i][startSidx2] = knock2newV[i][startSidx2];
							for(j=startSidx2+1; j<SM2; j++) 
								newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];					
							newV[i][SM2-1] = newV[i][SM2-1]/YMAd[i][SM2-1];
							for(j=SM2-2; j>=startSidx2; j--) 
								newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];					
				
							if(SmeshType2 == 0)				
								newV[i][SM2] = 2*newV[i][SM2-1] - newV[i][SM2-2];
							else
								newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * newV[i][SM2-2];
				
						} // for(j y축.. 고정 다 품

						for(j=startSidx2; j<=SM2; j++) {
							if(SmeshType1 == 0) 
								newV[SM1][j] = 2*newV[SM1-1][j] - newV[SM1-2][j];
							else
								newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * newV[SM1-2][j];
						} // newV 완결

						for(i=0; i<=startSidx1; i++)
							for(j=0; j<=startSidx2; j++) 								
									newV[i][j] = knocknewV[i][j];  // KI12

						for(i=0; i<=startSidx1; i++)
							for(j=startSidx2+1; j<=SM2; j++) 								
									newV[i][j] = knock1newV[i][j]; // KI1

						for(i=startSidx1+1; i<=SM1; i++)
							for(j=0; j<=startSidx2; j++) 								
									newV[i][j] = knock2newV[i][j]; // KI2
							
 					 
					}// if(m_kiFlag != 1)

				} //if(IsMidCompu == 0)
				else { // 조기상환일 계산 

					for(i=0; i<=SM1; i++) 
						tU[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;

					for(j=1; j<SM2; j++) { //y축 고정 풀기

						if(j<endSidx2) { // 행사 안되는 영역  if(j<endSidx2)
							for(i=1; i<SM1; i++) 
								vecV1[i] = knockoldV[i][j] + coeV1[i][j]*(-knockoldV[i+1][j-1]+knockoldV[i-1][j-1] + knockoldV[i+1][j+1] - knockoldV[i-1][j+1]); //herexxx

							tU[0][j] =knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(i=1; i<SM1; i++)
								tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];
							tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
							for(i=SM1-2; i>=0; i--) 				
								tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
										
							if(SmeshType1 == 0)				
								tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
							else
								tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				
						}
						else { // S2 행사가 이상 영역
							for(i=1; i<endSidx1; i++) // 행사 안되는 영역 
								vecV1[i] = knockoldV[i][j] + coeV1[i][j]*(-knockoldV[i+1][j-1]+knockoldV[i-1][j-1] + knockoldV[i+1][j+1] - knockoldV[i-1][j+1]); //herexxx
						
							tU[0][j] = knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(i=1; i<endSidx1; i++) 
								tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];
							for(i=endSidx1; i<=SM1; i++)
								tU[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
							for(i=endSidx1-1; i>=0; i--)
								tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
						}
					}// for(j=1

					for(i=0; i<=SM1; i++) {				
						if(SmeshType2 == 0) 					
							tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
						else					
							tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
					} // 일차 tU 완결 					 

					//////////////////////////////////////////////////////// newV 구하기!!
					for(j=0; j<=SM2; j++)
						knocknewV[0][j]=knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0; // D.C 적용
		 
					for(i=1; i<SM1; i++) {  //y축 고정시켜서.. 				 
						if(i<endSidx1) {
							for(j=1; j<SM2; j++) 				
								vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);  //herexxx					   				 		 		 

							knocknewV[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(j=1; j<SM2; j++) 
								knocknewV[i][j] = vecV2[j] - knockYMAl[i][j]*knocknewV[i][j-1];					
							knocknewV[i][SM2-1] = knocknewV[i][SM2-1]/knockYMAd[i][SM2-1];
							for(j=SM2-2; j>=0; j--) 
								knocknewV[i][j] = (knocknewV[i][j] - knockYMAu[i][j]*knocknewV[i][j+1])/knockYMAd[i][j];					
				
							if(SmeshType2 == 0)				
								knocknewV[i][SM2] = 2*knocknewV[i][SM2-1] - knocknewV[i][SM2-2];
							else
								knocknewV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knocknewV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knocknewV[i][SM2-2];
						}
						else {
							for(j=1; j<endSidx2; j++)
								vecV2[j] = tU[i][j] + coeV2[i][j]*(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]); //herexxx

							knocknewV[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(j=1; j<endSidx2; j++)
								knocknewV[i][j] = vecV2[j] - knockYMAl[i][j]*knocknewV[i][j-1];
							for(j=endSidx2; j<=SM2; j++)
								knocknewV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
							for(j=endSidx2-1; j>=0; j--) 
								knocknewV[i][j] = (knocknewV[i][j] - knockYMAu[i][j]*knocknewV[i][j+1])/knockYMAd[i][j];
						}				
					} // for(j y축.. 고정 다 품

					for(j=0; j<=SM2; j++) {
						if(SmeshType1 == 0) 
							knocknewV[SM1][j] = 2*knocknewV[SM1-1][j] - knocknewV[SM1-2][j];
						else
							knocknewV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knocknewV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knocknewV[SM1-2][j];
					} // knocknewV 완결
 			
														
					//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knocknewV[0][j] = 2*knocknewV[1][j] - knocknewV[2][j];
								 
							}
							else {
								knocknewV[0][j] = (h1[0]+h1[1])/h1[1] * knocknewV[1][j] - h1[0]/h1[1] * knocknewV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knocknewV[i][0] = 2*knocknewV[i][1] - knocknewV[i][2];
								 
							}
							else {
								knocknewV[i][0] = (h1[0]+h1[1])/h1[1] * knocknewV[i][1] - h1[0]/h1[1] * knocknewV[i][2];
								 
							}
						} // knocknewV 완결



					if(m_spreadXRate[curtN-1] > 0) {// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
						double spreadvalue;
						double   spmax;
		
						//for(i=spreadXSidx1; i<=SM1; i++) {
							//for(j=spreadXSidx2; j<=SM2; j++) {
						for(i=0; i<=SM1; i++) {
							for(j=0; j<=SM2; j++) {
								if(i<=xidx1 || j<=xidx2) {
	

						if(S1[i]<=S2[j]) {							 
							spmax = S1[xidx1];																				 
						}
						else {						 
							spmax = S2[xidx2];																			 
						}						 
						spreadvalue = AMOUNT * ( m_spreadXRate[curtN-1] * (MIN(S1[i],S2[j])-spmax) + 1+m_Cpn[curtN-1]);
									knocknewV[i][j] = MAX(knocknewV[i][j], spreadvalue);								
								}
							}
						}

						/*
											//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						double spvalueR;

						spvalueR =  AMOUNT * (1+m_Cpn[curtN-1]);
						for(i=xidx1; i<=SM1; i++) {
							maxindexR = 0;
							maxvalueR = knocknewV[i][0];
							for(j=1; j<=xidx2; j++) {
								if(maxvalueR<knocknewV[i][j]) {
									maxindexR = j;
									maxvalueR = knocknewV[i][j];
								}
							}
							if(maxindexR < xidx2) {
								tempindex = 2*xidx2 - spreadXSidx2;
								for(j=maxindexR+1; j<tempindex; j++) {
									spreadvalue = (maxvalueR-spvalueR)/(S2[maxindexR]-S2[tempindex]) * (S2[j]-S2[tempindex]) + spvalueR;
									knocknewV[i][j] = MAX(knocknewV[i][j],spreadvalue);
								}
							}
						}

						for(j=xidx2; j<=SM2; j++) {														
							maxindexR = 0;
							maxvalueR = knocknewV[0][j];
							for(i=1; i<=xidx1; i++) {
								if(maxvalueR<knocknewV[i][j]) {
									maxindexR = i;
									maxvalueR = knocknewV[i][j];
								}
							}
							if(maxindexR < xidx1) {
								tempindex = 2*xidx1 - spreadXSidx1;
								for(i=maxindexR+1; i<tempindex; i++) {
									spreadvalue = (maxvalueR-spvalueR)/(S1[maxindexR]-S1[tempindex]) * (S1[i]-S1[tempindex]) + spvalueR;
									knocknewV[i][j] = MAX(knocknewV[i][j],spreadvalue);
								}
							}
						}
						*/
							
					
					
					} // X spread
 	

						//////////////////////////// KI1				//herework!!					
					if(m_ki2Flag == 0) {						 
						 
						for(i=0; i<=SM1; i++)
							tU[i][startSidx2] = knocknewV[i][startSidx2];

						for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기	
							if(j<endSidx2) {
								for(i=1; i<SM1; i++) 
									vecV1[i] = knock1oldV[i][j] + coeV1[i][j]*(-knock1oldV[i+1][j-1]+knock1oldV[i-1][j-1] +knock1oldV[i+1][j+1]-knock1oldV[i-1][j+1]); //herexxx
                      
								tU[0][j] = knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1;
								for(i=1; i<SM1; i++)
									tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];

								tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
								for(i=SM1-2; i>=0; i--) 											
									tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
													
								if(SmeshType1 == 0)	{
									//tU[0][j] = 2*tU[1][j] - tU[2][j];
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
								}
								else {
									//tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
								}
							}
							else {															
								for(i=1; i<endSidx1; i++) 
									vecV1[i] = knock1oldV[i][j] + coeV1[i][j]*(-knock1oldV[i+1][j-1]+knock1oldV[i-1][j-1] +knock1oldV[i+1][j+1]-knock1oldV[i-1][j+1]); //herexxx
                      
								tU[0][j] = knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1;
								for(i=1; i<endSidx1; i++)
									tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];

								for(i=endSidx1; i<=SM1; i++)								
									tU[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(i=endSidx1-1; i>=0; i--) 											
									tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];

								/*					
								if(SmeshType1 == 0) 
									tU[0][j] = 2*tU[1][j] - tU[2][j]; 
								else 
									tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
									*/  
							}						 
						}// for(j=1

						for(i=0; i<=SM1; i++) {				
							if(SmeshType2 == 0) 					
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							else					
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
						} // 일차 tU 완결  
						
						//호
						for(j=startSidx2+1; j<=SM2; j++) {
							if(SmeshType1 == 0) 			
								tU[0][j] = 2*tU[1][j] - tU[2][j];				
							else					
								tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
						}
					//////////////////////////////////////////////////////// newV 구하기!!
						for(j=startSidx2+1; j<=SM2; j++)
							knock1newV[0][j]=knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1; //0.0;						 	

						for(i=1; i<SM1; i++) {  //x축 고정시켜서.. 				 
							if(i<endSidx1) {
								for(j=startSidx2+1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock1newV[i][startSidx2] =  knocknewV[i][startSidx2]; 
								for(j=startSidx2+1; j<SM2; j++) 								
									knock1newV[i][j] = vecV2[j] - YMAl[i][j]*knock1newV[i][j-1];					

								knock1newV[i][SM2-1] = knock1newV[i][SM2-1]/YMAd[i][SM2-1];
								for(j=SM2-2; j>=startSidx2; j--) 						
									knock1newV[i][j] = (knock1newV[i][j] - YMAu[i][j]*knock1newV[i][j+1])/YMAd[i][j];			
				
								if(SmeshType2 == 0)				
									knock1newV[i][SM2] = 2*knock1newV[i][SM2-1] - knock1newV[i][SM2-2];
								else
									knock1newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock1newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock1newV[i][SM2-2];								
							}
							else {
								for(j=startSidx2+1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock1newV[i][startSidx2] =  knocknewV[i][startSidx2]; 
								for(j=startSidx2+1; j<endSidx2; j++) 								
									knock1newV[i][j] = vecV2[j] - YMAl[i][j]*knock1newV[i][j-1];					

								for(j=endSidx2; j<=SM2; j++)								
									knock1newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(j=endSidx2-1; j>=startSidx2; j--) 						
									knock1newV[i][j] = (knock1newV[i][j] - YMAu[i][j]*knock1newV[i][j+1])/YMAd[i][j];
							}
						
				
						} // for(i x축.. 고정 다 품

						for(j=startSidx2+1; j<=SM2; j++) {
							if(SmeshType1 == 0) {
								knock1newV[0][j] = 2 * knock1newV[1][j] - knock1newV[2][j];
								knock1newV[SM1][j] = 2*knock1newV[SM1-1][j] - knock1newV[SM1-2][j];
							}
							else {
								knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
								knock1newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock1newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock1newV[SM1-2][j];
							}
						} // knocknewV 완결

								//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock1newV[0][j] = 2*knock1newV[1][j] - knock1newV[2][j];
								 
							}
							else {
								knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock1newV[i][0] = 2*knock1newV[i][1] - knock1newV[i][2];
								 
							}
							else {
								knock1newV[i][0] = (h1[0]+h1[1])/h1[1] * knock1newV[i][1] - h1[0]/h1[1] * knock1newV[i][2];
								 
							}
						} // knocknewV 완결


						for(i=0; i<=SM1; i++)
							for(j=0; j<=startSidx2; j++)  									
								knock1newV[i][j] = knocknewV[i][j];

									
						if(m_spreadXRate[curtN-1] > 0) { // && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
							double spreadvalue;
							double   spmax;
		
							
							//for(i=spreadXSidx1; i<=SM1; i++) {
								//for(j=spreadXSidx2; j<=SM2; j++) {
							for(i=0; i<=SM1; i++) {
								for(j=0; j<=SM2; j++) {
					
									if(i<=xidx1 || j<=xidx2) {

						if(S1[i]<=S2[j]) {							 
							spmax = S1[xidx1];																				 
						}
						else {						 
							spmax = S2[xidx2];																			 
						}						 
						spreadvalue = AMOUNT * ( m_spreadXRate[curtN-1] * (MIN(S1[i],S2[j])-spmax) + 1+m_Cpn[curtN-1]);
										knock1newV[i][j] = MAX(knock1newV[i][j], spreadvalue);								
									}
								}
							}
						
							/*
																		//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						double spvalueR;

						spvalueR =  AMOUNT * (1+m_Cpn[curtN-1]);
						for(i=xidx1; i<=SM1; i++) {
							maxindexR = 0;
							maxvalueR = knock1newV[i][0];
							for(j=1; j<=xidx2; j++) {
								if(maxvalueR<knock1newV[i][j]) {
									maxindexR = j;
									maxvalueR = knock1newV[i][j];
								}
							}
							if(maxindexR < xidx2) {
								tempindex = 2*xidx2 - spreadXSidx2;
								for(j=maxindexR+1; j<tempindex; j++) {
									spreadvalue = (maxvalueR-spvalueR)/(S2[maxindexR]-S2[tempindex]) * (S2[j]-S2[tempindex]) + spvalueR;
									knock1newV[i][j] = MAX(knock1newV[i][j],spreadvalue);
								}
							}
						}

						for(j=xidx2; j<=SM2; j++) {														
							maxindexR = 0;
							maxvalueR = knock1newV[0][j];
							for(i=1; i<=xidx1; i++) {
								if(maxvalueR<knock1newV[i][j]) {
									maxindexR = i;
									maxvalueR = knock1newV[i][j];
								}
							}
							if(maxindexR < xidx1) {
								tempindex = 2*xidx1 - spreadXSidx1;
								for(i=maxindexR+1; i<tempindex; i++) {
									spreadvalue = (maxvalueR-spvalueR)/(S1[maxindexR]-S1[tempindex]) * (S1[i]-S1[tempindex]) + spvalueR;
									knock1newV[i][j] = MAX(knock1newV[i][j],spreadvalue);
								}
							}
						}
							*/
						
						} // X spread
 					 
					} 
	//////////////////////////////// KI1 end
	 
 	//////////////////////////// KI2				//herework!!					
					if(m_ki1Flag == 0) {						 
						 
						for(i=startSidx1+1; i<=SM1; i++)
							tU[i][0] = knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); //DCValue2;	 						  
						 				 
						for(j=1; j<SM2; j++) { //y축 고정 풀기	
							if(j<endSidx2) {
								for(i=startSidx1+1; i<SM1; i++) 
									vecV1[i] = knock2oldV[i][j] + coeV1[i][j]*(-knock2oldV[i+1][j-1]+knock2oldV[i-1][j-1] +knock2oldV[i+1][j+1]-knock2oldV[i-1][j+1]); //herexxx
                      
								tU[startSidx1][j] = knocknewV[startSidx1][j];
								for(i=startSidx1+1; i<SM1; i++)
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];

								tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
								for(i=SM1-2; i>=startSidx1; i--) 											
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
													
								if(SmeshType1 == 0)				
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
								else
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
							}
							else {															
								for(i=startSidx1+1; i<endSidx1; i++) 
									vecV1[i] = knock2oldV[i][j] + coeV1[i][j]*(-knock2oldV[i+1][j-1]+knock2oldV[i-1][j-1] +knock2oldV[i+1][j+1]-knock2oldV[i-1][j+1]); //herexxx
                      
								tU[startSidx1][j] = knocknewV[startSidx1][j];
								for(i=startSidx1+1; i<endSidx1; i++)
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];

								for(i=endSidx1; i<=SM1; i++)								
									tU[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(i=endSidx1-1; i>=startSidx1; i--) 											
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
							}

						}// for(j=1

						for(i=startSidx1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								//tU[i][0] = 2 * tU[i][1] - tU[i][2];
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							}
							else  {
								//tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
							}
						} // 일차 tU 완결  
						
						//호
						for(i=startSidx1+1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								tU[i][0] = 2*tU[i][1] - tU[i][2];												
							}
							else  {
								tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];											
							}
						} // 일차 tU 완결  
						
					//////////////////////////////////////////////////////// newV 구하기!!
						for(j=0; j<=SM2; j++)
							knock2newV[startSidx1][j]=knocknewV[startSidx1][j];
		 
						for(i=startSidx1+1; i<SM1; i++) {  //x축 고정시켜서.. 				 
							if(i<endSidx1) {
								for(j=1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock2newV[i][0] =  knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); //DCValue2; 
								for(j=1; j<SM2; j++) 								
									knock2newV[i][j] = vecV2[j] - knockYMAl[i][j]*knock2newV[i][j-1];					

								knock2newV[i][SM2-1] = knock2newV[i][SM2-1]/knockYMAd[i][SM2-1];
								for(j=SM2-2; j>=0; j--) 						
									knock2newV[i][j] = (knock2newV[i][j] - knockYMAu[i][j]*knock2newV[i][j+1])/knockYMAd[i][j];			
				
								if(SmeshType2 == 0) {
									//knock2newV[i][0] = 2 * knock2newV[i][1] - knock2newV[i][2];
									knock2newV[i][SM2] = 2*knock2newV[i][SM2-1] - knock2newV[i][SM2-2];
								}
								else{
									//knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];
									knock2newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock2newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock2newV[i][SM2-2];
								}
							}
							else {															
								for(j=1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock2newV[i][0] =  knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); //DCValue2; 
								for(j=1; j<endSidx2; j++) 								
									knock2newV[i][j] = vecV2[j] - knockYMAl[i][j]*knock2newV[i][j-1];					
								
								for(j=endSidx2; j<=SM2; j++)								
									knock2newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(j=endSidx2-1; j>=0; j--) 						
									knock2newV[i][j] = (knock2newV[i][j] - knockYMAu[i][j]*knock2newV[i][j+1])/knockYMAd[i][j];			
								
								/*
								if(SmeshType2 == 0)  
									knock2newV[i][0] = 2 * knock2newV[i][1] - knock2newV[i][2];
								 else
									 knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];
									 */
							}

						} // for(i x축.. 고정 다 품

						for(j=1; j<=SM2; j++) {
							if(SmeshType1 == 0) 
								knock2newV[SM1][j] = 2*knock2newV[SM1-1][j] - knock2newV[SM1-2][j];
							else
								knock2newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock2newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock2newV[SM1-2][j];
						} // knocknewV 완결

						
						//호
						for(i=startSidx1+1; i<=SM1; i++) {				
							if(SmeshType2 == 0) {
								knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];												
							}
							else  {
								knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];											
							}
						} // 일차 tU 완결  

						
								//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock2newV[0][j] = 2*knock2newV[1][j] - knock2newV[2][j];
								 
							}
							else {
								knock2newV[0][j] = (h1[0]+h1[1])/h1[1] * knock2newV[1][j] - h1[0]/h1[1] * knock2newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];
								 
							}
							else {
								knock2newV[i][0] = (h1[0]+h1[1])/h1[1] * knock2newV[i][1] - h1[0]/h1[1] * knock2newV[i][2];
								 
							}
						} // knocknewV 완결



						for(i=0; i<=startSidx1; i++)
							for(j=0; j<=SM2; j++)  									
								knock2newV[i][j] = knocknewV[i][j];
									
						if(m_spreadXRate[curtN-1] > 0) {// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
							double spreadvalue;
							double   spmax;
		
								
							//for(i=spreadXSidx1; i<=SM1; i++) {
								//for(j=spreadXSidx2; j<=SM2; j++) {
							for(i=0; i<=SM1; i++) {
								for(j=0; j<=SM2; j++) {							
									if(i<=xidx1 || j<=xidx2) {

						if(S1[i]<=S2[j]) {							 
							spmax = S1[xidx1];																				 
						}
						else {						 
							spmax = S2[xidx2];																			 
						}						 
						spreadvalue = AMOUNT * ( m_spreadXRate[curtN-1] * (MIN(S1[i],S2[j])-spmax) + 1+m_Cpn[curtN-1]);
										knock2newV[i][j] = MAX(knock2newV[i][j], spreadvalue);								
									}
								}
							}
							
							/*
							//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						double spvalueR;

						spvalueR =  AMOUNT * (1+m_Cpn[curtN-1]);
						for(i=xidx1; i<=SM1; i++) {
							maxindexR = 0;
							maxvalueR = knock2newV[i][0];
							for(j=1; j<=xidx2; j++) {
								if(maxvalueR<knock2newV[i][j]) {
									maxindexR = j;
									maxvalueR = knock2newV[i][j];
								}
							}
							if(maxindexR < xidx2) {
								tempindex = 2*xidx2 - spreadXSidx2;
								for(j=maxindexR+1; j<tempindex; j++) {
									spreadvalue = (maxvalueR-spvalueR)/(S2[maxindexR]-S2[tempindex]) * (S2[j]-S2[tempindex]) + spvalueR;
									knock2newV[i][j] = MAX(knock2newV[i][j],spreadvalue);
								}
							}
						}

						for(j=xidx2; j<=SM2; j++) {														
							maxindexR = 0;
							maxvalueR = knock2newV[0][j];
							for(i=1; i<=xidx1; i++) {
								if(maxvalueR<knock2newV[i][j]) {
									maxindexR = i;
									maxvalueR = knock2newV[i][j];
								}
							}
							if(maxindexR < xidx1) {
								tempindex = 2*xidx1 - spreadXSidx1;
								for(i=maxindexR+1; i<tempindex; i++) {
									spreadvalue = (maxvalueR-spvalueR)/(S1[maxindexR]-S1[tempindex]) * (S1[i]-S1[tempindex]) + spvalueR;
									knock2newV[i][j] = MAX(knock2newV[i][j],spreadvalue);
								}
							}
						}
						*/
							

						
						
						} // X spread
 					 
					} 
	//////////////////////////////// KI2 end




					// no KI 영역 풀기 
					if(m_ki1Flag == 0 && m_ki2Flag == 0) {		
						for(i=startSidx1; i<=SM2; i++)
							tU[i][startSidx2] = knock2newV[i][startSidx2];

						for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기
							if(j<endSidx2) {
								for(i=startSidx1+1; i<SM1; i++) 
									vecV1[i] = oldV[i][j] + coeV1[i][j]*(-oldV[i+1][j-1]+oldV[i-1][j-1] + oldV[i+1][j+1] - oldV[i-1][j+1]); //herexxx
							
								tU[startSidx1][j] = knock1newV[startSidx1][j];
								for(i=startSidx1+1; i<SM1; i++)
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
								tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
								for(i=SM1-2; i>=startSidx1; i--) 				
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
										
								if(SmeshType1 == 0)				
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
								else
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				
							}
							else {													
								for(i=startSidx1+1; i<endSidx1; i++) // 행사 안되는 영역 
									vecV1[i] = oldV[i][j] + coeV1[i][j]*(-oldV[i+1][j-1]+oldV[i-1][j-1] + oldV[i+1][j+1] - oldV[i-1][j+1]); //herexxx
							
								tU[startSidx1][j] = knock1newV[startSidx1][j];
								for(i=startSidx1+1; i<endSidx1; i++) 
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
								for(i=endSidx1; i<=SM1; i++)
									tU[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(i=endSidx1-1; i>=startSidx1; i--)
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
							}
						}// for(j=1

						for(i=startSidx1; i<=SM1; i++) {				
							if(SmeshType2 == 0) 					
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							else					
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
						} // 일차 tU 완결 
 
						//////////////////////////////////////////////////////// newV 구하기!!
						for(j=startSidx2; j<=SM2; j++)
							newV[startSidx1][j]=knock1newV[startSidx1][j]; // D.C 적용
		 
						for(i=startSidx1+1; i<SM1; i++) {  //y축 고정시켜서.. 				 
							if(i<endSidx1) {
								for(j=startSidx2+1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);  //herexxx					   				 		 		

								newV[i][startSidx2] = knock2newV[i][startSidx2];
								for(j=startSidx2+1; j<SM2; j++) 
									newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];					
								newV[i][SM2-1] = newV[i][SM2-1]/YMAd[i][SM2-1];
								for(j=SM2-2; j>=startSidx2; j--) 
									newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];					
				
								if(SmeshType2 == 0)				
									newV[i][SM2] = 2*newV[i][SM2-1] - newV[i][SM2-2];
								else
									newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * newV[i][SM2-2];
							}
							else {													
								for(j=startSidx2+1; j<endSidx2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);  //herexxx					   				 		 		

								newV[i][startSidx2] = knock2newV[i][startSidx2];
								for(j=startSidx2+1; j<endSidx2; j++) 
									newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];					
								for(j=endSidx2; j<=SM2; j++)
									newV[i][j] = AMOUNT * (1+m_Cpn[curtN-1]);
								for(j=endSidx2-1; j>=startSidx2; j--) 
									newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];					
							}				
						} // for(j y축.. 고정 다 품

						for(j=startSidx2; j<=SM2; j++) {
							if(SmeshType1 == 0) 
								newV[SM1][j] = 2*newV[SM1-1][j] - newV[SM1-2][j];
							else
								newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * newV[SM1-2][j];
						} // newV 완결

						
						for(i=0; i<=startSidx1; i++)
							for(j=0; j<=startSidx2; j++) 								
									newV[i][j] = knocknewV[i][j];  // KI12

						for(i=0; i<=startSidx1; i++)
							for(j=startSidx2+1; j<=SM2; j++) 								
									newV[i][j] = knock1newV[i][j]; // KI1

						for(i=startSidx1+1; i<=SM1; i++)
							for(j=0; j<=startSidx2; j++) 								
									newV[i][j] = knock2newV[i][j]; // KI2
 

		
						if(m_spreadXRate[curtN-1] > 0) {// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
							double spreadvalue;
							double   spmax;
		
							
							//for(i=spreadXSidx1; i<=SM1; i++) {
								//for(j=spreadXSidx2; j<=SM2; j++) {
							for(i=0; i<=SM1; i++) {
								for(j=0; j<=SM2; j++) {
									if(i<=xidx1 || j<=xidx2) {
										 
							

						if(S1[i]<=S2[j]) {							 
							spmax = S1[xidx1];																				 
						}
						else {						 
							spmax = S2[xidx2];																			 
						}						 
						spreadvalue = AMOUNT * ( m_spreadXRate[curtN-1] * (MIN(S1[i],S2[j])-spmax) + 1+m_Cpn[curtN-1]);
										newV[i][j] = MAX(newV[i][j], spreadvalue);								
									}
								}
							}
						
							 /*
													//kzR here
							int maxindexR;
							int tempindex;
							double maxvalueR;
							double spvalueR;

							spvalueR =  AMOUNT * (1+m_Cpn[curtN-1]);
							for(i=xidx1; i<=SM1; i++) {
								maxindexR = 0;
								maxvalueR = newV[i][0];
								for(j=1; j<=xidx2; j++) {
									if(maxvalueR<newV[i][j]) {
										maxindexR = j;
										maxvalueR = newV[i][j];
									}
								}
								if(maxindexR < xidx2) {
									tempindex = 2*xidx2 - spreadXSidx2;
									for(j=maxindexR+1; j<tempindex; j++) {
										spreadvalue = (maxvalueR-spvalueR)/(S2[maxindexR]-S2[tempindex]) * (S2[j]-S2[tempindex]) + spvalueR;
										newV[i][j] = MAX(newV[i][j],spreadvalue);
									}
								}
							}

							for(j=xidx2; j<=SM2; j++) {														
								maxindexR = 0;
								maxvalueR = newV[0][j];
								for(i=1; i<=xidx1; i++) {
									if(maxvalueR<newV[i][j]) {
										maxindexR = i;
										maxvalueR = newV[i][j];
									}
								}
								if(maxindexR < xidx1) {
									tempindex = 2*xidx1 - spreadXSidx1;
									for(i=maxindexR+1; i<tempindex; i++) {
										spreadvalue = (maxvalueR-spvalueR)/(S1[maxindexR]-S1[tempindex]) * (S1[i]-S1[tempindex]) + spvalueR;
										newV[i][j] = MAX(newV[i][j],spreadvalue);
									}
								}
							}
							*/
							 
							
						} // X spread
 				 
 					 
					}// if(m_kiFlag != 1)

 
				} //if (IsMidCompu == 0) else				 

			}//for(intraLoop)

	   
			//이산배당 적용
			disc_div1 = 0.0;
			for(i=0; i<m_divN1; i++) {
				if(matuLoop == (m_divBDay1[i]-1)) {
					disc_div1 = m_divCash1[i];
					break;
				}
			}
		
			disc_div2 = 0.0;
			for(j=0; j<m_divN2; j++) {
				if(matuLoop == (m_divBDay2[j]-1)) {
					disc_div2 = m_divCash2[j];
					break;
				}
			}

			if(disc_div1 > 0.0 || disc_div2 > 0.0) {
				double divSlevel1[2], divSlevel2[2];
				double pole[]={0, 1};
				double alpha1, alpha2;						

				//memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); 
				for(i=0;i<=SM1; i++)
					memcpy(knockoldV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
				/*
				for(i=0; i<=SM1; i++)
					for(j=0; j<=SM2; j++)
						knockoldV[i][j] = knocknewV[i][j];
				*/

				divSlevel1[0] = (m_divbval1*m_divApy - disc_div1)/m_bval1;
				divSlevel1[1] = (m_divbval1*m_divApy + disc_div1)/m_bval1;

				divSlevel2[0] = (m_divbval2*m_divApy - disc_div2)/m_bval2;
				divSlevel2[1] = (m_divbval2*m_divApy + disc_div2)/m_bval2;

				for(i=0; i<=SM1; i++) {
					if(disc_div1 > 0.0) 				 
						alpha1 = interp1(S1[i],divSlevel1,pole,2,1,0);
					else
						alpha1 = 0.0;

					for(j=0; j<=SM2; j++) {
						if(disc_div2 > 0.0)					
							alpha2 = interp1(S2[j],divSlevel2,pole,2,1,0);
						else
							alpha2 = 0.0;
					
						knocknewV[i][j] = interp2(S1[i]-alpha1*disc_div1/m_bval1,S2[j]-alpha2*disc_div2/m_bval2,S1,S2,knockoldV,SM1,SM2); // SM1은.. S1 분할 수.. 0부터 시작.. 마지막 인덱스
					}
				}		

				
				if(m_ki2Flag == 0) {								 
					for(i=0; i<=SM1; i++)
						memcpy(knock1oldV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
				 

					divSlevel1[0] = (m_divbval1*m_divApy - disc_div1)/m_bval1;
					divSlevel1[1] = (m_divbval1*m_divApy + disc_div1)/m_bval1;

					divSlevel2[0] = (m_divbval2*m_divApy - disc_div2)/m_bval2;
					divSlevel2[1] = (m_divbval2*m_divApy + disc_div2)/m_bval2;

					for(i=0; i<=SM1; i++) {
						if(disc_div1 > 0.0) 				 
							alpha1 = interp1(S1[i],divSlevel1,pole,2,1,0);
						else
							alpha1 = 0.0;

						for(j=0; j<=SM2; j++) {
							if(disc_div2 > 0.0)					
								alpha2 = interp1(S2[j],divSlevel2,pole,2,1,0);
							else
								alpha2 = 0.0;
					
							knock1newV[i][j] = interp2(S1[i]-alpha1*disc_div1/m_bval1,S2[j]-alpha2*disc_div2/m_bval2,S1,S2,knock1oldV,SM1,SM2); // SM1은.. S1 분할 수.. 0부터 시작.. 마지막 인덱스
						}
					}			
				}  
							
				if(m_ki1Flag == 0) {								 
					for(i=0; i<=SM1; i++)
						memcpy(knock2oldV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
				 

					divSlevel1[0] = (m_divbval1*m_divApy - disc_div1)/m_bval1;
					divSlevel1[1] = (m_divbval1*m_divApy + disc_div1)/m_bval1;

					divSlevel2[0] = (m_divbval2*m_divApy - disc_div2)/m_bval2;
					divSlevel2[1] = (m_divbval2*m_divApy + disc_div2)/m_bval2;

					for(i=0; i<=SM1; i++) {
						if(disc_div1 > 0.0) 				 
							alpha1 = interp1(S1[i],divSlevel1,pole,2,1,0);
						else
							alpha1 = 0.0;

						for(j=0; j<=SM2; j++) {
							if(disc_div2 > 0.0)					
								alpha2 = interp1(S2[j],divSlevel2,pole,2,1,0);
							else
								alpha2 = 0.0;
					
							knock2newV[i][j] = interp2(S1[i]-alpha1*disc_div1/m_bval1,S2[j]-alpha2*disc_div2/m_bval2,S1,S2,knock2oldV,SM1,SM2); // SM1은.. S1 분할 수.. 0부터 시작.. 마지막 인덱스
						}
					}			
				}  


			
				if(m_ki1Flag == 0 && m_ki2Flag == 0) {								 
					for(i=0; i<=SM1; i++)
						memcpy(oldV[i],newV[i],(size_t)((SM2+1)*sizeof(double)));
					/*
					for(i=0; i<=SM1; i++)
						for(j=0; j<=SM2; j++)
							oldV[i][j] = newV[i][j];
					*/

					divSlevel1[0] = (m_divbval1*m_divApy - disc_div1)/m_bval1;
					divSlevel1[1] = (m_divbval1*m_divApy + disc_div1)/m_bval1;

					divSlevel2[0] = (m_divbval2*m_divApy - disc_div2)/m_bval2;
					divSlevel2[1] = (m_divbval2*m_divApy + disc_div2)/m_bval2;

					for(i=0; i<=SM1; i++) {
						if(disc_div1 > 0.0) 				 
							alpha1 = interp1(S1[i],divSlevel1,pole,2,1,0);
						else
							alpha1 = 0.0;

						for(j=0; j<=SM2; j++) {
							if(disc_div2 > 0.0)					
								alpha2 = interp1(S2[j],divSlevel2,pole,2,1,0);
							else
								alpha2 = 0.0;
					
							newV[i][j] = interp2(S1[i]-alpha1*disc_div1/m_bval1,S2[j]-alpha2*disc_div2/m_bval2,S1,S2,oldV,SM1,SM2); // SM1은.. S1 분할 수.. 0부터 시작.. 마지막 인덱스
						}
					}			
				} //if(m_knockFlag !=1 )

			} //if(disc_div > 0.0)
			//이산배당 적용 end
 			
		 
			
			// KI 종류에 따라..
			if(m_ki1Flag == 1 && m_ki2Flag == 1) { 
				for(i=0; i<=SM1; i++)
					memcpy(newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
				for(i=0; i<=SM1; i++)
					memcpy(knock1newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
				for(i=0; i<=SM1; i++)
					memcpy(knock2newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
			}
		
			if(m_ki1Flag == 1 && m_ki2Flag == 0) { 
				for(i=0; i<=SM1; i++)
					memcpy(newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
				for(i=0; i<=SM1; i++)
					memcpy(knock2newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
			}
		
			if(m_ki1Flag == 0 && m_ki2Flag == 1) { 
				for(i=0; i<=SM1; i++)
					memcpy(newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
				for(i=0; i<=SM1; i++)
					memcpy(knock1newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
			}
	

		 
			if(m_curtgreek == 1 && matuLoop == 1) { // 기본 그릭	Theta를 위한 +1일 저장
				for(i=0; i<outN + 4; i++) {
					for(j=0; j<outN + 4; j++) {
						sim_knock1Flag = 0;					
						sim_knock2Flag = 0;					
						if(SrTemp1[i] <= m_kihval1)						
							sim_knock1Flag = 1;					
						if(SrTemp2[j] <= m_kihval2)						
							sim_knock2Flag = 1;

						if(evalq<=0) {
							sleveltemp1 = SrTemp1[i]/m_bval1;
							sleveltemp2 = SrTemp2[j]/m_bval2;
						}
						else { // n영업일 평균 처리
							if(m_ctime>15) { // 장종료
								sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN[curtN-1]-evalq); //  = m_matu
								for(evali=0; evali<evalq; evali++)
									sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
								sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
														
								sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN[curtN-1]-evalq); //  = m_matu
								for(evali=0; evali<evalq; evali++)
									sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
								sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];
							}
							else { // 장중
								sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN[curtN-1]-evalq+1);
								for(evali=1; evali<evalq; evali++)
									sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
								sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
																					
								sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN[curtN-1]-evalq+1);
								for(evali=1; evali<evalq; evali++)
									sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
								sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];
							}
						}				
								
						if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
							SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knocknewV, SM1, SM2);  
						if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
							SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock1newV, SM1, SM2);
						if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
							SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock2newV, SM1, SM2);
						if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
							SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2);
 
					}// for(j			
				} //for(i
			} //(m_curtgreek == 1 && matuLoop == 1)
 	   
 	   

			if(m_curtgreek == 1 && mid_n > 0 && matuLoop == 0) { // 오늘이 조기 상환일 중에 하나
				if(m_ctime > 15 && slevel1 >= xlevel1 && slevel1 >= xlevel2) { // 오늘이 조기상환일 이면서 종가 이후 상환 된 때..
					for(i=0; i<outN + 4; i++) {	
						for(j=0; j<outN + 4; j++) {													
							sim_knock1Flag = 0;											
							sim_knock2Flag = 0;											
							if(SrTemp1[i] <= m_kihval1)			
								sim_knock1Flag = 1;									
							if(SrTemp2[j] <= m_kihval2)										
								sim_knock2Flag = 1;

							sleveltemp1 = SrTemp1[i]/m_bval1;
							for(evali=1; evali<evalq; evali++)
								sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
							sleveltemp1 = sleveltemp1 / m_evalN[curtN-1];
												
							sleveltemp2 = SrTemp2[j]/m_bval2;
							for(evali=1; evali<evalq; evali++)
								sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
							sleveltemp2 = sleveltemp2 / m_evalN[curtN-1];

							
							if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
								SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knocknewV, SM1, SM2);  
							if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
								SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock1newV, SM1, SM2);
							if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
								SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock2newV, SM1, SM2);
							if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
								SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2);
 
						}	
					}
				}
				if(m_ctime > 15 && (slevel1 < xlevel1 || slevel2 < xlevel2)) { //조기상환일에 조기 상환 안된 경우
					expiryq = 0;
					slevel1 = m_sval1/m_bval1;
					slevel2 = m_sval2/m_bval2;
				
					evalq = m_evalN[curtN] - m_matu[curtN]; // 조기 상환일에 조기상환 안된 경우 evalq  다음번으로 
				}
			}
	 
				
		}//for(matuLoop)	

		if(m_matu[m_matuN-1]>0) {
			free(dailyIFrf);
			free(dailyIFdc);

		}


		///////////////////////////////////////////--------------Debug	Mode s 

	 
		int debugFlag;
		debugFlag = 1;
		if(debugFlag) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/debug_sdbgnufdprice.txt",ios::out);	 
	 
			logsgnuplot << "S1= [" ;
			for(i=0; i<=SM1-1; i++) {
				logsgnuplot << S1[i]*m_bval1 << endl;
			}
			logsgnuplot << S1[SM1]*m_bval1 << "];" << endl;

		
			logsgnuplot << "S2= [";
			for(i=0; i<=SM2-1; i++) {
				logsgnuplot << S2[i]*m_bval2 << endl;
			}
			logsgnuplot << S2[SM2]*m_bval2 << "];" << endl;

			logsgnuplot << "newC=[" ;		
			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++)
					if(i == SM1 && j == SM2 )				
						logsgnuplot << newV[i][j] << "];" ;
					else				
						logsgnuplot << newV[i][j] << "  " ;
				logsgnuplot << endl;
			}
			logsgnuplot.close();


			logsgnuplot.open("c:/logland/debug_avesdbgnufdprice.txt",ios::out);	 		 
			logsgnuplot << "S1= [" ;
			for(i=0; i<=SM1-1; i++) {
				logsgnuplot << S1[i]*m_bval1 << endl;
			}
			logsgnuplot << S1[SM1]*m_bval1 << "];" << endl;
		
		
			logsgnuplot << "S2= [";
			for(i=0; i<=SM2-1; i++) {
				logsgnuplot << S2[i]*m_bval2 << endl;
			}
			logsgnuplot << S2[SM2]*m_bval2 << "];" << endl;

			logsgnuplot << "avenewC=[" ;

			double stemp1,stemp2, valuetemp;
			for(i=0; i<=SM1; i++) {						
				for(j=0; j<=SM2; j++) {
											
					sim_knock1Flag = 0;										
					sim_knock2Flag = 0;						 

					if(S1[i]*m_bval1 <= m_kihval1) 
						sim_knock1Flag = 1;
					if(S2[j]*m_bval2 <= m_kihval2)
						sim_knock2Flag = 1;

					if(evalq <=1) {
						stemp1 = S1[i]*m_bval1;
						stemp2 = S2[j]*m_bval2;
					}
					else {
						stemp1 = S1[i]*m_bval1 * (m_evalN[curtN-1]-(evalq-1));			
						for(evali=1;evali<evalq; evali++)				
							stemp1 = stemp1 + m_psval1[evali];
						stemp1 = stemp1 / m_evalN[curtN-1];

						stemp2 = S2[j]*m_bval2 * (m_evalN[curtN-1]-(evalq-1));			
						for(evali=1;evali<evalq; evali++)				
							stemp2 = stemp2 + m_psval2[evali];
						stemp2 = stemp2 / m_evalN[curtN-1];
					}

					
					if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
						valuetemp = interp2(stemp1/m_bval1, stemp2/m_bval2, S1, S2, knocknewV, SM1, SM2);  
					if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
						valuetemp = interp2(stemp1/m_bval1, stemp2/m_bval2, S1, S2, knock1newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
						valuetemp = interp2(stemp1/m_bval1, stemp2/m_bval2, S1, S2, knock2newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
						valuetemp = interp2(stemp1/m_bval1, stemp2/m_bval2, S1, S2, newV, SM1, SM2);

					 												
					if(i == SM1 && j == SM2 )				
						logsgnuplot << valuetemp << "];" ;
					else				
						logsgnuplot << valuetemp << "  " ;				 
				}
				logsgnuplot << endl;
			}
			logsgnuplot.close();
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	 
			
	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/kzhere3.txt",ios::out);	 
	 
			 
			logsgnuplot << "m_ki1Flag = " << m_ki1Flag << endl;
			
			logsgnuplot << "m_ki2Flag = " << m_ki2Flag << endl;

			logsgnuplot << "newV[49][100] = " << newV[49][100] << endl;
			logsgnuplot << "knocknewV[49][100] = " << knocknewV[49][100] << endl;


		 
			logsgnuplot.close();

 
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	 
		if(m_curtgreek == 1 || m_curtgreek == 2) {
			for(i=0; i<outN + 4; i++) {			
				for(j=0; j<outN + 4; j++) {
											
					sim_knock1Flag = 0;										
					sim_knock2Flag = 0;										
					if(SrTemp1[i] <= m_kihval1)								
						sim_knock1Flag = 1;										
					if(SrTemp2[j] <= m_kihval2)								
						sim_knock2Flag = 1;

					if(evalq <= 1) {
						sleveltemp1 = SrTemp1[i]/m_bval1;
						sleveltemp2 = SrTemp2[j]/m_bval2;
					}
					else {
						sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN[curtN-1] - (evalq-1));
						for(evali=1; evali<evalq; evali++)
							sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
						sleveltemp1 = sleveltemp1/m_evalN[curtN-1];
										
						sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN[curtN-1] - (evalq-1));
						for(evali=1; evali<evalq; evali++)
							sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
						sleveltemp2 = sleveltemp2/m_evalN[curtN-1];
					}
			
										
					if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
						SrVTemp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knocknewV, SM1, SM2);  
					if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
						SrVTemp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock1newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
						SrVTemp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock2newV, SM1, SM2);
					if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
						SrVTemp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2);

					 
				}
			}
		
			Baseidx = SoutRange + 2; // m_sval에 해당하는 idx 위 아래 2개씩 더 계산해서.. +2 필요
		}

		if(m_curtgreek == 1) { // 기본 그릭  
			value = interp2(slevel1,slevel2,S1,S2,newV,SM1,SM2);
  
	/*		
			greeks[twoPriceidx] = SrVTemp[Baseidx][Baseidx]; // [0] = price

			greeks[twoUpDeltaCidx1] = (SrVTemp[Baseidx+1][Baseidx] - SrVTemp[Baseidx][Baseidx]) / 0.01; // [1] = up delta cash = delta * s
			greeks[twoDownDeltaCidx1] = (SrVTemp[Baseidx][Baseidx] - SrVTemp[Baseidx-1][Baseidx]) / 0.01; // [2] = down delta cash		
			greeks[twoUpDeltaCidx2] = (SrVTemp[Baseidx][Baseidx+1] - SrVTemp[Baseidx][Baseidx]) / 0.01; // [3] = up delta cash = delta * s
			greeks[twoDownDeltaCidx2] = (SrVTemp[Baseidx][Baseidx] - SrVTemp[Baseidx][Baseidx-1]) / 0.01; // [4] = down delta cash

			greeks[twoUpGammaCidx1] = (SrVTemp[Baseidx+2][Baseidx] - 2*SrVTemp[Baseidx+1][Baseidx] + SrVTemp[Baseidx][Baseidx]) / 0.01; // [5] = up gamma cash = gamma * s^2 * 0.01
			greeks[twoDownGammaCidx1] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx-1][Baseidx] + SrVTemp[Baseidx-2][Baseidx]) / 0.01; // [6] = down gamma cash		
			greeks[twoUpGammaCidx2] = (SrVTemp[Baseidx][Baseidx+2] - 2*SrVTemp[Baseidx][Baseidx+1] + SrVTemp[Baseidx][Baseidx]) / 0.01; // [7] = up gamma cash = gamma * s^2 * 0.01
			greeks[twoDownGammaCidx2] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx][Baseidx-1] + SrVTemp[Baseidx][Baseidx-2]) / 0.01; // [8] = down gamma cash

			greeks[twoCrossGammaCidx] = (SrVTemp[Baseidx+1][Baseidx+1] - SrVTemp[Baseidx+1][Baseidx-1] - SrVTemp[Baseidx-1][Baseidx+1] + SrVTemp[Baseidx-1][Baseidx-1]) / (4* 0.01); // [9] = cross gamma cash

			greeks[twoThetaidx] = (SrV1Temp[Baseidx][Baseidx] - SrVTemp[Baseidx][Baseidx]); // [10] = 1 day theta
			greeks[twoCharmCidx1] = (SrV1Temp[Baseidx+1][Baseidx]-SrV1Temp[Baseidx-1][Baseidx]-SrVTemp[Baseidx+1][Baseidx]+SrVTemp[Baseidx-1][Baseidx])/(2*0.01); // [6] = 1 day charm cash
			greeks[twoCharmCidx2] = (SrV1Temp[Baseidx][Baseidx+1]-SrV1Temp[Baseidx][Baseidx-1]-SrVTemp[Baseidx][Baseidx+1]+SrVTemp[Baseidx][Baseidx-1])/(2*0.01); // [6] = 1 day charm cash
	*/	
		
			greeks[twoPriceidx] = SrVTemp[Baseidx][Baseidx]; // [0] = price
		
			//2012-05-23 1% delta cash로 수정
			greeks[twoUpDeltaCidx1] = (SrVTemp[Baseidx+1][Baseidx] - SrVTemp[Baseidx][Baseidx]); // [1] = up delta cash = delta * 0.01s
			greeks[twoDownDeltaCidx1] = (SrVTemp[Baseidx][Baseidx] - SrVTemp[Baseidx-1][Baseidx]); // [2] = down delta cash		
			greeks[twoUpDeltaCidx2] = (SrVTemp[Baseidx][Baseidx+1] - SrVTemp[Baseidx][Baseidx]); // [3] = up delta cash = delta * 0.01s
			greeks[twoDownDeltaCidx2] = (SrVTemp[Baseidx][Baseidx] - SrVTemp[Baseidx][Baseidx-1]); // [4] = down delta cash
						
			//델타 완화 2012-07-12
			greeks[twoUpDeltaCidx1] = MAX(MIN(greeks[twoUpDeltaCidx1],maxhgdelta),minhgdelta);  // deltaadd
			greeks[twoDownDeltaCidx1] = MAX(MIN(greeks[twoDownDeltaCidx1],maxhgdelta),minhgdelta);  // deltaadd;
			greeks[twoUpDeltaCidx2] = MAX(MIN(greeks[twoUpDeltaCidx2],maxhgdelta),minhgdelta);  // deltaadd;
			greeks[twoDownDeltaCidx2] = MAX(MIN(greeks[twoDownDeltaCidx2],maxhgdelta),minhgdelta);  // deltaadd;
		
			//2012-05-23 1% gamma cash로 수정
			greeks[twoUpGammaCidx1] = (SrVTemp[Baseidx+2][Baseidx] - 2*SrVTemp[Baseidx+1][Baseidx] + SrVTemp[Baseidx][Baseidx]); // [5] = up gamma cash = gamma * (0.01s)^2 
			greeks[twoDownGammaCidx1] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx-1][Baseidx] + SrVTemp[Baseidx-2][Baseidx]); // [6] = down gamma cash		
			greeks[twoUpGammaCidx2] = (SrVTemp[Baseidx][Baseidx+2] - 2*SrVTemp[Baseidx][Baseidx+1] + SrVTemp[Baseidx][Baseidx]); // [7] = up gamma cash = gamma * (0.01s)^2 
			greeks[twoDownGammaCidx2] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx][Baseidx-1] + SrVTemp[Baseidx][Baseidx-2]); // [8] = down gamma cash
		
			//2012-05-23 1% cross gamma cash로 수정
			greeks[twoCrossGammaCidx] = (SrVTemp[Baseidx+1][Baseidx+1] - SrVTemp[Baseidx+1][Baseidx-1] - SrVTemp[Baseidx-1][Baseidx+1] + SrVTemp[Baseidx-1][Baseidx-1]) /4; // [9] = 1% cross gamma cash 

			greeks[twoThetaidx] = (SrV1Temp[Baseidx][Baseidx] - SrVTemp[Baseidx][Baseidx]); // [10] = 1 day theta

			//2012-05-23 1% charm cash로 수정
			greeks[twoCharmCidx1] = (SrV1Temp[Baseidx+1][Baseidx]-SrV1Temp[Baseidx-1][Baseidx]-SrVTemp[Baseidx+1][Baseidx]+SrVTemp[Baseidx-1][Baseidx])/2; // [6] = 1 day 1% charm cash1
			greeks[twoCharmCidx2] = (SrV1Temp[Baseidx][Baseidx+1]-SrV1Temp[Baseidx][Baseidx-1]-SrVTemp[Baseidx][Baseidx+1]+SrVTemp[Baseidx][Baseidx-1])/2; // [6] = 1 day 1% charm cash2


 
			if(m_batchFlag) {

				using namespace std;
				using std::ofstream;

				ofstream batchsf;
			
				char *fname;
				fname = new char[kCodeN+1];

				strcpy_s(fname,kCodeN,m_kCode);
				strcat_s(fname,kCodeN,".txt");
		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 			
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp1[i] << " ";
				}
				batchsf << endl;

				
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp2[i] << " ";
				}
				batchsf << endl;
			
				for(i=0; i<outN + 4; i++)  {
					for(j=0; j<outN + 4; j++) {
						batchsf << SrVTemp[i][j] << " ";
					}
					batchsf << endl;
				}
				batchsf.close();


				strcpy_s(fname,kCodeN,m_kCode);
				strcat_s(fname,kCodeN,"_day+.txt");

				batchsf.open(fname,ios::out);
			
				batchsf << fixed;
				batchsf.precision(15);
				 
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp1[i] << " ";
				}
				batchsf << endl;
			
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp2[i] << " ";
				}
				batchsf << endl;

				for(i=0; i<outN + 4; i++)  {
					for(j=0; j<outN + 4; j++) {
						batchsf << SrV1Temp[i][j] << " ";
					}
					batchsf << endl;
				}
				batchsf.close();
	/////
		
				if(m_matu[curtN-1] == 0 && m_ctime < 15) {		
					char *tempstring;

					tempstring = new char[kCodeN+1];
					intraTN = 24;//m_meshTN1;						 
					dt = (1.0/YFactor)/intraTN;

					r = m_iRFrate[0]; 
					rd = m_iDCrate[0]; 

					if(m_voltype == 0) {
						volskew1[0] = m_vol1[0][0];
						volskew2[0] = m_vol2[0][0];
					}
					else {
						curttimeVol1(volskew1, 0);
						curttimeVol2(volskew2, 0);
					}

					cor = curttimeCorr(0);
		 
					// X 생성		 
					for(j=1; j<SM2; j++) {  //y축 고정시켜서.. 
						if(m_voltype == 2) {				
							Srv2 = S2[j]*m_bval2 / m_volbval2;				
							v2 = interp1(Srv2, m_volSmness,volskew2,m_volSN,1,0);
						}
						else {
							v2 = volskew2[0];
						}

						for(i=1; i<SM1; i++) { //x축 풀기 위하여 행렬 만들기
							if(m_voltype == 2) {
								Srv1 = S1[i]*m_bval1 / m_volbval1;
								v1 = interp1(Srv1, m_volSmness, volskew1, m_volSN,1,0);
							}
							else {
								v1 = volskew1[0];
							}

							if(SmeshType1 == 0) { //균등 mesh
								knockXMAl[i][j] = -0.5*(v1*v1)*(i*i)*dt + (r-q1)*i*dt/2.0 ;// + 0.5*cor*v1*v2*i*j*dt/2.0;
								knockXMAd[i][j] = v1*v1*i*i*dt + 0.5*rd*dt + 1;
								knockXMAu[i][j] = -0.5*(v1*v1)*(i*i)*dt - (r-q1)*i*dt/2.0 ;// - 0.5*cor*v1*v2*i*j*dt/2.0;

								coeV1[i][j] = 0.5*cor*v1*v2*i*j*dt/4.0;   //0.5*cor*v1*v2*i*j*dt/2.0;
							}
							else {						
								knockXMAl[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*(h1[i]+h1[i-1])) + (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]); // + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
								knockXMAd[i][j] = (v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*h1[i]) + 0.5*rd*dt + 1;
								knockXMAu[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i]*(h1[i]+h1[i-1])) - (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]); // - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);								
					
								coeV1[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*(h2[j]+h2[j-1]));	//0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
							}				
			
							if(m_ki1Flag == 0) {							
								XMAl[i][j] = knockXMAl[i][j];				
								XMAd[i][j] = knockXMAd[i][j];					
								XMAu[i][j] = knockXMAu[i][j];
							}

						} // 행렬 만들기					
						knockXMAl[0][j] = 0.0;				
						knockXMAd[0][j] = 1.0;				
						knockXMAu[0][j] = 0.0;		

						//경계 선형 조건
						if(SmeshType1 == 0) {
							knockXMAl[SM1-1][j] = knockXMAl[SM1-1][j] - knockXMAu[SM1-1][j];
							knockXMAd[SM1-1][j] = knockXMAd[SM1-1][j] + 2*knockXMAu[SM1-1][j];
						}
						else {			 
							knockXMAl[SM1-1][j] = knockXMAl[SM1-1][j] - knockXMAu[SM1-1][j]*h1[SM1-1]/h1[SM1-2];
							knockXMAd[SM1-1][j] = knockXMAd[SM1-1][j] + knockXMAu[SM1-1][j]*(h1[SM1-1]+h1[SM1-2])/h1[SM1-2];			
						}

						if(m_ki1Flag == 0) {
							XMAl[SM1-1][j] = knockXMAl[SM1-1][j];			
							XMAd[SM1-1][j] = knockXMAd[SM1-1][j];

							XMAl[startSidx1][j] = 0;		
							XMAd[startSidx1][j] = 1;			
							XMAu[startSidx1][j] = 0;	
						
							//LU decomposition			
							for(i=startSidx1+1; i<=SM1-1; i++) {				
								XMAl[i][j] = XMAl[i][j]/XMAd[i-1][j];				
								XMAd[i][j] = XMAd[i][j] - XMAl[i][j]*XMAu[i-1][j];			
							}			
						}

						//LU decomposition
						for(i=1; i<=SM1-1; i++) {
							knockXMAl[i][j] = knockXMAl[i][j]/knockXMAd[i-1][j];
							knockXMAd[i][j] = knockXMAd[i][j] - knockXMAl[i][j]*knockXMAu[i-1][j];
						}				
					} //X 행렬 끝
	 
					// Y 생성		 
					for(i=1; i<SM1; i++) {  //y축 고정시켜서.. 
						if(m_voltype == 2) {				
							Srv1 = S1[i]*m_bval1 / m_volbval1;				
							v1 = interp1(Srv1, m_volSmness,volskew1,m_volSN,1,0);
						}
						else {
							v1 = volskew1[0];
						}

						for(j=1; j<SM2; j++) { //x축 풀기 위하여 행렬 만들기
							if(m_voltype == 2) {
								Srv2 = S2[j]*m_bval2 / m_volbval2;
								v2 = interp1(Srv2, m_volSmness, volskew2, m_volSN,1,0);
							}
							else {
								v2 = volskew2[0];
							}

							if(SmeshType2 == 0) { //균등 mesh
								knockYMAl[i][j] = -0.5*(v2*v2)*(j*j)*dt + (r-q2)*j*dt/2.0; // + 0.5*cor*v1*v2*i*j*dt/2.0;
								knockYMAd[i][j] = v2*v2*j*j*dt + 0.5*rd*dt + 1;
								knockYMAu[i][j] = -0.5*(v2*v2)*(j*j)*dt - (r-q2)*j*dt/2.0; // - 0.5*cor*v1*v2*i*j*dt/2.0;

								coeV2[i][j] = 0.5*cor*v1*v2*i*j*dt/4.0; //0.5*cor*v1*v2*i*j*dt/2.0;
							}
							else {						
								knockYMAl[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*(h2[j]+h2[j-1])) + (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]); // + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);
								knockYMAd[i][j] = (v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*h2[j]) + 0.5*rd*dt + 1;
								knockYMAu[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j]*(h2[j]+h2[j-1])) - (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]); // - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);							
					
								coeV2[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*(h1[i]+h1[i-1]));	//0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);
							}				
			
							if(m_ki2Flag == 0) {							
								YMAl[i][j] = knockYMAl[i][j];				
								YMAd[i][j] = knockYMAd[i][j];					
								YMAu[i][j] = knockYMAu[i][j];
							}

						} // 행렬 만들기					
						knockYMAl[i][0] = 0.0;				
						knockYMAd[i][0] = 1.0;				
						knockYMAu[i][0] = 0.0;	

						//경계 선형 조건
						if(SmeshType2 == 0) {
							knockYMAl[i][SM2-1] = knockYMAl[i][SM2-1] - knockYMAu[i][SM2-1];
							knockYMAd[i][SM2-1] = knockYMAd[i][SM2-1] + 2*knockYMAu[i][SM2-1];
						}
						else {			 
							knockYMAl[i][SM2-1] = knockYMAl[i][SM2-1] - knockYMAu[i][SM2-1]*h2[SM2-1]/h2[SM2-2];
							knockYMAd[i][SM2-1] = knockYMAd[i][SM2-1] + knockYMAu[i][SM2-1]*(h2[SM2-1]+h2[SM2-2])/h2[SM2-2];			
						}

						if(m_ki2Flag == 0) {
							YMAl[i][SM2-1] = knockYMAl[i][SM2-1];			
							YMAd[i][SM2-1] = knockYMAd[i][SM2-1];

							YMAl[i][startSidx2] = 0;		
							YMAd[i][startSidx2] = 1;			
							YMAu[i][startSidx2] = 0;	
						
							//LU decomposition			
							for(j=startSidx2+1; j<=SM2-1; j++) {				
								YMAl[i][j] = YMAl[i][j]/YMAd[i][j-1];				
								YMAd[i][j] = YMAd[i][j] - YMAl[i][j]*YMAu[i][j-1];			
							}			
						}

						//LU decomposition
						for(j=1; j<=SM2-1; j++) {
							knockYMAl[i][j] = knockYMAl[i][j]/knockYMAd[i][j-1];									
							knockYMAd[i][j] = knockYMAd[i][j] - knockYMAl[i][j]*knockYMAu[i][j-1];
						}				
					} //Y 행렬 끝		 
	   
 
  
		
					for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 확인							 
						DCValue12 = DCValue12 * 1/(1+rd/(YFactor*intraTN));
						DCValue1 = DCValue1 * 1/(1+rd/(YFactor*intraTN));
						DCValue2 = DCValue2 * 1/(1+rd/(YFactor*intraTN));

						
						//memcpy(knockoldV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
						for(i=0; i<=SM1; i++)			
							memcpy(knockoldV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
			
						if(m_ki2Flag == 0)
							for(i=0; i<=SM1; i++)				
								memcpy(knock1oldV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));

						if(m_ki1Flag == 0)
							for(i=0; i<=SM1; i++)				
								memcpy(knock2oldV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));

						if(m_ki1Flag == 0 && m_ki2Flag == 0)
							for(i=0; i<=SM1; i++)										
								memcpy(oldV[i],newV[i],(size_t)((SM2+1)*sizeof(double)));
				  
										
						for(i=0; i<=SM1; i++) 						
							tU[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;

						for(j=1; j<SM2; j++) { //y축 고정 풀기

							for(i=1; i<SM1; i++) 
								vecV1[i] = knockoldV[i][j] + coeV1[i][j]*(-knockoldV[i+1][j-1]+knockoldV[i-1][j-1] +knockoldV[i+1][j+1]-knockoldV[i-1][j+1]); // herexxx

							tU[0][j] = knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(i=1; i<SM1; i++)
								tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];
							tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
							for(i=SM1-2; i>=0; i--) 				
								tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
										
							if(SmeshType1 == 0)				
								tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
							else
								tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				

						}// for(j=1

						for(i=0; i<=SM1; i++) {				
							if(SmeshType2 == 0) 					
								tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
							else					
								tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
						} // 일차 tU 완결 

						//////////////////////////////////////////////////////// newV 구하기!!
						for(j=0; j<=SM2; j++)
							knocknewV[0][j]=knockoldV[0][j]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
		 
						for(i=1; i<SM1; i++) {  //x축 고정시켜서.. 				 

							for(j=1; j<SM2; j++) 				
								vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

							knocknewV[i][0] = knockoldV[i][0]*1/(1+rd/(YFactor*intraTN));//DCValue12; //0.0;
							for(j=1; j<SM2; j++) 
								knocknewV[i][j] = vecV2[j] - knockYMAl[i][j]*knocknewV[i][j-1];					
							knocknewV[i][SM2-1] = knocknewV[i][SM2-1]/knockYMAd[i][SM2-1];
							for(j=SM2-2; j>=0; j--) 
								knocknewV[i][j] = (knocknewV[i][j] - knockYMAu[i][j]*knocknewV[i][j+1])/knockYMAd[i][j];					
				
							if(SmeshType2 == 0)				
								knocknewV[i][SM2] = 2*knocknewV[i][SM2-1] - knocknewV[i][SM2-2];
							else
								knocknewV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knocknewV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knocknewV[i][SM2-2];
				
						} // for(i x축.. 고정 다 품

						for(j=0; j<=SM2; j++) {
							if(SmeshType1 == 0) 
								knocknewV[SM1][j] = 2*knocknewV[SM1-1][j] - knocknewV[SM1-2][j];
							else
								knocknewV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knocknewV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knocknewV[SM1-2][j];
						} // knocknewV 완결

												//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knocknewV[0][j] = 2*knocknewV[1][j] - knocknewV[2][j];
								 
							}
							else {
								knocknewV[0][j] = (h1[0]+h1[1])/h1[1] * knocknewV[1][j] - h1[0]/h1[1] * knocknewV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knocknewV[i][0] = 2*knocknewV[i][1] - knocknewV[i][2];
								 
							}
							else {
								knocknewV[i][0] = (h1[0]+h1[1])/h1[1] * knocknewV[i][1] - h1[0]/h1[1] * knocknewV[i][2];
								 
							}
						} // knocknewV 완결




		//////////////////////////// KI1				//herework!!					
						if(m_ki2Flag == 0) {						 
						 
							for(i=0; i<=SM1; i++)
								tU[i][startSidx2] = knocknewV[i][startSidx2];

							for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기						
								for(i=1; i<SM1; i++) 
									vecV1[i] = knock1oldV[i][j] + coeV1[i][j]*(-knock1oldV[i+1][j-1]+knock1oldV[i-1][j-1] +knock1oldV[i+1][j+1]-knock1oldV[i-1][j+1]); //herexxx
                      
								tU[0][j] = knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1;
								for(i=1; i<SM1; i++)
									tU[i][j] = vecV1[i] - knockXMAl[i][j]*tU[i-1][j];

								tU[SM1-1][j] = tU[SM1-1][j]/knockXMAd[SM1-1][j];
								for(i=SM1-2; i>=0; i--) 											
									tU[i][j] = (tU[i][j] - knockXMAu[i][j]*tU[i+1][j])/knockXMAd[i][j];
													
								if(SmeshType1 == 0) {		
									//tU[0][j] = 2*tU[1][j] - tU[2][j];
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
								}
								else {
									//tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
								}
						 			

							}// for(j=1

							for(i=0; i<=SM1; i++) {				
								if(SmeshType2 == 0) 			
									tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
								else					
									tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
							} // 일차 tU 완결  
						
								//호
							for(j=startSidx2+1; j<=SM2; j++) {
								if(SmeshType1 == 0) 			
									tU[0][j] = 2*tU[1][j] - tU[2][j];				
								else					
									tU[0][j] = (h1[0]+h1[1])/h1[1] * tU[1][j] - h1[0]/h1[1] * tU[2][j];
							}
						//////////////////////////////////////////////////////// newV 구하기!!
							for(j=startSidx2+1; j<=SM2; j++)
								knock1newV[0][j] = knock1oldV[0][j]* 1/(1+rd/(YFactor*intraTN)); //DCValue1; //0.0;
		 
							for(i=1; i<SM1; i++) {  //x축 고정시켜서.. 				 

								for(j=startSidx2+1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock1newV[i][startSidx2] =  knocknewV[i][startSidx2]; 
								for(j=startSidx2+1; j<SM2; j++) 								
									knock1newV[i][j] = vecV2[j] - YMAl[i][j]*knock1newV[i][j-1];					

								knock1newV[i][SM2-1] = knock1newV[i][SM2-1]/YMAd[i][SM2-1];
								for(j=SM2-2; j>=startSidx2; j--) 						
									knock1newV[i][j] = (knock1newV[i][j] - YMAu[i][j]*knock1newV[i][j+1])/YMAd[i][j];			
				
								if(SmeshType2 == 0)				
									knock1newV[i][SM2] = 2*knock1newV[i][SM2-1] - knock1newV[i][SM2-2];
								else
									knock1newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock1newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock1newV[i][SM2-2];
				
							} // for(i x축.. 고정 다 품

							for(j=startSidx2+1; j<=SM2; j++) {
								if(SmeshType1 == 0)  {
									knock1newV[0][j] = 2*knock1newV[1][j] - knock1newV[2][j];
									knock1newV[SM1][j] = 2*knock1newV[SM1-1][j] - knock1newV[SM1-2][j];
								}
								else {
									knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
									knock1newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock1newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock1newV[SM1-2][j];
								}
							} // knocknewV 완결


							
												//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock1newV[0][j] = 2*knock1newV[1][j] - knock1newV[2][j];
								 
							}
							else {
								knock1newV[0][j] = (h1[0]+h1[1])/h1[1] * knock1newV[1][j] - h1[0]/h1[1] * knock1newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock1newV[i][0] = 2*knock1newV[i][1] - knock1newV[i][2];
								 
							}
							else {
								knock1newV[i][0] = (h1[0]+h1[1])/h1[1] * knock1newV[i][1] - h1[0]/h1[1] * knock1newV[i][2];
								 
							}
						} // knocknewV 완결




							for(i=0; i<=SM1; i++)
								for(j=0; j<=startSidx2; j++)  									
									knock1newV[i][j] = knocknewV[i][j];
 					 
						} 
		//////////////////////////////// KI1 end
	 
 		//////////////////////////// KI2				//herework!!					
						if(m_ki1Flag == 0) {						 
						 
							for(i=startSidx1+1; i<=SM1; i++)
								tU[i][0] = knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); //DCValue2;

							for(j=1; j<SM2; j++) { //y축 고정 풀기						
								for(i=startSidx1+1; i<SM1; i++) 
									vecV1[i] = knock2oldV[i][j] + coeV1[i][j]*(-knock2oldV[i+1][j-1]+knock2oldV[i-1][j-1] +knock2oldV[i+1][j+1]-knock2oldV[i-1][j+1]); //herexxx
                      
								tU[startSidx1][j] = knocknewV[startSidx1][j];
								for(i=startSidx1+1; i<SM1; i++)
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];

								tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
								for(i=SM1-2; i>=startSidx1; i--) 											
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
													
								if(SmeshType1 == 0)									
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];							
								else								
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
							
						 			

							}// for(j=1

							for(i=startSidx1; i<=SM1; i++) {				
								if(SmeshType2 == 0) {
									//tU[i][0] = 2*tU[i][1] - tU[i][2];
									tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
								}
								else  {
									//tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];
									tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
								}
							} // 일차 tU 완결  
								
							//호
							for(i=startSidx1+1; i<=SM1; i++) {				
								if(SmeshType2 == 0) {
									tU[i][0] = 2*tU[i][1] - tU[i][2];												
								}
								else  {
									tU[i][0] = (h2[0]+h2[1])/h2[1] * tU[i][1] - h2[0]/h2[1] * tU[i][2];											
								}
							} // 일차 tU 완결  
						
						//////////////////////////////////////////////////////// newV 구하기!!
							for(j=0; j<=SM2; j++)
								knock2newV[startSidx1][j]=knocknewV[startSidx1][j];
		 
							for(i=startSidx1+1; i<SM1; i++) {  //x축 고정시켜서.. 				 

								for(j=1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);		//herexxx		   				 		 		

								knock2newV[i][0] = knock2oldV[i][0]* 1/(1+rd/(YFactor*intraTN)); // DCValue2; 
								for(j=1; j<SM2; j++) 								
									knock2newV[i][j] = vecV2[j] - knockYMAl[i][j]*knock2newV[i][j-1];					

								knock2newV[i][SM2-1] = knock2newV[i][SM2-1]/knockYMAd[i][SM2-1];
								for(j=SM2-2; j>=0; j--) 						
									knock2newV[i][j] = (knock2newV[i][j] - knockYMAu[i][j]*knock2newV[i][j+1])/knockYMAd[i][j];			
				
								if(SmeshType2 == 0) {
									//knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];
									knock2newV[i][SM2] = 2*knock2newV[i][SM2-1] - knock2newV[i][SM2-2];
								}
								else {
									//knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];
									knock2newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * knock2newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * knock2newV[i][SM2-2];
								}
				
							} // for(i x축.. 고정 다 품

							for(j=1; j<=SM2; j++) {
								if(SmeshType1 == 0) 
									knock2newV[SM1][j] = 2*knock2newV[SM1-1][j] - knock2newV[SM1-2][j];
								else
									knock2newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * knock2newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * knock2newV[SM1-2][j];
							} // knocknewV 완결

							
							//호
							for(i=startSidx1+1; i<=SM1; i++) {				
								if(SmeshType2 == 0) {
									knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];												
								}
								else  {
									knock2newV[i][0] = (h2[0]+h2[1])/h2[1] * knock2newV[i][1] - h2[0]/h2[1] * knock2newV[i][2];											
								}
							} // 일차 tU 완결  

										
												//modi						
					//for(j=startSidx2+1; j<=SM2; j++) {
						for(j=5; j<=SM2; j++) {
							if(SmeshType1 == 0)  {
								knock2newV[0][j] = 2*knock2newV[1][j] - knock2newV[2][j];
								 
							}
							else {
								knock2newV[0][j] = (h1[0]+h1[1])/h1[1] * knock2newV[1][j] - h1[0]/h1[1] * knock2newV[2][j];
								 
							}
						} // knocknewV 완결
										
					//for(i=startSidx1+1; i<=SM1; i++) {
						for(i=5; i<=SM1; i++) {
							if(SmeshType1 == 0)  {
								knock2newV[i][0] = 2*knock2newV[i][1] - knock2newV[i][2];
								 
							}
							else {
								knock2newV[i][0] = (h1[0]+h1[1])/h1[1] * knock2newV[i][1] - h1[0]/h1[1] * knock2newV[i][2];
								 
							}
						} // knocknewV 완결

							for(i=0; i<=startSidx1; i++)
								for(j=0; j<=SM2; j++)  									
									knock2newV[i][j] = knocknewV[i][j];
 					 
						} 
		//////////////////////////////// KI2 end
	 
						// no KI 영역 풀기 
						if(m_ki1Flag == 0 && m_ki2Flag == 0) {				
					 
							for(i=startSidx1; i<=SM1; i++)
								tU[i][startSidx2] = knock2newV[i][startSidx2];

							for(j=startSidx2+1; j<SM2; j++) { //y축 고정 풀기						
								for(i=startSidx1+1; i<SM1; i++) 
									vecV1[i] = oldV[i][j] + coeV1[i][j]*(-oldV[i+1][j-1]+oldV[i-1][j-1] +oldV[i+1][j+1]-oldV[i-1][j+1]); //herexxx
                        
								tU[startSidx1][j] = knock1newV[startSidx1][j];
								for(i=startSidx1+1; i<SM1; i++)
									tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
								tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
								for(i=SM1-2; i>=startSidx1; i--) 				
									tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
										
								if(SmeshType1 == 0)				
									tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
								else
									tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];				

							}// for(j=1

							for(i=startSidx1; i<=SM1; i++) {				
								if(SmeshType2 == 0) 					
									tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];				
								else					
									tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];			
							} // 일차 tU 완결 

							//////////////////////////////////////////////////////// newV 구하기!!
							for(j=startSidx2; j<=SM2; j++)
								newV[startSidx1][j]=knock1newV[startSidx1][j]; // D.C 적용
		 
							for(i=startSidx1+1; i<SM1; i++) {  //y축 고정시켜서.. 				 

								for(j=startSidx2+1; j<SM2; j++) 				
									vecV2[j] =tU[i][j] + coeV2[i][j] *(-tU[i-1][j+1]+tU[i-1][j-1] + tU[i+1][j+1] - tU[i+1][j-1]);	//herexxx				   				 		 		

								newV[i][startSidx2] = knock2newV[i][startSidx2];
								for(j=startSidx2+1; j<SM2; j++) 
									newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];					
								newV[i][SM2-1] = newV[i][SM2-1]/YMAd[i][SM2-1];
								for(j=SM2-2; j>=startSidx2; j--) 
									newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];					
				
								if(SmeshType2 == 0)				
									newV[i][SM2] = 2*newV[i][SM2-1] - newV[i][SM2-2];
								else
									newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * newV[i][SM2-2];
				
							} // for(j y축.. 고정 다 품

							for(j=startSidx2; j<=SM2; j++) {
								if(SmeshType1 == 0) 
									newV[SM1][j] = 2*newV[SM1-1][j] - newV[SM1-2][j];
								else
									newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * newV[SM1-2][j];
							} // newV 완결

							for(i=0; i<=startSidx1; i++)
								for(j=0; j<=startSidx2; j++) 								
										newV[i][j] = knocknewV[i][j];  // KI12

							for(i=0; i<=startSidx1; i++)
								for(j=startSidx2+1; j<=SM2; j++) 								
										newV[i][j] = knock1newV[i][j]; // KI1

							for(i=startSidx1+1; i<=SM1; i++)
								for(j=0; j<=startSidx2; j++) 								
										newV[i][j] = knock2newV[i][j]; // KI2
							
 					 
						}// if(m_kiFlag != 1)
 
							
						// KI 종류에 따라..
						if(m_ki1Flag == 1 && m_ki2Flag == 1) { 
							for(i=0; i<=SM1; i++)
								memcpy(newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
							for(i=0; i<=SM1; i++)
								memcpy(knock1newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
							for(i=0; i<=SM1; i++)
								memcpy(knock2newV[i],knocknewV[i],(size_t)((SM2+1)*sizeof(double)));
						}
		
						if(m_ki1Flag == 1 && m_ki2Flag == 0) { 
							for(i=0; i<=SM1; i++)
								memcpy(newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
							for(i=0; i<=SM1; i++)
								memcpy(knock2newV[i],knock1newV[i],(size_t)((SM2+1)*sizeof(double)));
						}
		
						if(m_ki1Flag == 0 && m_ki2Flag == 1) { 
							for(i=0; i<=SM1; i++)
								memcpy(newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
							for(i=0; i<=SM1; i++)
								memcpy(knock1newV[i],knock2newV[i],(size_t)((SM2+1)*sizeof(double)));
						}

 
						
						if(m_curtgreek == 1) { // 만기에
							for(i=0; i<outN + 4; i++) {		
								for(j=0; j<outN + 4; j++) {
																				
									sim_knock1Flag = 0;										
									sim_knock2Flag = 0;										
									if(SrTemp1[i] <= m_kihval1)								
										sim_knock1Flag = 1;										
									if(SrTemp2[j] <= m_kihval2)								
										sim_knock2Flag = 1;

									if(evalq<=1) { 							
										sleveltemp1 = SrTemp1[i]/m_bval1;
										sleveltemp2 = SrTemp2[j]/m_bval2;
									}
									else {
										sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN[curtN-1]-(evalq-1));
										for(evali=1; evali<evalq; evali++)
											sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
										sleveltemp1 = sleveltemp1/m_evalN[curtN-1];

										sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN[curtN-1]-(evalq-1));
										for(evali=1; evali<evalq; evali++)
											sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
										sleveltemp2 = sleveltemp2/m_evalN[curtN-1];
									}

																
									if(sim_knock1Flag == 1 && sim_knock2Flag == 1) 
										SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knocknewV, SM1, SM2);  
									if(sim_knock1Flag == 1 && sim_knock2Flag == 0) 
										SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock1newV, SM1, SM2);
									if(sim_knock1Flag == 0 && sim_knock2Flag == 1) 
										SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, knock2newV, SM1, SM2);
									if(sim_knock1Flag == 0 && sim_knock2Flag == 0)			
										SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2);
									
								}
							}		 
						}
					

						strcpy_s(fname,kCodeN,m_kCode);
						strcat_s(fname,kCodeN,"_");
						_itoa_s(intraLoop-10,tempstring,kCodeN,10);
						strcat_s(fname,kCodeN,tempstring);
						strcat_s(fname,kCodeN,".txt");

						batchsf.open(fname,ios::out);
					
						batchsf << fixed;
						batchsf.precision(15);
									 
						for(i=0; i<outN + 4; i++) {
							batchsf << SrTemp1[i] << " ";
						}
						batchsf << endl;
			
						for(i=0; i<outN + 4; i++) {
							batchsf << SrTemp2[i] << " ";
						}
						batchsf << endl;

						for(i=0; i<outN + 4; i++)  {
							for(j=0; j<outN + 4; j++) {
								batchsf << SrV1Temp[i][j] << " ";
							}
							batchsf << endl;
						}

						batchsf.close();

					}//for(intraLoop)

					delete tempstring;
				} //if(m_matu == 0)
	 
				delete fname;
			
			}// if(m_batchFlag)

			free(SrTemp1);
			free(SrTemp2);
			for(i=0; i<outN+4; i++){		
				free(SrVTemp[i]);		
				free(SrV1Temp[i]);
			}
			free(SrVTemp);
			free(SrV1Temp);
 
		}

		if(m_curtgreek == 2) { // vega 
			value = interp2(slevel1,slevel2,S1,S2,newV,SM1,SM2);
			for(i=0; i<3; i++) {
				greeks[i] = SrVTemp[Baseidx + (i-1)][Baseidx]; // greeks[0]: s-, greeks[1]: s0, greeks[2]: s+
			}
			for(i=3; i<6; i++) {
				greeks[i] = SrVTemp[Baseidx][Baseidx + (i-1-3)]; // greeks[3]: s-, greeks[4]: s0, greeks[5]: s+
			}
		
			// free 전에.. file이나 memory로 출력 필요!!
			// 단, file 저장시  v+ 인지 v-인지, i얼마인지 필요


			if(m_batchFlag) {

				using namespace std;
				using std::ofstream;

				ofstream batchsf;
			
				char *fname;
				char *tempstring;
				fname = new char[kCodeN+1];
				tempstring = new char[kCodeN+1];


				strcpy_s(fname,kCodeN,m_kCode);			 			
				if(m_curtvegapertubation == 0) {
					strcat_s(fname,kCodeN,"_sigma1+_");
				}
				if(m_curtvegapertubation == 1) {			
					strcat_s(fname,kCodeN,"_sigma1-_");
				}			
				if(m_curtvegapertubation == 2) {
					strcat_s(fname,kCodeN,"_sigma2+_");
				}
				if(m_curtvegapertubation == 3) {
					strcat_s(fname,kCodeN,"_sigma2-_");
				}

				_itoa_s(m_curtvegaidx,tempstring,kCodeN,10);			
				strcat_s(fname,kCodeN,tempstring);			
				strcat_s(fname,kCodeN,".txt");
		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 
											 
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp1[i] << " ";
				}
				batchsf << endl;
			
				for(i=0; i<outN + 4; i++) {
					batchsf << SrTemp2[i] << " ";
				}
				batchsf << endl;

				for(i=0; i<outN + 4; i++)  {
					for(j=0; j<outN + 4; j++) {
						batchsf << SrVTemp[i][j] << " ";
					}
					batchsf << endl;
				}

				batchsf.close();

 				delete tempstring;
				delete fname;
			
			}// if(m_batchFlag)



			free(SrTemp1);
			free(SrTemp2);
			for(i=0; i<outN+4; i++){		
				free(SrVTemp[i]);		
				free(SrV1Temp[i]);
			}
			free(SrVTemp);
			free(SrV1Temp);
 
		}

		if(m_curtgreek == 3) { // rho
			value = interp2(slevel1,slevel2,S1,S2,newV,SM1,SM2); // value - m_Ogreeks[0] 을 이용하기 위하여
		}
	
		if(m_curtgreek == 4) { // CorrelationDelta
			value = interp2(slevel1,slevel2,S1,S2,newV,SM1,SM2); // value - m_Ogreeks[0] 을 이용하기 위하여
		}
	 
 
		free(S1);
		free(h1);
		free(S2);
		free(h2);
	
		for(i=0; i<=SM1; i++) {	
			free(oldV[i]);
			free(newV[i]);
			free(tU[i]);

			free(knockoldV[i]);
			free(knocknewV[i]);		 

			//SH
			free(knock1oldV[i]);
			free(knock1newV[i]);		 
			free(knock2oldV[i]);
			free(knock2newV[i]);		 

			free(XMAl[i]);
			free(XMAd[i]);
			free(XMAu[i]);

			free(YMAl[i]);
			free(YMAd[i]);
			free(YMAu[i]);			 

			free(knockXMAl[i]);
			free(knockXMAd[i]);
			free(knockXMAu[i]);

			free(knockYMAl[i]);
			free(knockYMAd[i]);
			free(knockYMAu[i]);

			free(coeV1[i]);
			free(coeV2[i]);
		}
		free(oldV);
		free(newV);
		free(tU);
		
		free(knockoldV);
		free(knocknewV);

		//SH
		free(knock1oldV);
		free(knock1newV);
		free(knock2oldV);
		free(knock2newV);

	
		free(XMAu);
		free(XMAd);
		free(XMAl);
		
		free(knockXMAu);
		free(knockXMAd);
		free(knockXMAl);
	
		free(coeV1);
		free(vecV1);	

		free(YMAu);
		free(YMAd);
		free(YMAl);
	
		free(knockYMAu);
		free(knockYMAd);
		free(knockYMAl);

		free(coeV2);
		free(vecV2);

		free(volskew1);
		free(volskew2);


		return value;
	 
} // fd_price();



void HStepDownB::curttimeVol1(double *vol, int day)
{
	int i;
	
	if(m_voltype == 1) { // term imp vol
		double v1, v2;
		int t1, t2;
		double *termvol;

		termvol = new double[m_volTN];
		for(i=0; i<m_volTN; i++) 
			termvol[i] = m_vol1[0][i];
		
		t1 = day;
		v1 = interp1(t1,m_volBDay,termvol,m_volTN,1,0);
		t2 = day+1;
		v2 = interp1(t2,m_volBDay,termvol,m_volTN,1,0);

		vol[0] = sqrt(MAX(minVol*minVol, (v2*v2*t2 - v1*v1*t1)/(t2-t1))); // implied forward vol

		delete termvol;
	}
	else { // step function 형 local vol
		int findidx;
		
		findidx = m_volTN-1;
		for(i=0; i<m_volTN; i++) {
			if(day <= m_volBDay[i]) {
				findidx = i;
				break;
			}
		}		 

		for(i=0; i<m_volSN; i++) 
			vol[i] = m_vol1[i][findidx];
	}
	
}

void HStepDownB::curttimeVol2(double *vol, int day)
{
	int i;
	
	if(m_voltype == 1) { // term imp vol
		double v1, v2;
		int t1, t2;
		double *termvol;

		termvol = new double[m_volTN];
		for(i=0; i<m_volTN; i++) 
			termvol[i] = m_vol2[0][i];
		
		t1 = day;
		v1 = interp1(t1,m_volBDay,termvol,m_volTN,1,0);
		t2 = day+1;
		v2 = interp1(t2,m_volBDay,termvol,m_volTN,1,0);

		vol[0] = sqrt(MAX(minVol*minVol, (v2*v2*t2 - v1*v1*t1)/(t2-t1))); // implied forward vol

		delete termvol;
	}
	else { // step function 형 local vol
		int findidx;
		
		findidx = m_volTN-1;
		for(i=0; i<m_volTN; i++) {
			if(day <= m_volBDay[i]) {
				findidx = i;
				break;
			}
		}		 

		for(i=0; i<m_volSN; i++) 
			vol[i] = m_vol2[i][findidx];
	}
	
}

double HStepDownB::curttimeCorr(int day)
{
	int i;		
	int findidx;
		
	findidx = m_corrN-1;
	for(i=0; i<m_corrN; i++) {
		if(day <= m_corrBDay[i]) {
			findidx = i;
			break;
		}
	}		 

	return m_corr[findidx];
	 
	
}

int HStepDownB::findinteridx(double x, double *Dx, int DxN)
{
	int i;		
	int findidx;
	
	findidx = DxN-1;
	for(i=0; i<DxN; i++) {
		if(Dx[i] < x && x <= Dx[i+1]) {
			findidx = i+1;
			break;
		}
	}

	return findidx;	
}

 