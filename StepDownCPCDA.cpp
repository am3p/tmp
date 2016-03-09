#include "StepDownCPCDA.h"
#include <iostream>
#include <fstream>

 

StepDownCPCDA::StepDownCPCDA(int kiFlag, double sval, double bval, double *xval, double kihval, 
	      double *Cpn,
	      int *matuN, int *matu, double ctime,
		  int maxevalN, int *evalN,  double *psval,
		  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,
		  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int batchFlag, char *kCode,
		  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
		  double shiftHRate, double spreadHRate, double *spreadXRate)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	int i, j;	 

	//입력 데이터 내부 전역 변수로 저장

	m_matuN = matuN[1];
	m_matu = new int[m_matuN];
	m_sval = sval;
	m_bval = bval;	
	m_xval = new double[m_matuN];
	
	m_kihval = kihval;
	m_kiFlag = kiFlag;
	m_Cpn = new double[m_matuN+2];


	
	m_paymatuN = matuN[0];
	m_paymatu = new int[m_paymatuN];
	m_payxval = new double[m_paymatuN];
	m_payCpn = new double[m_paymatuN];


	m_shiftHRate = shiftHRate;
	m_spreadHRate = abs(spreadHRate);

	m_spreadXRate = new double[m_matuN];
	//m_spreadXS = new double[m_matuN];

	for(i=0; i<m_paymatuN; i++) {
		m_payxval[i] = xval[i];
		m_payCpn[i] = Cpn[i];
		m_paymatu[i] = matu[i];
	}

	for(i=0; i<m_matuN; i++) {
		m_xval[i] = xval[m_paymatuN+i];
		m_Cpn[i] = Cpn[m_paymatuN+i];
		m_matu[i] = matu[m_paymatuN+i];
		m_spreadXRate[i] = abs(spreadXRate[i]);				
	}
	m_Cpn[m_matuN] = Cpn[m_paymatuN+m_matuN]; // dummy cpn
	m_Cpn[m_matuN+1] = Cpn[m_paymatuN+m_matuN+1]; // 원금보장형 조정 cpn
		 
	m_ctime = ctime;

	m_maxevalN = maxevalN;
	m_evalN = new int[m_matuN];
	for(i=0; i<m_matuN; i++) 
		m_evalN[i] = evalN[i];
	m_psval = new double[m_maxevalN];
	for(i=0; i<m_maxevalN; i++)
		m_psval[i] = psval[i];	 
	
	//cd
	m_CDN = CDN;
	m_CDadjval = CDadjval;
	m_ObBDay = new int[m_CDN];
	m_PayBDay = new int[m_CDN];
	m_FlegFactor = new double[m_CDN];
	m_ObCDrate = new double[m_CDN];

	for(i=0; i<m_CDN; i++) {
		m_ObBDay[i] = ObBDay[i];
		m_PayBDay[i] = PayBDay[i];
		m_FlegFactor[i] = FlegFactor[i];
		m_ObCDrate[i] = ObCDrate[i];
	}


	m_irateN = irateN;
	m_irateBDay = new int[m_irateN];
	 
	m_iRFrate = new double[m_irateN];
	m_iDCrate = new double[m_irateN];
	m_iIRSrate = new double[m_irateN];
	for(i=0; i<m_irateN; i++) {
		m_irateBDay[i] = irateBDay[i];
	 
		m_iRFrate[i] = iRFrate[i];
		m_iDCrate[i] = iDCrate[i];
		m_iIRSrate[i] = iIRSrate[i];
	}

	m_divrate = divrate;
	m_divbval = divbval;
	m_divN = divN;
	m_divBDay = new int[m_divN];
	m_divCash = new double[m_divN];
	for(i=0; i<m_divN; i++) {
		m_divBDay[i] = divBDay[i];
		m_divCash[i] = divCash[i];
	}
	m_divApy = divApy;

	m_voltype = voltype;
	m_volbval = volbval;
	m_volTN = volTN;
	m_volSN = volSN;
	m_volBDay = new int[m_volTN];
	for(i=0; i<m_volTN; i++)
		m_volBDay[i] = volBDay[i];
	m_volSmness = new double[m_volSN];
	for(i=0; i<m_volSN; i++)
		m_volSmness[i] = volSmness[i];
	m_vol = new double*[m_volSN]; // m_volSN은.. 1이상
	for(i=0; i<m_volSN; i++)
		m_vol[i] = new double[m_volTN];
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[i][j] = vol[i*m_volTN+j];

	m_batchFlag = batchFlag;
	m_kCode = new char[kCodeN+1];
	strcpy_s(m_kCode,kCodeN,kCode);


	m_Sminlevel = Sminlevel;
	m_Smaxlevel = Smaxlevel;
	
	
	m_meshSM1 = SmeshM[0];
	m_meshSM2 = SmeshM[1];
	m_meshSM3 = SmeshM[2];

	m_meshTN1 = TmeshN[0];
	m_meshTN2 = TmeshN[1];


	 
	m_kihval = m_kihval + m_bval*m_shiftHRate;
	m_spreadHS = m_kihval - m_bval*m_spreadHRate;

	/*
	for(i=0; i<m_matuN; i++) 	  
		m_spreadXS[i] = m_xval[i] - m_bval*m_spreadXRate[i];
		*/
	
	m_startS = m_spreadHS;
	//m_endS = m_xval[m_matuN-2]; // 마지막 조기상환일에 필요 조기상환일에 세팅!! 

	 

	// ki 체크
	if(m_sval <= m_kihval)			
		m_kiFlag = 1;
	 

	int debugFlag;
	debugFlag = 0;
	if(debugFlag) {
		using namespace std;
		using std::ofstream;

		ofstream logs;
		logs.open("c:/logland/stepdowncdainput.txt",ios::out);
			
		logs << "m_sval = " << m_sval << endl;	 
	 
		for(i=0; i<m_irateN; i++) {
			logs << "m_irateBDay[" << i <<"]" <<"= "<< m_irateBDay[i] << endl;
		} 

		logs.close();
	}

	  
}

StepDownCPCDA::~StepDownCPCDA()
{	
	int i;

	delete m_evalN;
	delete m_psval;

	//cd
	delete m_ObBDay;
	delete m_PayBDay;
	delete m_FlegFactor;
	delete m_ObCDrate;


	delete m_irateBDay;
	 
	delete m_iRFrate; 
	delete m_iDCrate;  
	delete m_iIRSrate; 	 
			
	delete m_divBDay;
	delete m_divCash;
			
	delete m_volBDay;	 
	delete m_volSmness;	 

	for(i=0; i<m_volSN; i++)
		delete m_vol[i];
	delete m_vol;

	delete m_matu;
	delete m_xval;
	delete m_Cpn;
	delete m_spreadXRate;
	//delete m_spreadXS;

	
	delete m_paymatu;
	delete m_payxval;
	delete m_payCpn;

	delete m_kCode;

}

double StepDownCPCDA::getValue() {

	double greeks[MAXGREEKS1];
	double value;
	int i;

	m_curtgreek = 1;	
	value = fd_price(greeks);
	

	for(i=0;i<BasicGRsN1; i++) {
		m_Ogreeks[i] = greeks[i];
	}
 

	return value;

}


double StepDownCPCDA::getVega(int vegaN, int *vegaidx) {

	double greeks[MAXGREEKS1];
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
			tempv[i][j] = m_vol[i][j]; // vol 저장
	

	m_Ogreeks[Vegaidx] = 0;
	m_Ogreeks[VannaCidx] = 0;
	m_Ogreeks[ZommaCidx] = 0;

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
						m_vol[k][j] = tempv[k][j] + curtpertubation;
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
							m_vol[k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[k][j] = tempv[k][j] + curtpertubation;
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
		m_Otvegas[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/

		//2012-05-23 1% vanna cash / 1% zomma cash 수정
		m_Otvegas[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/4; // term 1% vanna cash
		m_Otzommas[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/2; // term 1% zomma cash

		m_Ogreeks[Vegaidx] += m_Otvegas[i];
		m_Ogreeks[VannaCidx] += m_Otvannas[i];
		m_Ogreeks[ZommaCidx] += m_Otzommas[i];
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[i][j] = tempv[i][j]; // vol 복원
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;

	return m_Ogreeks[Vegaidx];

}


double StepDownCPCDA::getRho(int rhoN,int *rhoidx) {

	double greeks[MAXGREEKS1];
	double value;

	double *tempirate;
	double *tempirs;
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		 

	m_curtgreek = 3;
 


	tempirate = new double[m_irateN];
	tempirs = new double[m_irateN]; //cd
	//rf rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iRFrate[i]; // RF rate 저장
		tempirs[i] = m_iIRSrate[i]; //cd
	}
	 
	m_Ogreeks[Rhorfidx] = 0; //rf rho sum
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // 평행이동 rho
			for(j=0; j<m_irateN; j++) {
				m_iRFrate[j] = tempirate[j] + RhoPertubation;		
				m_iIRSrate[j] = tempirs[j] + RhoPertubation;
			}
		}
		else {
			if(i == 0) {		
				pertubationX[0] = m_irateBDay[rhoidx[i]];
				pertubationX[1] = m_irateBDay[rhoidx[i+1]];

				pertubationY[0] = RhoPertubation;
				pertubationY[1] = 0.0;

				for(j=0; j<m_irateN; j++) {			
					m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
					m_iIRSrate[j] = tempirs[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
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
						m_iIRSrate[j] = tempirs[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
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
						m_iIRSrate[j] = tempirs[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat
					}	

				}
			}
		}// if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otrfrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[Rhorfidx] += m_Otrfrhos[i];
	}
 

	for(i=0; i<m_irateN; i++) {
		m_iRFrate[i] = tempirate[i]; // RF rate 복원
		m_iIRSrate[i] = tempirs[i]; // RF rate 복원
	}
	

//////////////////////////////
	//dc rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iDCrate[i]; // DC rate 저장
	}
	 
	m_Ogreeks[Rhodcidx] = 0; // dc rho sum
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
		m_Ogreeks[Rhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate 복원
	}
 	


	delete tempirate;



	return m_Ogreeks[Rhorfidx] + m_Ogreeks[Rhodcidx];

}

double StepDownCPCDA::cf_price()
{	
	return 7777;
}

double StepDownCPCDA::fd_price(double *greeks)
{

	if(m_modeltype == 10 || m_modeltype == 11)
		m_kiFlag = 1;

	double DCValue, CDlegDCValue;	

	int i, j, k, matuLoop, intraLoop;
	int SM, SmeshType, intraTN;

	double ds, ds3, dt;

	double slevel, xlevel,hlevel, spreadHSlevel;//, spreadXSlevel;
	double *S, *h;
	double *oldV, *newV;
	double *knockoldV, *knocknewV;

	//cd
	double *CDoldV, *CDnewV;
	double *fdObCDrate, *FlegCpn, *CpnCDleg;

	double *MAu, *MAd, *MAl;
	double *knockMAu, *knockMAd, *knockMAl;

	double startSlevel, endSlevel;

	int xidx;
	int startSidx, endSidx;
	int spreadHSidx;//, spreadXSidx;
	int hidx;
	int evali;
	int evalq, expiryq;


	double value;

	double *SrTemp, *SrVTemp, *SrV1Temp; // greek 계산을 위한 변수 당일값, 다음날 값 저장
	double sleveltemp;		
	int outN;	
	int Baseidx;
	int curtN;
	 
	double cdpayvalue;
	double disc_div;
	double Srv;
	double *volskew;

	int IsMidCompu, mid_n;


	fdObCDrate = new double[m_CDN];
	FlegCpn = new double[m_CDN];
	CpnCDleg = new double[m_matuN]; //조기 상환일 CD 지급 cash
	
	double cdtempval;
	double IFRtemp;

	for(i=0; i<m_CDN; i++) {
		FlegCpn[i] = 0.0;
		if(m_ObBDay[i] <= 0 && m_ObCDrate[i] <= 0)
			fdObCDrate[i] = m_iIRSrate[0];
		else
			fdObCDrate[i] = m_ObCDrate[i];
		if(m_PayBDay[i] >= 0) {
			//if(i==0 || m_ObBDay[i] <= 0)
			if(m_ObBDay[i] <= 0)
				FlegCpn[i] = (fdObCDrate[i]+m_CDadjval)*m_FlegFactor[i];
			else {
				//cdtempval = (interp1(m_PayBDay[i],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_PayBDay[i] - interp1(m_PayBDay[i-1],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_PayBDay[i-1])/(m_PayBDay[i]-m_PayBDay[i-1]);
				if(i == m_CDN-1)
					cdtempval = (interp1(m_PayBDay[i],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_PayBDay[i] - interp1(m_ObBDay[i],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_ObBDay[i])/(m_PayBDay[i]-m_ObBDay[i]);
				else				
					cdtempval = (interp1(m_ObBDay[i+1],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_ObBDay[i+1] - interp1(m_ObBDay[i],m_irateBDay,m_iIRSrate,m_irateN,1,0)*m_ObBDay[i])/(m_ObBDay[i+1]-m_ObBDay[i]);
				
				FlegCpn[i] = (cdtempval+m_CDadjval)*m_FlegFactor[i];
			}
		}

		if(i != m_CDN-1 && m_PayBDay[i] == 0 && m_ctime > 15)
			FlegCpn[i] = 0; // cd 받은것 cash화 된..
	}

	CpnCDleg[m_matuN-1] = FlegCpn[m_CDN-1];
	for(i=0; i<m_matuN-1; i++) {
		CpnCDleg[i] = 0.0;
		if(m_matu[i] >= 0) {
			for(j=0; j<m_CDN-1; j++) {
				if(m_matu[i]<=m_PayBDay[j]) {
					IFRtemp = interp1(m_PayBDay[j],m_irateBDay,m_iDCrate,m_irateN,1,0)*m_PayBDay[j] - interp1(m_matu[i],m_irateBDay,m_iDCrate,m_irateN,1,0)*m_matu[i];
					//CpnCDleg[i] = FlegCpn[j]*exp(-IFRtemp/YFactor);
					if (m_FlegFactor[j] == 0)
						CpnCDleg[i] = 0.0;
					else
						CpnCDleg[i] = FlegCpn[j]/m_FlegFactor[j] * (m_matu[i]-m_ObBDay[j])/YFactor;
					break;
				}
			}
		}
	}
	//cd
	double *temppayCpn;
			
	temppayCpn = new double[m_paymatuN];

	for(i=0; i<m_paymatuN; i++) {
		temppayCpn[i] = m_payCpn[i];

		if(i != m_paymatuN-1 && m_paymatu[i] == 0 && m_ctime > 15)
			temppayCpn[i] = 0;
	}


	// 단순 초기화
	i = 0;
	j = 0;
	k = 0;
	matuLoop = 0;
	intraLoop = 0;
	

	curtN = m_matuN;

	evalq = m_evalN[curtN-1] - m_matu[curtN-1];
	
	if(evalq <=1 ) 	
		slevel = m_sval/m_bval;	
	else {
		slevel = m_sval/m_bval * (m_evalN[curtN-1] - (evalq-1));
		for(evali=1; evali<evalq; evali++)
			slevel = slevel + m_psval[evali]/m_bval;
		slevel = slevel/m_evalN[curtN-1];
	}	 
	 


	xlevel = m_xval[curtN-1]/m_bval;
	hlevel = m_kihval/m_bval;
	startSlevel = m_startS/m_bval;
	//endSlevel = m_endS/m_bval;     //<<-- 만기에서는 필요 없음
	spreadHSlevel = m_spreadHS/m_bval;
	//spreadXSlevel = m_spreadXS[curtN-1]/m_bval;

	if(xlevel == hlevel)
		m_kiFlag = 1;
	 
	double meshxlevel;
	for(i=0; i<m_matuN; i++) {
		if(m_matu[i]>=0) {
			meshxlevel = m_xval[i]/m_bval;			
			break;
		}
	}
	 
			
	ds = (m_Smaxlevel - m_Sminlevel)/m_meshSM1;
	xidx = (int)floor((meshxlevel-m_Sminlevel)/ds);
	 
		
	// Adaptive mesh 조건 만족 못할 때는 균등 분할 적용
	if((m_meshSM2 <=0) || (m_meshSM3 <=0) || ds <= AdmeshR ) { // || (xidx1 - m_meshSM2 < 0) || (m_meshSM1 < xlevel1 + m_meshSM2)) { 
		SmeshType = 0; // 균등
		SM = m_meshSM1;
		 
		S = (double *)malloc((size_t)(SM+1)*sizeof(double));
		//S mesh generate
		for(i=0; i<=SM; i++) 
			S[i] = m_Sminlevel + ds*i;
		h = (double *)malloc((size_t)(1)*sizeof(double)); // free(h)를 위한 더미 설정
		h[0] = ds;
		
	}
	else {		
		SmeshType = 1; // Adaptive mesh
		if(m_meshSM3 == 1) {
			double *tempgenS;
			tempgenS = new double[m_meshSM1 + int(4*Admesh/AdmeshR)+4];

			int cri;
			double temph, temps;
			
			cri = 0;
			tempgenS[cri] = m_Sminlevel;
			temph = ds;
			temps = tempgenS[cri] + temph;
			while (temps <= m_Smaxlevel + 2*myeps) {
				if(temps <= hlevel - Admesh + myeps)
					temph = ds;
				if(temps > hlevel - Admesh + myeps )
					temph = AdmeshR;
				if(temps > hlevel + Admesh + myeps)
					temph = ds;
				if(temps > meshxlevel - Admesh + myeps)
					temph = AdmeshR;
				if(temps > meshxlevel + Admesh + myeps)
					temph = ds;

				temps = tempgenS[cri] + temph;
				cri = cri + 1;
				tempgenS[cri] = temps;
				temps = tempgenS[cri] + temph;		
				if(tempgenS[cri] > m_Smaxlevel)
					tempgenS[cri] = m_Smaxlevel;
			}	

			SM = cri;
			S = (double *)malloc((size_t)(SM+1)*sizeof(double));
			for(i=0; i<=SM; i++)
				S[i] = tempgenS[i];

			delete tempgenS;			
		}
		else {
			SM = m_meshSM1 + 2*m_meshSM2*(m_meshSM3-1);
			ds3 = ds/m_meshSM3;

			S = (double *)malloc((size_t)(SM+1)*sizeof(double));
			//S mesh generate
			for(i=0; i<=m_meshSM1; i++) {
				if(i<=xidx-m_meshSM2) {
					S[i] = m_Sminlevel + ds*i;
				}
				else {
					if(xidx-m_meshSM2 < i && i<= xidx + m_meshSM2) {
						for(j=1; j<=m_meshSM3; j++) {					
							S[xidx-m_meshSM2 + m_meshSM3*(i-(xidx-m_meshSM2+1)) +j] = S[xidx-m_meshSM2] + ds3*(m_meshSM3*(i-(xidx-m_meshSM2+1))+j);
						}
					}
					else {
						S[i+2*m_meshSM2*(m_meshSM3-1)] = S[xidx+m_meshSM2*(2*m_meshSM3-1)] + ds*(i-(xidx+m_meshSM2));
					}
				}
			}
			
		}

		h = (double *)malloc((size_t)(SM)*sizeof(double));
		for(i=0; i<SM; i++) 
			h[i] = S[i+1] - S[i];

	}
 

	for(i=0;i<SM; i++) {
		if(S[i] < startSlevel && startSlevel <= S[i+1]) {
			startSidx = i+1;
			break;
		}
	}

	/*
	for(i=0; i<SM; i++) {
		if(S[i] < endSlevel && endSlevel <= S[i+1]) {
			endSidx = i+1;
			break;
		}
	}
	*/
	
	for(i=0; i<SM; i++) {
		if(S[i] < spreadHSlevel && spreadHSlevel <= S[i+1]) {
			spreadHSidx = i+1;
			break;
		}
	}

		
	/*
	for(i=0; i<SM; i++) {
		if(S[i] < spreadXSlevel && spreadXSlevel <= S[i+1]) {
			spreadXSidx = i+1;
			break;
		}
	}
	*/
	
	for(i=0; i<SM; i++) {
		if(S[i] < hlevel && hlevel <= S[i+1]) {
			hidx = i+1;
			break;
		}
	}

	
	for(i=0; i<SM; i++) {
		if(S[i] < xlevel && xlevel <= S[i+1]) {
			xidx = i+1;
			break;
		}
	}
	

	// greek 계산을 위한 작업
	if(m_curtgreek == 1 || m_curtgreek == 2) { // 기본 그릭		 		 

		outN = 2*SoutRange + 1; // +- 20%출력
		SrTemp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // up down gamma를 위해 위 아래 2개씩 더 
		SrVTemp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // 당일 저장
		SrV1Temp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // 다음날 저장

		for(i=0; i<outN + 4; i++) {
			SrTemp[i] = m_sval*(1-(SoutRange+2)*0.01 + 0.01*i); //m_sval의 78%~122%까지 1%씩 저장 0~44 45개 중심 22			
		}		 
	}




	
	CDoldV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	CDnewV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	
	oldV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	newV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAu = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAd = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAl = (double *)malloc((size_t)(SM+1)*sizeof(double));
	
	knockoldV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	knocknewV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	knockMAu = (double *)malloc((size_t)(SM+1)*sizeof(double));
	knockMAd = (double *)malloc((size_t)(SM+1)*sizeof(double));
	knockMAl = (double *)malloc((size_t)(SM+1)*sizeof(double));
	
	volskew = (double*)malloc((size_t)(m_volSN)*sizeof(double));
 
	 

	// 만기 payoff setting
	
	expiryq = 1;

	for(i=0; i<=SM; i++) {
		if(i>=xidx) {
			knocknewV[i] = AMOUNT * (m_Cpn[curtN-1]); //(1+m_Cpn[curtN-1]);
			newV[i] = knocknewV[i];
		}
		else {
			if(m_modeltype == 0 || m_modeltype == 10) 			
				knocknewV[i] = AMOUNT * (-1.0 + MIN(S[i],1));   //S[i];			AMOUNT * (-1.0 + S[i]);

			if(m_modeltype == 1 || m_modeltype == 11)
				knocknewV[i] = AMOUNT * (m_Cpn[curtN+1]); //(1+m_Cpn[curtN+1]);

			if(i>hidx) 
				newV[i] = AMOUNT * (m_Cpn[curtN]);   //(1+m_Cpn[curtN]); // dummy cpn 저장
			else
				newV[i] = knocknewV[i];
		}

		CDnewV[i] = AMOUNT * CpnCDleg[curtN-1]; // CD 는 전구간 동일
	}

	/*
	if(m_spreadHRate > 0) { // KI 쪽 spread 만기에만 
		double spreadvalue;

		for(i=spreadHSidx; i<=hidx; i++) {
			//spreadvalue = AMOUNT * (((1+m_Cpn[curtN])-0)/(hlevel - spreadHSlevel) * (S[i]-spreadHSlevel) + 0);
			spreadvalue = AMOUNT * (((1+m_Cpn[curtN])-0)/(S[hidx] - S[spreadHSidx]) * (S[i]-S[spreadHSidx]) + 0);
			//newV[i] = MAX(newV[i],spreadvalue);
			newV[i] = MAX(newV[i],spreadvalue-AMOUNT);
		}
	}
	*/

	if(m_spreadXRate[curtN-1] > 0){// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
		double spreadvalue;

		for(i=0; i<=xidx; i++) {
			//spreadvalue = AMOUNT * (((1+m_Cpn[curtN-1])-0)/(xlevel - spreadXSlevel) * (S[i]-spreadXSlevel) + 0);
			spreadvalue = AMOUNT * (m_spreadXRate[curtN-1] * (S[i]-S[xidx]) + 1+m_Cpn[curtN-1]);
			//knocknewV[i] = MAX(knocknewV[i],spreadvalue);
			//newV[i] = MAX(newV[i],spreadvalue);
			knocknewV[i] = MAX(knocknewV[i],spreadvalue-AMOUNT);
			newV[i] = MAX(newV[i],spreadvalue-AMOUNT);
		}
	}	 
	
	
		double payxlevel;
		payxlevel = m_payxval[m_paymatuN-1]/m_bval;
		for(i=0; i<=SM; i++) {
			if(S[i] >= payxlevel) {
				knocknewV[i] += AMOUNT * (temppayCpn[m_paymatuN-1]);
				newV[i] += AMOUNT * (temppayCpn[m_paymatuN-1]);
			}
		}
	 
	 
 
	
	if(m_kiFlag == 1) {	// ki or ko 되었으면..
		memcpy(newV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
	}
	
	int sim_knockFlag;
	if(m_curtgreek == 1 && (m_matu[curtN-1] == 1 || m_matu[curtN-1] == 0)) { // 기본 그릭	Theta를 위한 +1일 저장
		for(i=0; i<outN + 4; i++) {												
			sim_knockFlag = 0;
			if(SrTemp[i] <= m_kihval)
				sim_knockFlag = 1;

			if(m_matu[curtN-1] == 1) {
				if(m_ctime>15) { // 장 종료 후 는.. 오늘 종가가 내일 psval의 [0]에 들어 가는 경우
					sleveltemp = SrTemp[i]/m_bval; // *(m_evalN[curtN-1] - evalq) = 1 = m_matu[curtN-1] 해당
					for(evali=0; evali<evalq; evali++)
						sleveltemp = sleveltemp + m_psval[evali]/m_bval;
					sleveltemp = sleveltemp / m_evalN[curtN-1];
				}
				else { // 장중
					sleveltemp = (1+MIN(evalq,1))*SrTemp[i]/m_bval;					
					for(evali=1; evali<evalq; evali++)
						sleveltemp = sleveltemp + m_psval[evali]/m_bval;
					sleveltemp = sleveltemp / m_evalN[curtN-1];
				}
			}
			else {
				sleveltemp = SrTemp[i]/m_bval;
				for(evali=1; evali<evalq; evali++)
					sleveltemp = sleveltemp + m_psval[evali]/m_bval;
				sleveltemp = sleveltemp / m_evalN[curtN-1];
			}

			if(sim_knockFlag == 1) {
				SrV1Temp[i] = interp1(sleveltemp,S,knocknewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
			}
			else {			 
				SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
			}
		}		 
	}
 	 

	double r, q, v, rd;
	double *dailyIFrf, *dailyIFdc;	
	// intra Day에서는 입력 파라미터 고정
	// 하루마다 해당 파라미터 계산
	q = m_divrate;

	if(m_matu[curtN-1]>0) { // 현재 curtN 는 m_matuN <-- 만기해당..
		dailyIFrf = (double *)malloc((size_t)(m_matu[curtN-1])*sizeof(double));
		dailyIFdc = (double *)malloc((size_t)(m_matu[curtN-1])*sizeof(double));
		
		for(i=0; i<m_matu[curtN-1]; i++) {
			dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
			dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
		}
	} //if(m_matu[curtN-1]>0)
	
	 


	// time.. 계산 시작

	DCValue = knocknewV[0];
	CDlegDCValue = CDnewV[0];
	
	for(matuLoop=m_matu[curtN-1]-1; matuLoop>=0; matuLoop--) {
		//time Adaptive mesh
		mid_n = IsinArr(matuLoop, m_matuN, m_matu); // 조기상환일 중에 하나면 1 ~ m_matuN-1 값 리턴
		if(mid_n > 0) { //
			curtN = mid_n;
			xlevel = m_xval[curtN-1]/m_bval;
 
			//n 영업일 평균 관련 처리
			evalq = m_evalN[curtN-1] - m_matu[curtN-1];
	
			if(evalq <=1 ) 	
				slevel = m_sval/m_bval;	
			else {
				slevel = m_sval/m_bval * (m_evalN[curtN-1] - (evalq-1));
				for(evali=1; evali<evalq; evali++)
					slevel = slevel + m_psval[evali]/m_bval;
				slevel = slevel/m_evalN[curtN-1];
			}	 

			m_endS = m_xval[curtN-1]; 
			endSlevel = m_endS/m_bval;

			//spreadXSlevel = m_spreadXS[curtN-1]/m_bval;

			for(i=0; i<SM; i++) {
				if(S[i] < endSlevel && endSlevel <= S[i+1]) {
					endSidx = i+1;
					break;
				}
			}
		
			/*
			for(i=0; i<SM; i++) {
				if(S[i] < spreadXSlevel && spreadXSlevel <= S[i+1]) {
					spreadXSidx = i+1;
					break;
				}
			}
			*/
	
			for(i=0; i<SM; i++) {
				if(S[i] < xlevel && xlevel <= S[i+1]) {
					xidx = i+1;
					break;
				}
			}
		} // 조기상환일 관련 setting


		if(matuLoop == m_matu[curtN-1]-1) 
			intraTN = m_meshTN1;		
		else  
			intraTN = m_meshTN2;		
		dt = (1.0/YFactor)/intraTN;

		r = dailyIFrf[matuLoop]; // 금리: One Day implied forward rate 를 이용
		rd = dailyIFdc[matuLoop];

		if(m_voltype == 0) {
			volskew[0] = m_vol[0][0];
		}
		else {
			curttimeVol(volskew, matuLoop);				
		}

		// A*newV = oldV A 생성
		for(i=1; i<=SM-1; i++) {
			if(m_voltype == 2) { //local vol	
				Srv = S[i]*m_bval / m_volbval;			
				v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
			}
			else { // const vol or imp vol
				v = volskew[0];
			}

			if(SmeshType == 0) { // u mesh
				knockMAl[i] = (r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;
				knockMAd[i] = 1 + rd*dt + (v*v)*(i*i)*dt;
				knockMAu[i] = -(r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;				 
			}
			else { // A mesh
				knockMAl[i] = (r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i-1]*(h[i]+h[i-1]));
				knockMAd[i] = 1 + rd*dt + (v*v)*(S[i]*S[i])*dt/(h[i-1]*h[i]);
				knockMAu[i] = -(r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i]*(h[i]+h[i-1]));
			}

			if(m_kiFlag != 1) {			
				MAl[i] = knockMAl[i];			
				MAd[i] = knockMAd[i];			
				MAu[i] = knockMAu[i];
			}
		} //for(i) A만들기

		//경계 : Smax N 조건
		if(SmeshType == 0) {
			/*
			knockMAd[1] = knockMAd[1] + 2*knockMAl[1];
			knockMAu[1] = knockMAu[1] - knockMAl[1];
			*/

			knockMAl[SM-1] = knockMAl[SM-1] - knockMAu[SM-1];
			knockMAd[SM-1] = knockMAd[SM-1] + 2*knockMAu[SM-1];		 		
		}
		else {
			/*
			knockMAd[1] = knockMAd[1] + knockMAl[1]*(h[0]+h[1])/h[1];
			knockMAu[1] = knockMAu[1] - knockMAl[1]*h[0]/h[1];
			*/

			knockMAl[SM-1] = knockMAl[SM-1] - knockMAu[SM-1]*h[SM-1]/h[SM-2];
			knockMAd[SM-1] = knockMAd[SM-1] + knockMAu[SM-1]*(h[SM-1]+h[SM-2])/h[SM-2];			
		}

		// Smin D 조건
		knockMAd[0] = 1;			
		knockMAu[0] = 0;			
		knockMAl[0] = 0;

		if(m_kiFlag != 1) {
			/*
			MAd[1] = knockMAd[1];
			MAu[1] = knockMAu[1];
			*/
			MAl[SM-1] = knockMAl[SM-1];
			MAd[SM-1] = knockMAd[SM-1];

			
			MAd[startSidx] = 1;
			MAu[startSidx] = 0;
			MAl[startSidx] = 0;

			//LU decomposition
			for(i=startSidx+1; i<=SM-1; i++) {
				MAl[i] = MAl[i]/MAd[i-1];
				MAd[i] = MAd[i] - MAl[i]*MAu[i-1];
			}			
		}//if(m_kiFlag != 1)
		
		//LU decomposition
		//for(i=2; i<=SM-1; i++) { //Smin N 조건일 때
		for(i=1; i<=SM-1; i++) { //Smin D 조건일 때
			knockMAl[i] = knockMAl[i]/knockMAd[i-1];
			knockMAd[i] = knockMAd[i] - knockMAl[i]*knockMAu[i-1];
		}
		
		for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {
			DCValue = DCValue * 1/(1+rd/(YFactor*intraTN));
			CDlegDCValue = CDnewV[0] * 1/(1+rd/(YFactor*intraTN));

			memcpy(knockoldV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
			memcpy(CDoldV,CDnewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
			
			IsMidCompu = 0;

			if(intraLoop == 1 && mid_n > 0) { // 조기 상환일 중에 하나
				if(matuLoop == 0 && m_ctime > 15 && slevel < xlevel)   // 오늘이 조기상환일 이면서 종가 이후 상환 안된 때..
					IsMidCompu = 0; // 조기상환 안된 next 헤지를 위해서..														
				else
					IsMidCompu = 1;
			}

			if(IsMidCompu == 0) { //조기 상환일이 아닌 평시..형태로 풀기 (조기상환일에 상환 안된 경우도 해당)
				//변수명 절약? 하며서 A에 LU저장 newV가 tempV 역할..
				//끝점 선형 조건으로 계산 0~SM 중에.. 1부터 SM-1 까지 풀고 0과 SM은 경계조건으로 해결

				// L*tempV = oldV forward로 풀기
				 
				knocknewV[0] =  DCValue; //0.0; // D 0 조건
				//for(i=2; i<=SM-1; i++) // Smin N 조건
				for(i=1; i<=SM-1; i++) // Smin D 조건
					knocknewV[i] = knockoldV[i] - knockMAl[i]*knocknewV[i-1];

				//U*newV = tempV backward로 풀기
				knocknewV[SM-1] = knocknewV[SM-1]/knockMAd[SM-1];
				for(i=SM-2; i>=1; i--)
					knocknewV[i] = (knocknewV[i] - knockMAu[i]*knocknewV[i+1])/knockMAd[i];

				//경계 조건 적용
				if(SmeshType == 0) { // 균등 분할 			
					//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
					knocknewV[SM] = 2*knocknewV[SM-1] - knocknewV[SM-2];
				}
				else { // A mesh
					//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
					knocknewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * knocknewV[SM-1] - h[SM-1]/h[SM-2] * knocknewV[SM-2];
				}

				///////////////////////////////////cd
							 
				CDnewV[0] =  CDlegDCValue; //0.0; // D 0 조건
				//for(i=2; i<=SM-1; i++) // Smin N 조건
				for(i=1; i<=SM-1; i++) // Smin D 조건
					CDnewV[i] = CDoldV[i] - knockMAl[i]*CDnewV[i-1];

				//U*newV = tempV backward로 풀기
				CDnewV[SM-1] = CDnewV[SM-1]/knockMAd[SM-1];
				for(i=SM-2; i>=1; i--)
					CDnewV[i] = (CDnewV[i] - knockMAu[i]*CDnewV[i+1])/knockMAd[i];

				//경계 조건 적용
				if(SmeshType == 0) { // 균등 분할 			
					//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
					CDnewV[SM] = 2*CDnewV[SM-1] - CDnewV[SM-2];
				}
				else { // A mesh
					//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
					CDnewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * CDnewV[SM-1] - h[SM-1]/h[SM-2] * CDnewV[SM-2];
				}
				//////////////////////////////////////////////////////////////cd

				// no KI 영역 풀기 
				if(m_kiFlag != 1) {
					memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
							 
					//변수명 절약? 하며서 A에 LU저장 newV가 tempV 역할..
					//끝점 선형 조건으로 계산 0~SM 중에.. 1부터 SM-1 까지 풀고 0과 SM은 경계조건으로 해결
					//startS, endS는 확정 조건

					oldV[startSidx] = knocknewV[startSidx];

					// L*tempV = oldV forward로 풀기
					newV[startSidx] = oldV[startSidx]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
					for(i=startSidx+1; i<=SM-1; i++)
						newV[i] = oldV[i] - MAl[i]*newV[i-1];

					//U*newV = tempV backward로 풀기
					newV[SM-1] = newV[SM-1]/MAd[SM-1];
					for(i=SM-2; i>=startSidx; i--)
						newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

					//경계 조건 적용
					if(SmeshType == 0) { // 균등 분할 								
						newV[SM] = 2*newV[SM-1] - newV[SM-2];
					}
					else { // A mesh					
						newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
					}

					for(i=0; i<startSidx; i++) 
						newV[i] = knocknewV[i];
				}// if(m_kiFlag != 1)

			} //if(IsMidCompu == 0)
			else { // 조기상환일 계산

				// L*tempV = oldV forward로 풀기
				//knocknewV[1] = knockoldV[1]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
							 
				knocknewV[0] =  DCValue; //0.0; // D 0 조건
				//for(i=2; i<=SM-1; i++) // Smin N 조건
				for(i=1; i<endSidx; i++) // Smin D 조건
					knocknewV[i] = knockoldV[i] - knockMAl[i]*knocknewV[i-1];

				//U*newV = tempV backward로 풀기
				for(i=endSidx; i<=SM; i++)				
					knocknewV[i] =  AMOUNT * (m_Cpn[curtN-1]);// (1+m_Cpn[curtN-1]); // D 행사 조건으로 처리

				for(i=endSidx-1; i>=0; i--)
					knocknewV[i] = (knocknewV[i] - knockMAu[i]*knocknewV[i+1])/knockMAd[i];

				/*
				//경계 조건 적용
				if(SmeshType == 0) { // 균등 분할 			
					//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
					//knocknewV[SM] = 2*knocknewV[SM-1] - knocknewV[SM-2];
				}
				else { // A mesh
					//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
					//knocknewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * knocknewV[SM-1] - h[SM-1]/h[SM-2] * knocknewV[SM-2];
				}
				*/
				
				if(m_spreadXRate[curtN-1] > 0){// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
					double spreadvalue;

					for(i=0; i<=xidx; i++) {
						//spreadvalue = AMOUNT * (((1+m_Cpn[curtN-1])-0)/(xlevel - spreadXSlevel) * (S[i]-spreadXSlevel) + 0);
						spreadvalue = AMOUNT * (m_spreadXRate[curtN-1] * (S[i]-S[xidx]) + 1+m_Cpn[curtN-1]);
						//knocknewV[i] = MAX(knocknewV[i],spreadvalue);
						knocknewV[i] = MAX(knocknewV[i],spreadvalue-AMOUNT);
					}

					/*
					//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						maxindexR = SM;
						maxvalueR = knocknewV[SM];
						for(i=SM-1; i>=xidx; i--) {
							if(maxvalueR<knocknewV[i]) {
								maxindexR = i;
								maxvalueR = knocknewV[i];
							}
						}

						if(maxindexR > xidx) {
							tempindex = 2*xidx - spreadXSidx;

							for(i=maxindexR-1; i>tempindex; i--) {
								spreadvalue = (knocknewV[tempindex]-maxvalueR)/(S[tempindex]-S[maxindexR])*(S[i]-S[maxindexR]) + maxvalueR;
								knocknewV[i] = MAX(knocknewV[i],spreadvalue);
							}
						}
						*/

				
				}	 
	
	 
				////////////////////////////////////////cd
				
				// L*tempV = oldV forward로 풀기
				//knocknewV[1] = knockoldV[1]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
							 
				CDnewV[0] =  CDlegDCValue; //0.0; // D 0 조건
				//for(i=2; i<=SM-1; i++) // Smin N 조건
				for(i=1; i<endSidx; i++) // Smin D 조건
					CDnewV[i] = CDoldV[i] - knockMAl[i]*CDnewV[i-1];

				//U*newV = tempV backward로 풀기
				for(i=endSidx; i<=SM; i++)				
					CDnewV[i] =  AMOUNT * (CpnCDleg[curtN-1]);  //(1+m_Cpn[curtN-1]); // D 행사 조건으로 처리

				for(i=endSidx-1; i>=0; i--)
					CDnewV[i] = (CDnewV[i] - knockMAu[i]*CDnewV[i+1])/knockMAd[i];

				/*
				//경계 조건 적용
				if(SmeshType == 0) { // 균등 분할 			
					//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
					//knocknewV[SM] = 2*knocknewV[SM-1] - knocknewV[SM-2];
				}
				else { // A mesh
					//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
					//knocknewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * knocknewV[SM-1] - h[SM-1]/h[SM-2] * knocknewV[SM-2];
				}
				*/
				
				if(m_spreadXRate[curtN-1] > 0){// && CpnCDleg[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
					double spreadvalue;

					for(i=0; i<=xidx; i++) {
						//spreadvalue = AMOUNT * (((1+m_Cpn[curtN-1])-0)/(xlevel - spreadXSlevel) * (S[i]-spreadXSlevel) + 0);
						spreadvalue = AMOUNT * (m_spreadXRate[curtN-1] * (S[i]-S[xidx]) + 1+CpnCDleg[curtN-1]);
						//knocknewV[i] = MAX(knocknewV[i],spreadvalue);	
						CDnewV[i] = MAX(CDnewV[i],spreadvalue-AMOUNT);
					}
			 /*
	//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						maxindexR = SM;
						maxvalueR = CDnewV[SM];
						for(i=SM-1; i>=xidx; i--) {
							if(maxvalueR<CDnewV[i]) {
								maxindexR = i;
								maxvalueR = CDnewV[i];
							}
						}

						if(maxindexR > xidx) {
							tempindex = 2*xidx - spreadXSidx;

							for(i=maxindexR-1; i>tempindex; i--) {
								spreadvalue = (CDnewV[tempindex]-maxvalueR)/(S[tempindex]-S[maxindexR])*(S[i]-S[maxindexR]) + maxvalueR;
								CDnewV[i] = MAX(CDnewV[i],spreadvalue);
							}
						}
						*/
				
				}	 
	
	 
				//////////////////////////////////////////////////////////////////cd

				// no KI 영역 풀기 
				if(m_kiFlag != 1) {
					memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
							 
					//변수명 절약? 하며서 A에 LU저장 newV가 tempV 역할..
					//끝점 선형 조건으로 계산 0~SM 중에.. 1부터 SM-1 까지 풀고 0과 SM은 경계조건으로 해결
					//startS, endS는 확정 조건

					oldV[startSidx] = knocknewV[startSidx];

					// L*tempV = oldV forward로 풀기
					newV[startSidx] = oldV[startSidx]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
					for(i=startSidx+1; i<=SM-1; i++)
						newV[i] = oldV[i] - MAl[i]*newV[i-1];

					//U*newV = tempV backward로 풀기
					for(i=endSidx; i<=SM; i++) 
						newV[i] = AMOUNT *  (m_Cpn[curtN-1]); //(1+m_Cpn[curtN-1]); // D 행사 조건으로 처리

					for(i=endSidx-1; i>=0; i--)
						newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

					/*
					//경계 조건 적용
					if(SmeshType == 0) { // 균등 분할 								
						newV[SM] = 2*newV[SM-1] - newV[SM-2];
					}
					else { // A mesh					
						newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
					}
					*/

					for(i=0; i<startSidx; i++) 
						newV[i] = knocknewV[i];
								
					if(m_spreadXRate[curtN-1] > 0){// && m_Cpn[curtN-1]>0) { // X spread는 만기와 조기상환일에 가능
						double spreadvalue;

						for(i=0; i<=xidx; i++) {
							//spreadvalue = AMOUNT * (((1+m_Cpn[curtN-1])-0)/(xlevel - spreadXSlevel) * (S[i]-spreadXSlevel) + 0);							
							spreadvalue = AMOUNT * (m_spreadXRate[curtN-1] * (S[i]-S[xidx]) + 1+m_Cpn[curtN-1]);
							newV[i] = MAX(newV[i],spreadvalue-AMOUNT);
						}
					/*
	//kzR here
						int maxindexR;
						int tempindex;
						double maxvalueR;
						maxindexR = SM;
						maxvalueR = newV[SM];
						for(i=SM-1; i>=xidx; i--) {
							if(maxvalueR<newV[i]) {
								maxindexR = i;
								maxvalueR = newV[i];
							}
						}

						if(maxindexR > xidx) {
							tempindex = 2*xidx - spreadXSidx;

							for(i=maxindexR-1; i>tempindex; i--) {
								spreadvalue = (newV[tempindex]-maxvalueR)/(S[tempindex]-S[maxindexR])*(S[i]-S[maxindexR]) + maxvalueR;
								newV[i] = MAX(newV[i],spreadvalue);
							}
						}
						*/

					
					}	 

				}// if(m_kiFlag != 1)

			} //if (IsMidCompu == 0) else				 

		}//for(intraLoop)

		
			for(j=0; j<m_paymatuN; j++) {
				if(matuLoop == m_paymatu[j]) {
					payxlevel = m_payxval[j]/m_bval;
					for(i=0; i<=SM; i++) {
						if(S[i] >= payxlevel) {						
							knocknewV[i] += AMOUNT*(temppayCpn[j]);
							newV[i] += AMOUNT * (temppayCpn[j]);
						}
					}
					break;
				}
			}

		
		//cd payment Day!!
		cdpayvalue = 0.0;
		for(i=0; i<m_CDN; i++) {
			if(matuLoop == m_PayBDay[i]) {
				cdpayvalue = FlegCpn[i];
				break;
			}
		}

		if(cdpayvalue > 0.0) {			
			for(i=0; i<=SM; i++) 				
					CDnewV[i] += cdpayvalue;						
		}	




		//이산배당 적용
		disc_div = 0.0;
		for(i=0; i<m_divN; i++) {
			if(matuLoop == (m_divBDay[i]-1)) {
				disc_div = m_divCash[i];
				break;
			}
		}

		if(disc_div > 0.0) {
			double divSlevel[2];
			double divV[2];
			memcpy(knockoldV,knocknewV,(size_t)((SM+1)*sizeof(double))); 

			divSlevel[0] = (m_divbval*m_divApy - disc_div)/m_bval;
			divSlevel[1] = (m_divbval*m_divApy + disc_div)/m_bval;

			divV[0] = interp1(divSlevel[0],S,knockoldV,SM+1,SmeshType,1);
			divV[1] = interp1((divSlevel[1]+divSlevel[0])/2.0,S,knockoldV,SM+1,SmeshType,1);

			for(i=0; i<=SM; i++) {
				if(S[i] < divSlevel[0]) {
					knocknewV[i] = knockoldV[i];
				}
				else {
					if(S[i] > divSlevel[1]) {
						knocknewV[i] = interp1(S[i] - disc_div/m_bval, S,knockoldV,SM+1,SmeshType,1);
					}
					else {
						knocknewV[i] = interp1(S[i],divSlevel,divV,2,1,0); //2개이고 내부 점일 땐 1,0) 상관없음 
					}
				}
			}//for(i)

						
			memcpy(CDoldV,CDnewV,(size_t)((SM+1)*sizeof(double))); 

			divSlevel[0] = (m_divbval*m_divApy - disc_div)/m_bval;
			divSlevel[1] = (m_divbval*m_divApy + disc_div)/m_bval;

			divV[0] = interp1(divSlevel[0],S,CDoldV,SM+1,SmeshType,1);
			divV[1] = interp1((divSlevel[1]+divSlevel[0])/2.0,S,CDoldV,SM+1,SmeshType,1);

			for(i=0; i<=SM; i++) {
				if(S[i] < divSlevel[0]) {
					CDnewV[i] = CDoldV[i];
				}
				else {
					if(S[i] > divSlevel[1]) {
						CDnewV[i] = interp1(S[i] - disc_div/m_bval, S,CDoldV,SM+1,SmeshType,1);
					}
					else {
						CDnewV[i] = interp1(S[i],divSlevel,divV,2,1,0); //2개이고 내부 점일 땐 1,0) 상관없음 
					}
				}
			}//for(i)



			if(m_kiFlag != 1) {			
				memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); 

				divSlevel[0] = (m_divbval*m_divApy - disc_div)/m_bval;
				divSlevel[1] = (m_divbval*m_divApy + disc_div)/m_bval;

				divV[0] = interp1(divSlevel[0],S,oldV,SM+1,SmeshType,1);
				divV[1] = interp1((divSlevel[1]+divSlevel[0])/2.0,S,oldV,SM+1,SmeshType,1);

				for(i=0; i<=SM; i++) {
					if(S[i] < divSlevel[0]) {
						newV[i] = oldV[i];
					}
					else {
						if(S[i] > divSlevel[1]) {
							newV[i] = interp1(S[i] - disc_div/m_bval, S,oldV,SM+1,SmeshType,1);
						}
						else {
							newV[i] = interp1(S[i],divSlevel,divV,2,1,0); //2개이고 내부 점일 땐 1,0) 상관없음 
						}
					}
				}//for(i)
			} //if(m_knockFlag !=1 )

		} //if(disc_div > 0.0)
		//이산배당 적용 end

		if(m_kiFlag == 1) { // KI hit 시
			memcpy(newV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
		}
	 
		if(m_curtgreek == 1 && matuLoop == 1) { // 기본 그릭	Theta를 위한 +1일 저장
			for(i=0; i<outN + 4; i++) {	
							
				sim_knockFlag = 0;
				if(SrTemp[i] <= m_kihval)
					sim_knockFlag = 1;
				
				if(evalq<=0) {
						sleveltemp = SrTemp[i]/m_bval;
					}
				else { // n영업일 평균 처리
					if(m_ctime>15) { // 장종료
						sleveltemp = SrTemp[i]/m_bval * (m_evalN[curtN-1]-evalq); //  = m_matu[curtN-1]
						for(evali=0; evali<evalq; evali++)
							sleveltemp = sleveltemp + m_psval[evali]/m_bval;
						sleveltemp = sleveltemp / m_evalN[curtN-1];
					}
					else { // 장중
						sleveltemp = SrTemp[i]/m_bval * (m_evalN[curtN-1]-evalq+1);
						for(evali=1; evali<evalq; evali++)
							sleveltemp = sleveltemp + m_psval[evali]/m_bval;
						sleveltemp = sleveltemp / m_evalN[curtN-1];
					}
				}
				
				if(sim_knockFlag == 1) {				
					SrV1Temp[i] = interp1(sleveltemp,S,knocknewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
					SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				}
				else {
					SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
					SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				}
			}		 
		}
		
 	   

		if(m_curtgreek == 1 && mid_n > 0 && matuLoop == 0) { // 오늘이 조기 상환일 중에 하나
			if(m_ctime > 15 && slevel >= xlevel) { // 오늘이 조기상환일 이면서 종가 이후 상환 된 때..
				for(i=0; i<outN + 4; i++) {										
					sim_knockFlag = 0;
					if(SrTemp[i] <= m_kihval)
						sim_knockFlag = 1;

					sleveltemp = SrTemp[i]/m_bval;
					for(evali=1; evali<evalq; evali++)
						sleveltemp = sleveltemp + m_psval[evali]/m_bval;
					sleveltemp = sleveltemp / m_evalN[curtN-1];

					if(sim_knockFlag == 1) {					
						SrV1Temp[i] = interp1(sleveltemp,S,knocknewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
						SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
					}
					else {
						SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
						SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
					}
				}	
			}
			if(m_ctime > 15 && slevel < xlevel) { //조기상환일에 조기 상환 안된 경우
				expiryq = 0;
				slevel = m_sval/m_bval;
				
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
		logsgnuplot.open("c:/logland/debug_sdcdagnufdprice.txt",ios::out);	 
		for(i=0; i<=SM; i++) {
			logsgnuplot << S[i]*m_bval <<"     "<< newV[i] - CDnewV[i] << endl;
		}		 
		logsgnuplot.close();

		logsgnuplot.open("c:/logland/debug_avesdcdagnufdprice.txt",ios::out);	 
		double stemp, valuetemp;
		for(i=0; i<=SM; i++) {						
			sim_knockFlag = 0;
			if(S[i]*m_bval <= m_kihval)
				sim_knockFlag = 1;

			if(evalq <=1)
				stemp = S[i]*m_bval;
			else {
				stemp = S[i]*m_bval * (m_evalN[curtN-1]-(evalq-1));			
				for(evali=1;evali<evalq; evali++)				
					stemp = stemp + m_psval[evali];
				stemp = stemp / m_evalN[curtN-1];
			}

			if(sim_knockFlag == 1)	{			
				valuetemp = interp1(stemp/m_bval,S,knocknewV,SM+1,SmeshType,1);
				valuetemp -= interp1(stemp/m_bval,S,CDnewV,SM+1,SmeshType,1);
			}
			else {
				valuetemp = interp1(stemp/m_bval,S,newV,SM+1,SmeshType,1);
				valuetemp -= interp1(stemp/m_bval,S,CDnewV,SM+1,SmeshType,1);
			}

			logsgnuplot << S[i]*m_bval <<"     "<< valuetemp << endl;
		}		 
		logsgnuplot.close();
	}
	///////////////////////////////////////////--------------Debug	Mode e	
	 
	 
	if(m_curtgreek == 1 || m_curtgreek == 2) {
		for(i=0; i<outN + 4; i++) {			
			sim_knockFlag = 0;
			if(SrTemp[i] <= m_kihval)
				sim_knockFlag = 1;

			if(evalq <= 1)
				sleveltemp = SrTemp[i]/m_bval;
			else {
				sleveltemp = SrTemp[i]/m_bval * (m_evalN[curtN-1] - (evalq-1));
				for(evali=1; evali<evalq; evali++)
					sleveltemp = sleveltemp + m_psval[evali]/m_bval;
				sleveltemp = sleveltemp/m_evalN[curtN-1];
			}
			
			if(sim_knockFlag == 1)	{		
				SrVTemp[i] = interp1(sleveltemp,S,knocknewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				SrVTemp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
			}
			else {
				SrVTemp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
				SrVTemp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
			}
		}
		
		Baseidx = SoutRange + 2; // m_sval에 해당하는 idx 위 아래 2개씩 더 계산해서.. +2 필요
	}

	if(m_curtgreek == 1) { // 기본 그릭  
		value = interp1(slevel,S,newV,SM+1,SmeshType,1);
		value -= interp1(slevel,S,CDnewV,SM+1,SmeshType,1);
/*  
		greeks[Priceidx] = SrVTemp[Baseidx]; // [0] = price
		greeks[UpDeltaCidx] = (SrVTemp[Baseidx+1] - SrVTemp[Baseidx]) / 0.01; // [1] = up delta cash = delta * s
		greeks[DownDeltaCidx] = (SrVTemp[Baseidx] - SrVTemp[Baseidx-1]) / 0.01; // [2] = down delta cash
		greeks[UpGammaCidx] = (SrVTemp[Baseidx+2] - 2*SrVTemp[Baseidx+1] + SrVTemp[Baseidx]) / 0.01; // [3] = up gamma cash = gamma * s^2 * 0.01
		greeks[DownGammaCidx] = (SrVTemp[Baseidx] - 2*SrVTemp[Baseidx-1] + SrVTemp[Baseidx-2]) / 0.01; // [4] = down gamma cash
		greeks[Thetaidx] = (SrV1Temp[Baseidx] - SrVTemp[Baseidx]); // [5] = 1 day theta
		greeks[CharmCidx] = (SrV1Temp[Baseidx+1]-SrV1Temp[Baseidx-1]-SrVTemp[Baseidx+1]+SrVTemp[Baseidx-1])/(2*0.01); // [6] = 1 day charm cash
*/		
		greeks[Priceidx] = SrVTemp[Baseidx]; // [0] = price
		//greeks[UpDeltaCidx] = (SrVTemp[Baseidx+1] - SrVTemp[Baseidx]) / 0.01; // [1] = up delta cash = delta * s
		greeks[UpDeltaCidx] = (SrVTemp[Baseidx+1] - SrVTemp[Baseidx]); // [1] = up delta cash = delta * 0.01s  //2012-05-23 수정 1% delta cash로 수정
		//greeks[DownDeltaCidx] = (SrVTemp[Baseidx] - SrVTemp[Baseidx-1]) / 0.01; // [2] = down delta cash
		greeks[DownDeltaCidx] = (SrVTemp[Baseidx] - SrVTemp[Baseidx-1]); // [2] = down delta cash   // 2012-05-23 수정
		//greeks[UpGammaCidx] = (SrVTemp[Baseidx+2] - 2*SrVTemp[Baseidx+1] + SrVTemp[Baseidx]) / 0.01; // [3] = up gamma cash = gamma * s^2 * 0.01
		greeks[UpGammaCidx] = (SrVTemp[Baseidx+2] - 2*SrVTemp[Baseidx+1] + SrVTemp[Baseidx]); // [3] = up gamma cash = gamma * (0.01S)^2 // 2012-05-23 수정
		//greeks[DownGammaCidx] = (SrVTemp[Baseidx] - 2*SrVTemp[Baseidx-1] + SrVTemp[Baseidx-2]) / 0.01; // [4] = down gamma cash
		greeks[DownGammaCidx] = (SrVTemp[Baseidx] - 2*SrVTemp[Baseidx-1] + SrVTemp[Baseidx-2]); // [4] = down gamma cash // 2012-05-23 수정
		greeks[Thetaidx] = (SrV1Temp[Baseidx] - SrVTemp[Baseidx]); // [5] = 1 day theta
		//greeks[CharmCidx] = (SrV1Temp[Baseidx+1]-SrV1Temp[Baseidx-1]-SrVTemp[Baseidx+1]+SrVTemp[Baseidx-1])/(2*0.01); // [6] = 1 day charm cash
		greeks[CharmCidx] = (SrV1Temp[Baseidx+1]-SrV1Temp[Baseidx-1]-SrVTemp[Baseidx+1]+SrVTemp[Baseidx-1])/2; // [6] = 1 day charm cash //2012-05-23 1 day 1% charm cash 수정
		

		
		//델타 완화 2012-07-12
		greeks[UpDeltaCidx] = MAX(MIN(greeks[UpDeltaCidx],maxhgdelta),minhgdelta); //deltaadd
		greeks[DownDeltaCidx] = MAX(MIN(greeks[DownDeltaCidx],maxhgdelta),minhgdelta); //deltaadd
 
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
				batchsf << SrTemp[i] << " " << SrVTemp[i] << endl;
			}
			batchsf.close();



			strcpy_s(fname,kCodeN,m_kCode);
			strcat_s(fname,kCodeN,"_day+.txt");

			batchsf.open(fname,ios::out);
			
			batchsf << fixed;
			batchsf.precision(15);
			for(i=0; i<outN + 4; i++) {
				batchsf << SrTemp[i] << " " << SrV1Temp[i] << endl;
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
					volskew[0] = m_vol[0][0];
				}
				else {
					curttimeVol(volskew, 0);				
				}

				// A*newV = oldV A 생성
				for(i=1; i<=SM-1; i++) {
					if(m_voltype == 2) { //local vol	
						Srv = S[i]*m_bval / m_volbval;			
						v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
					}
					else { // const vol or imp vol
						v = volskew[0];
					}

					if(SmeshType == 0) { // u mesh
						knockMAl[i] = (r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;
						knockMAd[i] = 1 + rd*dt + (v*v)*(i*i)*dt;
						knockMAu[i] = -(r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;				 
					}
					else { // A mesh
						knockMAl[i] = (r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i-1]*(h[i]+h[i-1]));
						knockMAd[i] = 1 + rd*dt + (v*v)*(S[i]*S[i])*dt/(h[i-1]*h[i]);
						knockMAu[i] = -(r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i]*(h[i]+h[i-1]));
					}

					if(m_kiFlag != 1) {			
						MAl[i] = knockMAl[i];			
						MAd[i] = knockMAd[i];			
						MAu[i] = knockMAu[i];
					}
				} //for(i) A만들기

				//경계 : Smax N 조건
				if(SmeshType == 0) {
					/*
					knockMAd[1] = knockMAd[1] + 2*knockMAl[1];
					knockMAu[1] = knockMAu[1] - knockMAl[1];
					*/

					knockMAl[SM-1] = knockMAl[SM-1] - knockMAu[SM-1];
					knockMAd[SM-1] = knockMAd[SM-1] + 2*knockMAu[SM-1];		 		
				}
				else {
					/*
					knockMAd[1] = knockMAd[1] + knockMAl[1]*(h[0]+h[1])/h[1];
					knockMAu[1] = knockMAu[1] - knockMAl[1]*h[0]/h[1];
					*/

					knockMAl[SM-1] = knockMAl[SM-1] - knockMAu[SM-1]*h[SM-1]/h[SM-2];
					knockMAd[SM-1] = knockMAd[SM-1] + knockMAu[SM-1]*(h[SM-1]+h[SM-2])/h[SM-2];			
				}

				// Smin D 조건
				knockMAd[0] = 1;			
				knockMAu[0] = 0;			
				knockMAl[0] = 0;

				if(m_kiFlag != 1) {
					/*
					MAd[1] = knockMAd[1];
					MAu[1] = knockMAu[1];
					*/
					MAl[SM-1] = knockMAl[SM-1];
					MAd[SM-1] = knockMAd[SM-1];

			
					MAd[startSidx] = 1;
					MAu[startSidx] = 0;
					MAl[startSidx] = 0;

					//LU decomposition
					for(i=startSidx+1; i<=SM-1; i++) {
						MAl[i] = MAl[i]/MAd[i-1];
						MAd[i] = MAd[i] - MAl[i]*MAu[i-1];
					}			
				}//if(m_kiFlag != 1)
		
				//LU decomposition
				//for(i=2; i<=SM-1; i++) { //Smin N 조건일 때
				for(i=1; i<=SM-1; i++) { //Smin D 조건일 때
					knockMAl[i] = knockMAl[i]/knockMAd[i-1];
					knockMAd[i] = knockMAd[i] - knockMAl[i]*knockMAu[i-1];
				}
				///////////////

  
		
				for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 확인		
					DCValue = DCValue * 1/(1+rd/(YFactor*intraTN));
					CDlegDCValue = CDnewV[0] * 1/(1+rd/(YFactor*intraTN));
								 
					memcpy(knockoldV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
					memcpy(CDoldV,CDnewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..

					// L*tempV = oldV forward로 풀기
					 
					knocknewV[0] =  DCValue; //0.0;// D 0 조건
					//for(i=2; i<=SM-1; i++) // Smin N 조건
					for(i=1; i<=SM-1; i++) // Smin D 조건
						knocknewV[i] = knockoldV[i] - knockMAl[i]*knocknewV[i-1];

					//U*newV = tempV backward로 풀기
					knocknewV[SM-1] = knocknewV[SM-1]/knockMAd[SM-1];
					for(i=SM-2; i>=1; i--)
						knocknewV[i] = (knocknewV[i] - knockMAu[i]*knocknewV[i+1])/knockMAd[i];

					//경계 조건 적용
					if(SmeshType == 0) { // 균등 분할 			
						//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
						knocknewV[SM] = 2*knocknewV[SM-1] - knocknewV[SM-2];
					}
					else { // A mesh
						//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
						knocknewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * knocknewV[SM-1] - h[SM-1]/h[SM-2] * knocknewV[SM-2];
					}

					///////////////////////////////////cd
							 
					CDnewV[0] =  CDlegDCValue; //0.0;// D 0 조건
					//for(i=2; i<=SM-1; i++) // Smin N 조건
					for(i=1; i<=SM-1; i++) // Smin D 조건
						CDnewV[i] = CDoldV[i] - knockMAl[i]*CDnewV[i-1];

					//U*newV = tempV backward로 풀기
					CDnewV[SM-1] = CDnewV[SM-1]/knockMAd[SM-1];
					for(i=SM-2; i>=1; i--)
						CDnewV[i] = (CDnewV[i] - knockMAu[i]*CDnewV[i+1])/knockMAd[i];

					//경계 조건 적용
					if(SmeshType == 0) { // 균등 분할 			
						//knocknewV[0] = 2*knocknewV[1] - knocknewV[2];
						CDnewV[SM] = 2*CDnewV[SM-1] - CDnewV[SM-2];
					}
					else { // A mesh
						//knocknewV[0] = (h[0]+h[1])/h[1] * knocknewV[1] - h[0]/h[1] * knocknewV[2];
						CDnewV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * CDnewV[SM-1] - h[SM-1]/h[SM-2] * CDnewV[SM-2];
					}

					///////////////////////////////////////////////////////////////cd
					// no KI 영역 풀기 
					if(m_kiFlag != 1) {
						memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
							 
						//변수명 절약? 하며서 A에 LU저장 newV가 tempV 역할..
						//끝점 선형 조건으로 계산 0~SM 중에.. 1부터 SM-1 까지 풀고 0과 SM은 경계조건으로 해결
						//startS, endS는 확정 조건

						oldV[startSidx] = knocknewV[startSidx];

						// L*tempV = oldV forward로 풀기
						newV[startSidx] = oldV[startSidx]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
						for(i=startSidx+1; i<=SM-1; i++)
							newV[i] = oldV[i] - MAl[i]*newV[i-1];

						//U*newV = tempV backward로 풀기
						newV[SM-1] = newV[SM-1]/MAd[SM-1];
						for(i=SM-2; i>=startSidx; i--)
							newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

						//경계 조건 적용
						if(SmeshType == 0) { // 균등 분할 								
							newV[SM] = 2*newV[SM-1] - newV[SM-2];
						}
						else { // A mesh					
							newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
						}

						for(i=0; i<startSidx; i++) 
							newV[i] = knocknewV[i];
					}// if(m_kiFlag != 1)	  
					else	
						memcpy(newV,knocknewV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
						
					if(m_curtgreek == 1) { // 만기에
						for(i=0; i<outN + 4; i++) {										
							sim_knockFlag = 0;			
							if(SrTemp[i] <= m_kihval)				
								sim_knockFlag = 1;

							if(evalq<=1) 							
								sleveltemp = SrTemp[i]/m_bval;
							else {
								sleveltemp = SrTemp[i]/m_bval * (m_evalN[curtN-1]-(evalq-1));
								for(evali=1; evali<evalq; evali++)
									sleveltemp = sleveltemp + m_psval[evali]/m_bval;
								sleveltemp = sleveltemp/m_evalN[curtN-1];
							}

							if(sim_knockFlag == 1) { 							
								SrV1Temp[i] = interp1(sleveltemp,S,knocknewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
								SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
							}
							else {
								SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
								SrV1Temp[i] -= interp1(sleveltemp,S,CDnewV,SM+1,SmeshType,1); // (0: flat, 1: 선형)
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
						batchsf << SrTemp[i] << " " << SrV1Temp[i] << endl;
					}
					batchsf.close();

				}//for(intraLoop)

				delete tempstring;
			} //if(m_matu == 0)
			
	 


/////

			delete fname;
			
		}// if(m_batchFlag)

		free(SrTemp);
		free(SrVTemp);
		free(SrV1Temp);
 
	}

	if(m_curtgreek == 2) { // vega 
		value = interp1(slevel,S,newV,SM+1,SmeshType,1);
		value -= interp1(slevel,S,CDnewV,SM+1,SmeshType,1);
		for(i=0; i<3; i++) {
			greeks[i] = SrVTemp[Baseidx + (i-1)]; // greeks[0]: s-, greeks[1]: s0, greeks[2]: s+
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
				strcat_s(fname,kCodeN,"_sigma+_");
			}
			else {
				strcat_s(fname,kCodeN,"_sigma-_");
			}			
			_itoa_s(m_curtvegaidx,tempstring,kCodeN,10);			
			strcat_s(fname,kCodeN,tempstring);			
			strcat_s(fname,kCodeN,".txt");


		 
			batchsf.open(fname,ios::out);	  	 
			
			batchsf << fixed;
			batchsf.precision(15);
	 
			for(i=0; i<outN + 4; i++) {
				batchsf << SrTemp[i] << " " << SrVTemp[i] << endl;
			}
			batchsf.close();

 			delete tempstring;
			delete fname;
			
		}// if(m_batchFlag)



		free(SrTemp);
		free(SrVTemp);
		free(SrV1Temp);
 
	}

	if(m_curtgreek == 3) { // rho

		value = interp1(slevel,S,newV,SM+1,SmeshType,1); // value - m_Ogreeks[0] 을 이용하기 위하여
		value -= interp1(slevel,S,CDnewV,SM+1,SmeshType,1); // value - m_Ogreeks[0] 을 이용하기 위하여
	}
	 

 
	free(CDoldV);
	free(CDnewV);
	
	free(S);
	free(h);
	free(oldV);
	free(newV);
	free(MAu);
	free(MAd);
	free(MAl);
		
	free(knockoldV);
	free(knocknewV);
	free(knockMAu);
	free(knockMAd);
	free(knockMAl);
 
	free(volskew);

	
	delete fdObCDrate;
	delete FlegCpn;
	delete CpnCDleg;
			
	delete temppayCpn;

	return value;
		

} // fd_price();



void StepDownCPCDA::curttimeVol(double *vol, int day)
{
	int i;
	
	if(m_voltype == 1) { // term imp vol
		double v1, v2;
		int t1, t2;
		double *termvol;

		termvol = new double[m_volTN];
		for(i=0; i<m_volTN; i++) 
			termvol[i] = m_vol[0][i];
		
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
			vol[i] = m_vol[i][findidx];
	}
	
}
