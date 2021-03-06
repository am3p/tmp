#include "TomCDA.h"
#include <iostream>
#include <fstream>

extern odysseyRandom *opr;

 

TomCDA::TomCDA(double sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,
	      int matuN, int *matu, double *prate, double ctime,	
		  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,
		  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,			   
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  long SimulNo)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	
			
	 
	
	
	int i, j;	 
	 
	//입력 데이터 내부 전역 변수로 저장   
	  
	m_matuN = matuN;
	m_matu = new int[m_matuN];
	m_prate = new double[m_matuN+1];

	m_STDayN = STDayN;
	m_STStartDay = new int[m_STDayN];
	m_IXT = new double[m_STDayN];
	m_IndexT = new double[m_STDayN];
	m_STEndDay = new int[m_STDayN];
	m_IXTheta = new double[m_STDayN];
		
	for(i=0; i<AsNUM1; i++) {
		m_sval[i] = sval; 
	}
	  
	m_stom = stom;
	m_btom = btom;

	  
 

	for(i=0; i<m_matuN; i++) {	 		 
		m_matu[i] = matu[i];	
		m_prate[i] = prate[i];
	}
	m_prate[m_matuN] = prate[m_matuN];

	for(i=0; i<m_STDayN; i++) {
		m_STStartDay[i] = STStartDay[i];
		m_IXT[i] = IXT[i];
		m_IndexT[i] = IndexT[i];
		m_STEndDay[i] = STEndDay[i];
		m_IXTheta[i] = IXTheta[i];
	}
  

	
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

	for(i=0; i<AsNUM1; i++) {
		m_divrate[i] = divrate[i];
		m_divbval[i] = divbval[i];
		m_divN[i] = divN[i];
	}

	m_MAXdivN = m_divN[0];
	for(i=1; i<AsNUM1; i++) {
		if(m_MAXdivN < m_divN[i])
			m_MAXdivN = m_divN[i];
	}

	m_divBDay = new int*[AsNUM1];
	m_divCash = new double*[AsNUM1];
	for(i=0; i<AsNUM1; i++) {
		m_divBDay[i] = new int[m_MAXdivN];
		m_divCash[i] = new double[m_MAXdivN];
	}

	int sumdivN;
	sumdivN = 0;
	for(i=0; i<AsNUM1; i++) {		
		for(j=0; j<m_divN[i]; j++) {
			m_divBDay[i][j] = divBDay[sumdivN + j];
			m_divCash[i][j] = divCash[sumdivN + j];
		}
		sumdivN = sumdivN + m_divN[i];
	}
	
	m_divApy = divApy;

	m_voltype = voltype;
	for(i=0; i<AsNUM1; i++) {
		m_volbval[i] = volbval[i];
	}
	 
	m_volTN = volTN;
	m_volSN = volSN;
	m_volBDay = new int[m_volTN];
	for(i=0; i<m_volTN; i++)
		m_volBDay[i] = volBDay[i];
	m_volSmness = new double[m_volSN];
	for(i=0; i<m_volSN; i++)
		m_volSmness[i] = volSmness[i];

	m_vol = new double**[AsNUM1];
	for(i=0; i<AsNUM1; i++) {
		m_vol[i] = new double*[m_volSN];
		for(j=0; j<m_volSN; j++) {
			m_vol[i][j] = new double[m_volTN];
		}
	}

	int k;

	for(k=0; k<AsNUM1; k++) {
		for(i=0; i<m_volSN; i++) {
			for(j=0; j<m_volTN; j++) {
				m_vol[k][i][j] = vol[m_volSN*m_volTN*k + m_volTN*i + j];
			}
		}
	}
	 
	 
	  

	
	//rnd = opr;
	//m_noiters = rnd->SimulNo;
	rnd = new odysseyRandom();
	rnd->SimulNo = SimulNo;
	if(rnd->randseq(0,SimulNo,MAXDAYS,NULL)!=0)
		exit(1);

	if(rnd->randnum(0,SimulNo,MAXAST,NULL,7)!=0)
		exit(1);
	 
	 m_noiters = rnd->SimulNo;
}

TomCDA::~TomCDA()
{	
	
	int i,k;	 
	//delete m_evalN;
	delete m_matu;
	delete m_prate;

	
	//cd
	delete m_ObBDay;
	delete m_PayBDay;
	delete m_FlegFactor;
	delete m_ObCDrate;




	delete m_STStartDay;
	delete m_IXT;
	delete m_IndexT;
	delete m_STEndDay;
	delete m_IXTheta;


	delete m_irateBDay;	
	delete m_iDCrate;  
	delete m_iIRSrate; 	 			
	delete m_volBDay;	 
	delete m_volSmness;	 	 
	 
	 
 

	//
	for(k=0; k<AsNUM1; k++) {
		for(i=0; i<m_volSN; i++) {
			delete m_vol[k][i];
		}	 
	}

	for(k=0; k<AsNUM1; k++) {
		delete m_vol[k];
 
		 
		delete m_divBDay[k];
		delete m_divCash[k];
	 
	}
	
	 
	
	delete m_vol;
	  
	delete m_iRFrate;
	delete m_divBDay;
	delete m_divCash;
	 
	rnd->randseq(2,m_noiters,MAXDAYS,NULL);
	rnd->randnum(2,m_noiters,MAXAST,NULL,7);
	delete rnd;
}

double TomCDA::getValue() { 
	
	double value;	 
	
	value = mc_price(); 
		
	m_Ogreeks[Priceidx] = value;
	 
  

	return value;

}


double TomCDA::getBatch(double tsval) {
 
	 
	double value;
	int i;

	double tempsval[AsNUM1];
 
	for(i=0; i<AsNUM1; i++) {
		tempsval[i] = m_sval[i];
		m_sval[i] = tsval;
	}
	
	value =  mc_price();

	for(i=0; i<AsNUM1; i++) {		 
		m_sval[i] = tempsval[i];
	}
	 
	return value;

}

 

double TomCDA::getDeltaGammaTheta() {
 
	double tempvalue[AsNUM1][2];
	double value;
	int i,k;

	double tempsval[AsNUM1];
 
	for(i=0; i<AsNUM1; i++) {
		tempsval[i] = m_sval[i];
		 
	}

	for(i=0; i<AsNUM1; i++) {
		m_sval[i] = tempsval[i]*(1+0.01);
		value = mc_price();
		tempvalue[i][0] = value;

		m_sval[i] = tempsval[i]*(1-0.01);
		value = mc_price();
		tempvalue[i][1] = value;

		m_sval[i] = tempsval[i];
	}
	
	m_Ogreeks[UpDeltaCidx] = tempvalue[0][0] - m_Ogreeks[Priceidx];
	m_Ogreeks[DownDeltaCidx] = m_Ogreeks[Priceidx] - tempvalue[0][1];
	 
	
	m_Ogreeks[UpGammaCidx] = tempvalue[0][0] - 2*m_Ogreeks[Priceidx] + tempvalue[0][1];
	m_Ogreeks[DownGammaCidx] = tempvalue[0][0] - 2*m_Ogreeks[Priceidx] + tempvalue[0][1];
	 

	 
	if(m_matu[m_matuN-1] <= 0) {		
		value = m_Ogreeks[Priceidx];
	}
	else {
	
		for(i=0; i<m_matuN; i++)		
			m_matu[i] = m_matu[i] - 1;

		for(i=0; i<m_STDayN; i++) {
			m_STStartDay[i] = m_STStartDay[i] - 1;
			m_STEndDay[i] = m_STEndDay[i] - 1;
			if(m_STStartDay[i] == 0) {
				m_IXT[i] = m_stom;
				m_IndexT[i] = m_sval[0];
			}
			if(m_STEndDay[i] == 0)
				m_IXTheta[i] = m_stom;
		}

		
			

		for(k=0; k<AsNUM1; k++)
			for(i=1;i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] - 1;
		 
			
		for(i=0; i<m_CDN; i++) { //cd
			m_ObBDay[i] = m_ObBDay[i] - 1;				
			m_PayBDay[i] = m_PayBDay[i] - 1;
		}

		value = mc_price();

		for(i=0; i<m_matuN; i++)
			m_matu[i] = m_matu[i] +1;

		for(i=0; i<m_STDayN; i++) {
			m_STStartDay[i] = m_STStartDay[i] +1;
			m_STEndDay[i] = m_STEndDay[i] +1;
			if(m_STStartDay[i] == 1) {
				m_IXT[i] = 0;
				m_IndexT[i] = 0;
			}
			if(m_STEndDay[i] == 1) {
				m_IXTheta[i] = 0;
			}
		}

				
		for(i=0; i<m_CDN; i++) {
			m_ObBDay[i] = m_ObBDay[i] + 1;				
			m_PayBDay[i] = m_PayBDay[i] + 1;
		}


		for(k=0; k<AsNUM1; k++)
			for(i=1; i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] + 1;
		 

	}

		 
	 
 
	m_Ogreeks[Thetaidx] = value - m_Ogreeks[Priceidx];

	return value;

}

double TomCDA::getVega(int vegaN, int *vegaidx) {
		 
	double value;
	int i,j,k,udi;

	double **tempv;
	
	int pertubationX[3];
	double pertubationY[3];		
	double curtpertubation;
	double Price[3][2];
	double VegaPertubationR;
	
	//m_curtgreek = 2;
		
	tempv = new double*[m_volSN]; // m_volSN은.. 1이상
	for(i=0; i<m_volSN; i++)
		tempv[i] = new double[m_volTN];

	// v1 시작
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[0][i][j]; // vol 저장
	

	m_Ogreeks[Vegaidx] = 0;
	 

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
						m_vol[0][k][j] = tempv[k][j] + curtpertubation;
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
							m_vol[0][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[0][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[0][k][j] = tempv[k][j] + curtpertubation;
						}
					} //else 3 if(i == vegaN-1)
				} //else 2 if(i == 0)
			}// else 1 if(m_voltype == 0)

			value = mc_price();
			
			for(j=0; j<3; j++) 			
				Price[j][udi] = value; //greeks[j]; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+
			
			
		} //for(vi)
		//
		/*
		m_Otvegas1[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas1[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas1[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/
		
		//2012-05-23 1% vanna cash / 1% zomma cash
		m_Otvegas[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
	 

		m_Ogreeks[Vegaidx] += m_Otvegas[i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[0][i][j] = tempv[i][j]; // vol 복원

	/////////////////////v1 end

	
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	 
	return m_Ogreeks[Vegaidx];

}


double TomCDA::getRho(int rhoN,int *rhoidx) {
	
	 
	double value;

	double *tempirate, *tempirate1;  
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		  

	 

	tempirate = new double[m_irateN];
	tempirate1 = new double[m_irateN];
	 
	 


	m_Ogreeks[Rhorfidx] = 0; //rf rho sum
 
	 
	for(i=0; i<MAXRHOS; i++)
		m_Otrfrhos[i] = 0;

	//rf1 rho 시작
	 

	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iRFrate[i];
		tempirate1[i] = m_iIRSrate[i];
	}
 
		 
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // 평행이동 rho
			for(j=0; j<m_irateN; j++) {
				m_iRFrate[j] = tempirate[j] + RhoPertubation;
				m_iIRSrate[j] = tempirate1[j] + RhoPertubation;
 
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
					m_iIRSrate[j] = tempirate1[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
	 
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
						m_iIRSrate[j] = tempirate1[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
	 
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
						m_iIRSrate[j] = tempirate1[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat
	 					
					}	

				}
			}
		}// if(rhoN == 0) else

		value = mc_price();

		m_Otrfrhos[i] = value - m_Ogreeks[Priceidx]; // 10bp term rf rho 
		m_Ogreeks[Rhorfidx] += m_Otrfrhos[i];
			 	 
	}
 

	for(i=0; i<m_irateN; i++) {
		m_iRFrate[i] = tempirate[i]; // RF rate 복원		
		m_iIRSrate[i] = tempirate1[i]; // RF rate 복원		
	}
 
		
	 
 

//////////////////////////////
	//dc rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iDCrate[i]; // DC rate 저장
	}
	 
	m_Ogreeks[Rhodcidx] = 0; // dc rho sum

	 
	for(i=0; i<MAXRHOS; i++)
		m_Otdcrhos[i] = 0;



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

		value = mc_price();

		m_Otdcrhos[i] = value - m_Ogreeks[Priceidx]; // 10bp term rf rho 
		m_Ogreeks[Rhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate 복원
	}
 	
//////////////////////////////
	 


	delete tempirate;
	delete tempirate1;
 

		
	return m_Ogreeks[Rhorfidx] + m_Ogreeks[Rhodcidx] ;  // 의미 없음 

}


 
 

double TomCDA::mc_price()
{
	long *rsequns;
	double **randoms;
	long nostep;

	double simsval[AsNUM1];
 
	double simprice;
	double price;
	double dT, sqdT;
		
	int i,j,k;

	dT = 1.0/YFactor;
	sqdT = sqrt(dT);
	double mudt[AsNUM1];
	double vldt[AsNUM1];
	//m_modeltype = 0;

	double disc_div[AsNUM1];
	double Srv ;
	double *volskew ;

	double r[AsNUM1], q[AsNUM1], v[AsNUM1];
 
	double *dailyIFrf, *dailyZdc;

	
	//cd
	double *fdObCDrate, *FlegCpn, *CpnCDleg;
 
	 
	

	//데이터 체크
	for(i=0; i<m_STDayN; i++) {
		if(m_STStartDay[i] == 0 && m_IXT[i] == 0)
			m_IXT[i] = m_stom;
		if(m_STStartDay[i] == 0 && m_IndexT[i] == 0)
			m_IndexT[i] = m_sval[0];
		if(m_STEndDay[i] == 0 && m_IXTheta[i] == 0)
			m_IXTheta[i] = m_stom;
	}

  

	nostep = MAX(m_matu[m_matuN-1],1);
	 

	rsequns = (long*)malloc((size_t)((nostep)*sizeof(long)));
	randoms =(double **)malloc((size_t)((AsNUM1)*sizeof(double *)));	
	if(!rsequns || !randoms) return -100;

	for ( i = 0; i < AsNUM1; i++ ) {
		randoms[i] = (double *)malloc((size_t)((m_noiters)*sizeof(double)));		 
		if(!randoms[i]) {
			for(j=0; j<i; j++)
				free(randoms[j]);
			return -200;
		}
	}

	if ( rnd->randnum(1,m_noiters,AsNUM1,randoms,0) != 0 ) {
		free(rsequns);
		for(i=0; i<AsNUM1; i++)
			free(randoms[i]);
		free(randoms);
		return(-600);
	}

	dailyIFrf = new double[nostep];		
	
	dailyZdc = new double[nostep+1];

	 
	for(i=0; i<nostep; i++) { //사용할 때 i=1:nostep i 대응 [i-1]로
		dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
	}
 

	for(i=0; i<=nostep; i++) { //사용할때 i 대응 [i]로 사용
		dailyZdc[i] = interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0);
	}
	 

	volskew  = new double[m_volSN];
	
		 
	fdObCDrate = new double[m_CDN];
	FlegCpn = new double[m_CDN];
	CpnCDleg = new double[m_matuN];

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

		if(i!= m_CDN-1 && m_PayBDay[i] == 0 ) //&& m_ctime > 15)
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


	 
	for(i=0; i<AsNUM1; i++)
		q[i] = m_divrate[i];

	int si,l,sk,disci;
	
	double cdprice;
	int cddayi;
	double simcdprice;
	 

	  
	price = 0; 
	cdprice = 0;
 
	double simstom;

	double *simIXT, *simIndexT, *simIXTheta;

	simIXT = new double[m_STDayN];
	simIndexT = new double[m_STDayN];
	simIXTheta = new double[m_STDayN];

	 
	
	for(si=0; si<m_noiters; si++) { //MC SIM		 
		simprice = 0;
		simcdprice = 0;
		if(rnd->randseq(1,si,nostep,rsequns) != 0) {
			for(i=0; i<AsNUM1; i++) {
				free(randoms[i]);					 			 
			}
			free(randoms);
			delete dailyIFrf;
			free(rsequns);
			delete dailyZdc;
			delete volskew; 
			
			delete fdObCDrate;
			delete FlegCpn;	
			delete CpnCDleg;

			return -900;
		}

		l = 0;
		for(sk=0; sk<m_matuN; sk++) {
			if(m_matu[sk] < 0) 
				continue;

			for(;l<=m_matu[sk]; l++) {
				if(l==0) {
					for(k=0; k<AsNUM1; k++)	{				
						simsval[k] = m_sval[k];						 
					}
					simstom = m_stom;				
					for(k=0; k<m_STDayN; k++) {
						simIXT[k] =m_IXT[k];
						simIndexT[k] = m_IndexT[k];
						simIXTheta[k] = m_IXTheta[k];
					}
				}
				else { //if l
					for(k=0; k<AsNUM1; k++) {
						r[k] = dailyIFrf[l-1];

						//vol
						if(m_voltype == 0) {
							volskew[0] = m_vol[k][0][0];							
						}
						else {
							curttimeVol(volskew, l, k);
						}

						if(m_voltype == 2) {
							Srv = simsval[k]/m_volbval[k];
							v[k] = interp1(Srv, m_volSmness,volskew,m_volSN,1,0);
						}
						else {
							v[k] = volskew[0];
						}
 
						  

						//disc_div					 
						for(disci=1; disci<m_divN[k]; disci++) {
							if(l==m_divBDay[k][disci] && m_divCash[k][disci]>0) {	
								double divS[2];
								double pole[]={0,1};
								double alpha;
					 
								disc_div[k] = m_divCash[k][disci];

								divS[0] = (m_divbval[k]*m_divApy-disc_div[k]);
								divS[1] = (m_divbval[k]*m_divApy+disc_div[k]);
								alpha = interp1(simsval[k],divS,pole,2,1,0);
								 
								simsval[k] = simsval[k]-alpha*disc_div[k]; 
								break;
							}
						}

						mudt[k] = (r[k]-q[k]-0.5*v[k]*v[k])*dT;
						vldt[k] = v[k]*sqdT;
					}

					 
					for(k=0; k<AsNUM1; k++) {							 
						simsval[k] = simsval[k] * exp(mudt[k] + vldt[k]*randoms[k][rsequns[l-1]]);						 
					}

					for(k=0; k<m_STDayN; k++) {
						if(l == m_STStartDay[k]) {
							if(k==0) {
								simIXT[k] = m_stom;								 
							}
							else {
								simIXT[k] = simIXTheta[k-1];								 
							}
							simIndexT[k] = simsval[0];
						}

						if(l == m_STEndDay[k]) {
							simIXTheta[k] = simIXT[k] * simsval[0]/simIndexT[k];
						}
					}

					
					if(l <= m_STStartDay[0]) {
						simstom = m_stom;
					}
					
					if(l > m_STEndDay[m_STDayN-1]) {
						simstom = simIXTheta[m_STDayN-1];
					}
										
					for(k=0; k<m_STDayN; k++) {
						if(m_STStartDay[k] < l && l < m_STEndDay[k])
							simstom = simIXT[k] * simsval[0]/simIndexT[k];
					}

					for(k=0; k<m_STDayN-1; k++) {
						if(m_STEndDay[k] <= l && l <= m_STStartDay[k+1])
							simstom = simIXTheta[k];
					}
 

				}//if l==0 else

				
				for(cddayi=0; cddayi<m_CDN-1; cddayi++) { // 만기는 따로 처리 하기 위해서
					if(l == m_PayBDay[cddayi]) {
						simcdprice = simcdprice + AMOUNT * (FlegCpn[cddayi]) * exp(-dailyZdc[l]*l/YFactor);
						break;
					}
				}

			}//for l
			
			simprice = simprice + AMOUNT * MAX(0, ((simstom-m_btom)/m_btom /(sk+1)-m_prate[m_matuN]) * m_prate[sk]) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);

			if(sk==m_matuN-1) {
				//simprice = simprice + AMOUNT * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);

				simcdprice = simcdprice + AMOUNT * CpnCDleg[sk] * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
			}
				
		} // for sk

		price = price + simprice/m_noiters;
		cdprice = cdprice + simcdprice/m_noiters;


	} // MC SIM
	 
 
	
		 
		
	//메모리 해제
	for(i=0; i<AsNUM1; i++) {
		free(randoms[i]);
		 
	}
	free(randoms);
	delete dailyIFrf;
	free(rsequns);
	delete dailyZdc;
	delete volskew; 
	 

	delete simIXT;
	delete simIndexT;
	delete simIXTheta;
		
	delete fdObCDrate;
	delete FlegCpn;
	delete CpnCDleg; 

	// 메모리 해제
 
		 
			
	return price-cdprice;
	 		



} //mc_price();



void TomCDA::curttimeVol(double *vol, int day, int Asidx)
{
	int i;
	
	if(m_voltype == 1) { // term imp vol
		double v1, v2;
		int t1, t2;
		double *termvol;

		termvol = new double[m_volTN];
		for(i=0; i<m_volTN; i++) 
			termvol[i] = m_vol[Asidx][0][i];
		
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
			vol[i] = m_vol[Asidx][i][findidx];
	}
	
}
 

  


/*
double TomCDA::curttimeFXCorr2(int day)
{
	int i;		
	int findidx;
		
	findidx = m_FXcorrN-1;
	for(i=0; i<m_FXcorrN; i++) {
		if(day <= m_FXcorrBDay[i]) {
			findidx = i;
			break;
		}
	}		 

	return m_FXcorr2[findidx];
	 
	
}
*/

//heresky




int TomCDA::findinteridx(double x, double *Dx, int DxN)
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

 
 