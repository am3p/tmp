#include "StabilityA.h"
#include <iostream>
#include <fstream>

extern odysseyRandom *opr;

 

StabilityA::StabilityA(int erFlag, double sval, double psval, double erPerf, double xval, double multiplier,
		  double *CpnFactor, double *Cpn,
	      int matuN, int *matu, double ctime,		   
		  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  double SLimit, long SimulNo)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	
			
	 
	
	
	int i, j;	 
	 
	//입력 데이터 내부 전역 변수로 저장   
	  
	m_matuN = matuN;
	m_matu = new int[m_matuN];

	
	for(i=0; i<AsNUM1; i++) {
		m_sval[i] = sval;
		m_psval[i] = psval;
		m_xval[i] = xval;
	}
	 
	
	m_erFlag = erFlag;
	m_erPerf = erPerf;
	m_multiplier = multiplier;
	 

	m_CpnFactor = new double[m_matuN];
	m_Cpn = new double[m_matuN];
 

	for(i=0; i<m_matuN; i++) {	 
		m_CpnFactor[i] = CpnFactor[i];
		m_Cpn[i] = Cpn[i];
		m_matu[i] = matu[i];		 		
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

	m_SLimit = SLimit;
	 

	 
	 

	// ki 체크
	for(i=0; i<AsNUM1; i++) {
		if(m_sval[i]/m_psval[i] <= m_xval[i]) {
			m_erFlag = 1; // 하나라도 치면..
			if(m_erPerf < 0)
				m_erPerf = m_sval[i]/m_psval[i];
			break;
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

StabilityA::~StabilityA()
{	
	
	int i,k;	 
	//delete m_evalN;
	delete m_matu;
	delete m_Cpn;
	delete m_CpnFactor;

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

double StabilityA::getValue() { 
	
	double value;	 
	
	value = mc_price(); 
		
	m_Ogreeks[Priceidx] = value;
	 
  

	return value;

}



 

double StabilityA::getDeltaGammaTheta() {
 
	double tempvalue[AsNUM1][2];
	double value;
	int i,k;

	double tempsval[AsNUM1];
	double temppsval[AsNUM1];

	for(i=0; i<AsNUM1; i++) {
		tempsval[i] = m_sval[i];
		temppsval[i] = m_psval[i];
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
		for(k=0; k<AsNUM1; k++)
			for(i=1;i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] - 1;
		for(k=0; k<AsNUM1; k++) {
			m_psval[k] = tempsval[k];
		}

		value = mc_price();

		for(i=0; i<m_matuN; i++)
			m_matu[i] = m_matu[i] +1;
		for(k=0; k<AsNUM1; k++)
			for(i=1; i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] + 1;
		for(k=0; k<AsNUM1; k++) {
			m_psval[k] = temppsval[k];
		}

	}

		 

	if(m_erFlag == 1) {
		for(k=0; k<m_matuN-1; k++) {
			if(m_matu[k] == 0) {
				value = m_Ogreeks[Priceidx];
			}
		}
	}
 
	m_Ogreeks[Thetaidx] = value - m_Ogreeks[Priceidx];

	return value;

}

double StabilityA::getVega(int vegaN, int *vegaidx) {
		 
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
	 
	return m_Ogreeks[twoVegaidx1];

}


double StabilityA::getRho(int rhoN,int *rhoidx) {
	
	 
	double value;

	double *tempirate;  
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		  

	 

	tempirate = new double[m_irateN];
	 
	 


	m_Ogreeks[Rhorfidx] = 0; //rf rho sum
 
	 
	for(i=0; i<MAXRHOS; i++)
		m_Otrfrhos[i] = 0;

	//rf1 rho 시작
	 

	for(i=0; i<m_irateN; i++) 
		tempirate[i] = m_iRFrate[i];
 
		 
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // 평행이동 rho
			for(j=0; j<m_irateN; j++) {
				m_iRFrate[j] = tempirate[j] + RhoPertubation;	
 
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

		value = mc_price();

		m_Otrfrhos[i] = value - m_Ogreeks[Priceidx]; // 10bp term rf rho 
		m_Ogreeks[Rhorfidx] += m_Otrfrhos[i];
			 	 
	}
 

	for(i=0; i<m_irateN; i++) {
		m_iRFrate[i] = tempirate[i]; // RF rate 복원						 
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
 

		
	return m_Ogreeks[Rhorfidx] + m_Ogreeks[Rhodcidx] ;  // 의미 없음 

}


 
 

double StabilityA::mc_price()
{
	long *rsequns;
	double **randoms;
	long nostep;

	double simsval[AsNUM1], simpsval[AsNUM1];
 
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
	 
	 
	for(i=0; i<AsNUM1; i++)
		q[i] = m_divrate[i];

	int si,l,sk,disci;
	int simerFlag;
	double simknockPerf;
	 
	  
	price = 0;
 
 
	if(m_erFlag == 1) {
		for(sk=0; sk<m_matuN; sk++) {
			if(m_matu[sk] < 0) 
				continue;

			price = AMOUNT * MAX(0,1-(m_multiplier*(m_xval[0]-m_erPerf))) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
			break;
		}

	}
	else {
	
		for(si=0; si<m_noiters; si++) { //MC SIM
			simerFlag = 0;
			simprice = 0;
			if(rnd->randseq(1,si,nostep,rsequns) != 0) {
				for(i=0; i<AsNUM1; i++) {
					free(randoms[i]);					 			 
				}
				free(randoms);
				delete dailyIFrf;
				free(rsequns);
				delete dailyZdc;
				delete volskew; 

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
							simpsval[k] = m_psval[k];
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
							simpsval[k] = simsval[k];
							simsval[k] = simsval[k] * exp(mudt[k] + vldt[k]*randoms[k][rsequns[l-1]]);

							if(m_SLimit > 0)							
								simsval[k] = MAX(simpsval[k]*(1-m_SLimit), MIN(simpsval[k]*(1+m_SLimit),simsval[k]));

							simknockPerf = simsval[k]/simpsval[k];

							if(simknockPerf < m_xval[k]) {
								simerFlag = 1;				
								break;
							}
						}

						if(simerFlag == 1)
							break;

					}//if l==0 else

				}//for l

				if(simerFlag == 1) {
					simprice = simprice + AMOUNT * MAX(0,1-(m_multiplier*(m_xval[0]-simknockPerf))) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
					break;
				}
				else {
					simprice = simprice + AMOUNT * m_Cpn[k] * m_CpnFactor[k] *  exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				}				 

				if(sk==m_matuN-1) {
					simprice = simprice + AMOUNT * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				}
				
			} // for sk

			price = price + simprice/m_noiters;

		} // MC SIM
	} // else
		
 
	
		 
		
	//메모리 해제
	for(i=0; i<AsNUM1; i++) {
		free(randoms[i]);
		 
	}
	free(randoms);
	delete dailyIFrf;
	free(rsequns);
	delete dailyZdc;
	delete volskew; 
	 

	// 메모리 해제
 
		 
			
	return price;
	 		



} //mc_price();



void StabilityA::curttimeVol(double *vol, int day, int Asidx)
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
double StabilityA::curttimeFXCorr2(int day)
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




int StabilityA::findinteridx(double x, double *Dx, int DxN)
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

 
 