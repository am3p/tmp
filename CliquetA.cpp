#include "CliquetA.h"
#include <iostream>
#include <fstream>

 

CliquetA::CliquetA(double *sval, double *bval, double prate, double GFloor, double LCap, double LFloor, double NAQ,
	      int matuN, int *matu, double ctime,		  
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int batchFlag, char *kCode,
		  double Xmaxlevel, double Zmesh, int XmeshM, int *TmeshN)		  
		  // bval[0] 에.. 최초 기준가 bval[1]에 첫번째 reset day 종가
		  // payoff = NA * prate * max(sum(TRi), GF)
		  // TRi = max(min(S_i/S_(i-1) -1, C), F)
		  // 
{		
	int i, j;	 

	//입력 데이터 내부 전역 변수로 저장

	m_matuN = matuN;
	m_matu = new int[m_matuN];
	
	m_sval = sval[0];
	 

	m_bval = new double[m_matuN+1]; //S0 ~ SN : N+1개 저장
	  
	m_prate = abs(prate);
	if(prate >=0) {
		m_modeltype = 1;  // 보통	
		m_xalpha = 1.0;  // 추가!! plain Cliquet은 현재 행사가 100% 만!!
	}
	else {
		m_modeltype = -1; // digital
		m_xalpha = sval[1];  // 추가!!
	}

	m_GFloor = GFloor;
	m_LCap = LCap;
	m_LFloor = LFloor;	 
	m_NAQ = NAQ;

	for(i=0; i<m_matuN; i++) {		 
		m_matu[i] = matu[i];
		m_bval[i] = bval[i];				 
	}
	m_bval[m_matuN] = bval[m_matuN];

	m_ctime = ctime;
	/*
	for(i=0; i<m_matuN; i++) {		 		 
		if(m_ctime >= 15 && m_matu[i] == 0)
			m_bval[i+1] = m_sval;  // 오늘이 reset day인데.. 종가를 입력하지 않았을 경우 혹시나 하는 마음에..
	}
	*/

	 
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

	m_Xmaxlevel = Xmaxlevel; // x = log(S/S_ ) 으로 치환해서 풀고 있음 Xmaxlevel = 1인 경우는 S_*exp(1) = S_의 2.7배정도로 푸는 개념 
	                         // C,F가 0.1 수준에서는 충분
	 

	//m_ZFactor = ZFactor; // 0~C, 0~F 까지를 몇등분하느냐의 숫자
	m_Zmesh = Zmesh;
	m_meshXM = XmeshM;
	
	m_meshTN1 = TmeshN[0]; // reset day만기근처
	m_meshTN2 = TmeshN[1]; // 평시



	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/debug_veribook.txt",ios::out);	 

			logsgnuplot << " sval =    "<< sval[0] << endl;

			logsgnuplot << " bval[0] =    "<< bval[0] << endl;
			logsgnuplot << " bval[1] =    "<< bval[1] << endl;
			logsgnuplot << " bval[matuN] =    "<< bval[matuN] << endl;
			logsgnuplot << " prate =    "<< prate << endl; 
			logsgnuplot << " GFloor =    "<< GFloor << endl;
			logsgnuplot << " LCap =    "<< LCap << endl;
			logsgnuplot << " LFloor =    "<< LFloor << endl;
			logsgnuplot << " NAQ =    "<< NAQ << endl;
			logsgnuplot << " matuN =    "<< matuN << endl;
			logsgnuplot << " matu[0] =    "<< matu[0] << endl;
			logsgnuplot << " matu[matuN-1] =    "<< matu[matuN-1] << endl;
			logsgnuplot << " ctime =    "<< ctime << endl;
			logsgnuplot << " irateN =    "<< irateN << endl;
			logsgnuplot << " irateBDay[0] =    "<< irateBDay[0] << endl;
			logsgnuplot << " irateBDay[1] =    "<< irateBDay[1] << endl;
			logsgnuplot << " irateCDay[0] =    "<< irateCDay[0] << endl;
			logsgnuplot << " irateCDay[1] =    "<< irateCDay[1] << endl;
			logsgnuplot << " iRFrate[0] =    "<< iRFrate[0] << endl;
			logsgnuplot << " iRFrate[1] =    "<< iRFrate[1] << endl;
			logsgnuplot << " iDCrate[0] =    "<< iDCrate[0] << endl;
			logsgnuplot << " iDCrate[1] =    "<< iDCrate[1] << endl;
			logsgnuplot << " iIRSrate[0] =    "<< iIRSrate[0] << endl;
			logsgnuplot << " iIRSrate[1] =    "<< iIRSrate[1] << endl;
			logsgnuplot << " divrate =    "<< divrate << endl;
			logsgnuplot << " divbval =    "<< divbval << endl;
			logsgnuplot << " divN =    "<< divN << endl;
			logsgnuplot << " divBDay[0] =    "<< divBDay[0] << endl;
			logsgnuplot << " divBDay[1] =    "<< divBDay[1] << endl;
			logsgnuplot << " divCash[0] =    "<< divCash[0] << endl;
			logsgnuplot << " divCash[1] =    "<< divCash[1] << endl;
			logsgnuplot << " divApy =    "<< divApy << endl;
			logsgnuplot << " voltype =    "<< voltype << endl;
			logsgnuplot << " volbval =    "<< volbval << endl;
			logsgnuplot << " volTN =    "<< volTN << endl;
			logsgnuplot << " volSN =    "<< volSN << endl;
			logsgnuplot << " volBDay[0] =    "<< volBDay[0] << endl;
			logsgnuplot << " volSmness[0] =    "<< volSmness[0] << endl;
			logsgnuplot << " vol[0] =    "<< vol[0] << endl;
			logsgnuplot << " vol[1] =    "<< vol[1] << endl;
			logsgnuplot << " batchFlag =    "<< batchFlag << endl;
			logsgnuplot << " kCode =    "<< kCode << endl;
			logsgnuplot << " Xmaxlevel =    "<< Xmaxlevel << endl;
			logsgnuplot << " Zmesh =    "<< Zmesh << endl;
			logsgnuplot << " XmeshM =    "<< XmeshM << endl;
			logsgnuplot << " TmeshN[0] =    "<< TmeshN[0] << endl;
			logsgnuplot << " TmeshN[1] =    "<< TmeshN[1] << endl;
		 





			logsgnuplot.close();
		 
		}
	  
	 
	  
}

CliquetA::~CliquetA()
{	
	int i;
 

	delete m_irateBDay;
	delete m_irateCDay;  
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
	delete m_bval;

	delete m_kCode;

}

double CliquetA::getValue() {

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


double CliquetA::getVega(int vegaN, int *vegaidx) {

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


double CliquetA::getRho(int rhoN,int *rhoidx) {

	double greeks[MAXGREEKS1];
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
	 
	m_Ogreeks[Rhorfidx] = 0; //rf rho sum
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

double CliquetA::cf_price()
{	
	return 7777;
}

double CliquetA::fd_price(double *greeks)
{
 
	 

	if(m_modeltype == 1) { // 보통
		int i, j, k,  intraLoop;
		int s_idx, term;
		int intraTN, Nx, Nz; // m_meshZL, m_meshXM;

		double X_min, X_max, Z_min, Z_max;
		double **V_Old, **V_New;
		double *X, *Z;

		double dx, dz, dt; 
		double **MAu, **MAd, **MAl; 

		double sum_rate;
		double cur_x;
		int s_z_idx, e_z_idx, s_t_idx, e_t_idx;

		int xcen_idx, zcen_idx;
 
		double value;

		double *SrTemp, *SrVTemp, *SrV1Temp; // greek 계산을 위한 변수 당일값, 다음날 값 저장
 
		int outN;	
		int Baseidx;
		 
		double disc_div;
		double Srv;
		double *volskew;
 
		// 단순 초기화
		i = 0;
		j = 0;
		k = 0;
	 
		intraLoop = 0;
	
		for(s_idx=0; s_idx < m_matuN; s_idx++) {
			if(m_matu[s_idx]>=0)
				break;
		}

		if(m_matu[m_matuN-1] > 0) {
			if(m_matu[s_idx] == 0 && m_ctime >= 15 ) {
				s_idx++;	 
				if(m_bval[s_idx] <= 0.0)				
					m_bval[s_idx] = m_sval;
			}
		}


		cur_x = log(m_sval/m_bval[s_idx]);

		sum_rate = 0.0;
		for(i=0; i<s_idx; i++) {
			sum_rate += MAX(MIN(m_bval[i+1]/m_bval[i] - 1.0, m_LCap), m_LFloor);
		}

		X_min = -1.0*m_Xmaxlevel;
		X_max = 1.0*m_Xmaxlevel; 
		Nx = m_meshXM;
		xcen_idx = Nx/2;

		dx = (X_max - X_min)/((double)Nx);
		X = (double *)malloc((size_t)((Nx+1)*sizeof(double)));

		for(i=0; i<=Nx; i++)
			X[i] = X_min + i*dx;

	
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
	
		Z_min = sum_rate + m_LFloor*(m_matuN-s_idx); 
		Z_max = sum_rate + m_LCap*(m_matuN-s_idx);
		Nz = (int)floor((Z_max-Z_min)/m_Zmesh + 0.5); //반올림 
		//Nz = 2*m_ZFactor*(m_matuN-s_idx);
			
		dz = (Z_max-Z_min)/((double)Nz);
	 

		Z = (double *)malloc((size_t)((Nz+1)*sizeof(double)));

		for(i=0; i<=Nz; i++) 
			Z[i] = Z_min + i*dz;

		zcen_idx = (int)floor((sum_rate - Z_min)/dz);
 
		MAu = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		MAd = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		MAl = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		for(i=0; i<=Nx; i++) {
			MAu[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
			MAd[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
			MAl[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
		}

	 

		V_Old = (double **)malloc((size_t)((Nx+1)*sizeof(double *)));
		V_New = (double **)malloc((size_t)((Nx+1)*sizeof(double *)));
		for(i=0; i<=Nx; i++) {
			V_Old[i] = (double *)malloc((size_t)((Nz+1)*sizeof(double)));
			V_New[i] = (double *)malloc((size_t)((Nz+1)*sizeof(double)));
		}
		
		volskew = (double*)malloc((size_t)(m_volSN)*sizeof(double));
 
	  
	 	
 
		double tempvalue;
		int NGP;
		double alpha, alphaz;
		double xtemp;

		// 만기 payoff setting
 
		for(j=0; j<=Nz; j++)
			V_New[xcen_idx][j] = AMOUNT * m_prate * (m_NAQ + MAX(MIN(MAX(Z[j],m_LFloor*m_matuN),m_LCap*m_matuN),m_GFloor));
	 
	/*		
		if(m_matu[m_matuN-1] == 0) {
			double tempvalue, alpha;
			int NGP;
			for(i=0; i<=Nx; i++) {
				tempvalue = sum_rate + MAX(MIN(exp(X[i])-1,m_LCap),m_LFloor);
				NGP = (int)floor((tempvalue-Z_min)/dz);
				if(NGP<0) NGP = 0;
				if(NGP>=Nz) NGP = Nz-1;

				alpha = (tempvalue-Z[NGP])/(Z[NGP+1]-Z[NGP]);
				V_Old[i][zcen_idx] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
			}

			for(i=0; i<=Nx; i++)
				V_New[i][zcen_idx] = V_Old[i][zcen_idx];

		}


		if(m_curtgreek == 1 && m_matu[m_matuN-1] == 0) { // 기본 그릭	Theta를 위한 +1일 저장	 		
		
			for(i=0; i<outN + 4; i++) {		
				xtemp = log(SrTemp[i]/m_bval[s_idx]);					 
			
				NGP = (int)floor((xtemp-X_min)/dx);			
				if(NGP<0) NGP = 0;			
				if(NGP>=Nx) NGP = Nx-1;
			
				alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);			
				SrV1Temp[i] = (1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx];
			} //for(i)
		}//if(m_curgreek)
	*/  

		double r, q, v, rd;
		double *dailyIFrf, *dailyIFdc;	
		// intra Day에서는 입력 파라미터 고정
		// 하루마다 해당 파라미터 계산
		q = m_divrate;

		if(m_matu[m_matuN-1]>0) { // 현재 curtN 는 m_matuN <-- 만기해당..
			dailyIFrf = (double *)malloc((size_t)(m_matu[m_matuN-1])*sizeof(double));
			dailyIFdc = (double *)malloc((size_t)(m_matu[m_matuN-1])*sizeof(double));
		
			for(i=0; i<m_matu[m_matuN-1]; i++) {
				dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
				dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
			}
		} //if(m_matu[curtN-1]>0)
	
	   


		// time.. 계산 시작
	
		for(term=m_matuN-1; term>=s_idx; term--) {
			//s_z_idx = m_ZFactor*(m_matuN-term); ? m_ZFactor*(term - m_matuN) 이상하넹..
			//e_z_idx = Nz-s_z_idx;
			s_z_idx = (int)floor((sum_rate + m_LFloor*(term-s_idx) - Z_min)/dz);
			e_z_idx = (int)floor((sum_rate + m_LCap*(term-s_idx) - Z_min)/dz) +1; // add

			for(j=s_z_idx; j<=e_z_idx; j++) {
				for(i=0; i<=Nx; i++) {
					tempvalue = Z[j] + MAX(MIN(exp(X[i])-1,m_LCap),m_LFloor);
					NGP = (int)floor((tempvalue-Z_min)/dz);
					if(NGP<0) NGP = 0;
					if(NGP>=Nz) NGP = Nz-1;

					alpha = (tempvalue-Z[NGP])/(Z[NGP+1]-Z[NGP]);		
					V_Old[i][j] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
				} // for(i)
			} // for(j)
	
		///////////////////////////////////////////--------------Debug	Mode e	
			for(i=0; i<=Nx; i++) 
				for(j=s_z_idx; j<=e_z_idx; j++)
					V_New[i][j] = V_Old[i][j]; 
	 	
 
			   
			//if(m_curtgreek == 1 && (m_matu[m_matuN-1] == 0 || m_matu[m_matuN-1] == 1)) { // 기본 그릭	Theta를 위한 +1일 저장			 	
			if(m_curtgreek == 1 && (m_matu[m_matuN-1] == 0 || m_matu[term] == 1)) { // 기본 그릭	Theta를 위한 +1일 저장	 
				for(i=0; i<outN + 4; i++) {		
					xtemp = log(SrTemp[i]/m_bval[s_idx]);
					NGP = (int)floor((xtemp - X_min)/dx);
					if(NGP<0) NGP = 0;
					if(NGP>=Nx) NGP = Nx-1;

					alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
					alphaz = (sum_rate - Z[zcen_idx])/dz;
					SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
				}
			}

			if(term>s_idx) {
				e_t_idx = m_matu[term-1];
			}
			else {
				e_t_idx = 0;
			}

			for(s_t_idx = m_matu[term]-1; s_t_idx >= e_t_idx; s_t_idx--) {

				if(s_t_idx == m_matu[term]-1)
					intraTN = m_meshTN1;
				else
					intraTN = m_meshTN2;
				dt = (1.0/YFactor)/intraTN;

				r = dailyIFrf[s_t_idx];
				rd = dailyIFdc[s_t_idx];

				if(m_voltype == 0) 
					volskew[0] = m_vol[0][0];
				else
					curttimeVol(volskew, s_t_idx);
			
				//A*V_new = V_old 생성 Z[j]마다.. 
				for(j=s_z_idx; j<=e_z_idx; j++) {
					for(i=1; i<=Nx-1; i++) {
						if(m_voltype == 2) {
							Srv = realS(X[i],Z[j],m_bval[s_idx],sum_rate)/m_volbval;
							v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0);
						}
						else {
							v = volskew[0];
						}

						MAl[i][j] = (r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
						MAd[i][j] = 1 + rd*dt + v*v*dt/(dx*dx);
						MAu[i][j] = -(r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
					} // for(i)
									
					//B.C
					//MAd[1][j] = MAd[1][j] + 2*MAl[1][j];				
					//MAu[1][j] = MAu[1][j] - MAl[1][j];
					//log변환 B.C
					MAd[1][j] = MAd[1][j] + 4*(1+dx)/(2+3*dx) * MAl[1][j];
					MAu[1][j] = MAu[1][j] - (2+dx)/(2+3*dx) * MAl[1][j];
				
					//MAl[Nx-1][j] = MAl[Nx-1][j] - MAu[Nx-1][j];				
					//MAd[Nx-1][j] = MAd[Nx-1][j] + 2*MAu[Nx-1][j];
					MAl[Nx-1][j] = MAl[Nx-1][j] - (2-dx)/(2-3*dx) * MAu[Nx-1][j];				
					MAd[Nx-1][j] = MAd[Nx-1][j] + 4*(1-dx)/(2-3*dx) *MAu[Nx-1][j];

					//LU
					for(i=2; i<=Nx-1; i++) {
						MAl[i][j] = MAl[i][j]/MAd[i-1][j];
						MAd[i][j] = MAd[i][j] - MAl[i][j]*MAu[i-1][j];
					} 
				}//for(j)  // 행렬 완료
				 

				for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {

					//Vold에 Vnew 복사	
					for(i=0; i<=Nx; i++) 			
						for(j=s_z_idx; j<=e_z_idx; j++)				
							V_Old[i][j] = V_New[i][j]; 

					for(j=s_z_idx; j<=e_z_idx; j++) {

						//L*y = b forward 풀기
						V_New[1][j] = V_Old[1][j]; // Vold 복사 했으므로.. 꼭 필요한건 아니지만.. 로직상..
						for(i=2; i<=Nx-1; i++) 
							V_New[i][j] = V_Old[i][j] - MAl[i][j]*V_New[i-1][j];

						//U*x = y back풀기
						V_New[Nx-1][j] = V_New[Nx-1][j]/MAd[Nx-1][j];
						for(i=Nx-2; i>=1; i--)
							V_New[i][j] = (V_New[i][j] - MAu[i][j]*V_New[i+1][j])/MAd[i][j];

						//V_New[0][j] = 2*V_New[1][j] - V_New[2][j];
						//V_New[Nx][j] = 2*V_New[Nx-1][j] - V_New[Nx-2][j];
						//log 변환 B.C 적용
						V_New[0][j] = 4*(1+dx)/(2+3*dx) * V_New[1][j] - (2+dx)/(2+3*dx) * V_New[2][j];
						V_New[Nx][j] = 4*(1-dx)/(2-3*dx) * V_New[Nx-1][j] - (2-dx)/(2-3*dx) * V_New[Nx-2][j];
					
					} //for(j) z풀기
				}//for(intraLoop)					 

				//이산배당 체크
				disc_div = 0.0;
				for(i=0; i<m_divN; i++) {
					if(s_t_idx == (m_divBDay[i]-1)) {
						disc_div = m_divCash[i];
						break;
					}
				}

				// 적용
				if(disc_div > 0.0) {
					double divS[2];
					double pole[] = {0, 1};
					double alpha, divalpha;
					double Srtemp;
					int NGP;

					//Vold에 Vnew 복사	
					for(i=0; i<=Nx; i++) 			
						for(j=s_z_idx; j<=e_z_idx; j++)				
							V_Old[i][j] = V_New[i][j]; 

					divS[0] = (m_divbval*m_divApy - disc_div);
					divS[1] = (m_divbval*m_divApy + disc_div);

					for(j=s_z_idx; j<=e_z_idx; j++) {
						for(i=0; i<=Nx; i++) {
							Srtemp = realS(X[i],Z[j],m_bval[s_idx],sum_rate);
							divalpha = interp1(Srtemp,divS,pole,2,1,0);
							xtemp = log((Srtemp - divalpha*disc_div)/realS(0,Z[j],m_bval[s_idx],sum_rate));

							NGP = (int)floor((xtemp - X_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nx) NGP = Nx-1;

							alpha = (xtemp - X[NGP])/(X[NGP+1]-X[NGP]);

							V_New[i][j] = (1.0-alpha)*V_Old[NGP][j] + alpha*V_Old[NGP+1][j];
						}//for(i)
					}//for(j)

				}//if(disc_div)

				if(m_curtgreek == 1 && s_t_idx == 1) { // 기본 그릭 Theta를 위한 +1일 저장							
			
					if(m_matu[s_idx] == 0) {	
						double ztemp;

						for(i=0; i<outN + 4; i++) {		
							ztemp = Z[zcen_idx] + MAX(MIN(SrTemp[i]/m_bval[s_idx]-1,m_LCap),m_LFloor);
							NGP = (int)floor((ztemp - Z_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nz) NGP = Nz-1;

							alpha = (ztemp-Z[NGP])/(Z[NGP+1]-Z[NGP]);
							SrV1Temp[i] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
						}// for(i)
					}
					else {			 
						for(i=0; i<outN + 4; i++) {		
							xtemp = log(SrTemp[i]/m_bval[s_idx]);
							NGP = (int)floor((xtemp - X_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nx) NGP = Nx-1;

							alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
							alphaz = (sum_rate - Z[zcen_idx])/dz;				 
							//SrV1Temp[i] = (1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx];
							SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
						}// for(i)
					}

				}// if(m_curtgreek)

			}//for(s_t_idx)
		}//for(term)

	 

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
			logsgnuplot.open("c:/logland/debug_clqagnufdprice.txt",ios::out);	 
			for(i=0; i<=Nx; i++) {
				alphaz = (sum_rate - Z[zcen_idx])/dz;			
				logsgnuplot << exp(X[i])*m_bval[s_idx] <<"     "<< (1-alphaz)*V_New[i][zcen_idx] + alphaz*V_New[i][zcen_idx+1] << endl;
			}		 
			logsgnuplot.close();
		 
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	 
	 
		if(m_curtgreek == 1 || m_curtgreek == 2) {
				
	 		
			 
			for(i=0; i<outN + 4; i++) {		
				xtemp = log(SrTemp[i]/m_bval[s_idx]);
				NGP = (int)floor((xtemp - X_min)/dx);
				if(NGP<0) NGP = 0;
				if(NGP>=Nx) NGP = Nx-1;

				alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
				alphaz = (sum_rate-Z[zcen_idx])/dz;
				SrVTemp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
			}// for(i)

		 
		
			Baseidx = SoutRange + 2; // m_sval에 해당하는 idx 위 아래 2개씩 더 계산해서.. +2 필요
		}

		if(m_curtgreek == 1) { // 기본 그릭  
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);
			alphaz = (sum_rate - Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + (alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
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
		
				if(m_matu[s_idx] == 0 && m_ctime < 15) {		
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

						
					//A*V_new = V_old 생성 Z[j]마다.. 
					for(j=s_z_idx; j<=e_z_idx; j++) {
						for(i=1; i<=Nx-1; i++) {
							if(m_voltype == 2) {
								Srv = realS(X[i],Z[j],m_bval[s_idx],sum_rate)/m_volbval;
								v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0);
							}
							else {
								v = volskew[0];
							}

							MAl[i][j] = (r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
							MAd[i][j] = 1 + rd*dt + v*v*dt/(dx*dx);
							MAu[i][j] = -(r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
						} // for(i)
									
						//B.C
						//MAd[1][j] = MAd[1][j] + 2*MAl[1][j];				
						//MAu[1][j] = MAu[1][j] - MAl[1][j];
								
						//MAl[Nx-1][j] = MAl[Nx-1][j] - MAu[Nx-1][j];				
						//MAd[Nx-1][j] = MAd[Nx-1][j] + 2*MAu[Nx-1][j];

						MAd[1][j] = MAd[1][j] + 4*(1+dx)/(2+3*dx) * MAl[1][j];				
						MAu[1][j] = MAu[1][j] - (2+dx)/(2+3*dx) * MAl[1][j];
								
						MAl[Nx-1][j] = MAl[Nx-1][j] - (2-dx)/(2-3*dx) * MAu[Nx-1][j];								
						MAd[Nx-1][j] = MAd[Nx-1][j] + 4*(1-dx)/(2-3*dx) *MAu[Nx-1][j];

						//LU
						for(i=2; i<=Nx-1; i++) {
							MAl[i][j] = MAl[i][j]/MAd[i-1][j];
							MAd[i][j] = MAd[i][j] - MAl[i][j]*MAu[i-1][j];
						} 
					}//for(j)  // 행렬 완료
 
					///////////////

  
		
					for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 확인		
 
						//Vold에 Vnew 복사	
						for(i=0; i<=Nx; i++) 			
							for(j=s_z_idx; j<=e_z_idx; j++)				
								V_Old[i][j] = V_New[i][j]; 

						for(j=s_z_idx; j<=e_z_idx; j++) {

							//L*y = b forward 풀기
							V_New[1][j] = V_Old[1][j]; // Vold 복사 했으므로.. 꼭 필요한건 아니지만.. 로직상..
							for(i=2; i<=Nx-1; i++) 
								V_New[i][j] = V_Old[i][j] - MAl[i][j]*V_New[i-1][j];

							//U*x = y back풀기
							V_New[Nx-1][j] = V_New[Nx-1][j]/MAd[Nx-1][j];
							for(i=Nx-2; i>=1; i--)
								V_New[i][j] = (V_New[i][j] - MAu[i][j]*V_New[i+1][j])/MAd[i][j];

							//V_New[0][j] = 2*V_New[1][j] - V_New[2][j];
							//V_New[Nx][j] = 2*V_New[Nx-1][j] - V_New[Nx-2][j];							
							V_New[0][j] = 4*(1+dx)/(2+3*dx) * V_New[1][j] - (2+dx)/(2+3*dx) * V_New[2][j];					
							V_New[Nx][j] = 4*(1-dx)/(2-3*dx) * V_New[Nx-1][j] - (2-dx)/(2-3*dx) * V_New[Nx-2][j];
					
						} //for(j) z풀기


						if(m_curtgreek == 1) { // 만기에
								 
							for(i=0; i<outN + 4; i++) {		
								xtemp = log(SrTemp[i]/m_bval[s_idx]);
								NGP = (int)floor((xtemp - X_min)/dx);
								if(NGP<0) NGP = 0;
								if(NGP>=Nx) NGP = Nx-1;

								alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
								alphaz = (sum_rate-Z[zcen_idx])/dz;
								SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) +(alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
							}// for(i)
		 
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
				
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);
			alphaz = (sum_rate - Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + (alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
			//value = interp1(slevel,S,newV,SM+1,SmeshType,1);	

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
				
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);

			alphaz = (sum_rate-Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) +(alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;

			//value = interp1(slevel,S,newV,SM+1,SmeshType,1); // value - m_Ogreeks[0] 을 이용하기 위하여
		}
	 

		free(X);
		free(Z);

		for(i=0;i<=Nx; i++) {
			free(V_New[i]);
			free(V_Old[i]);
			free(MAu[i]);
			free(MAd[i]);
			free(MAl[i]);
		}	 	 
		free(V_New);
		free(V_Old);
		free(MAu);
		free(MAd);
		free(MAl);
	 
 
		free(volskew);

		return value;
	}
	else { // if m_modeltype   // here digital
		int i, j, k,  intraLoop;
		int s_idx, term;
		int intraTN, Nx, Nz; // m_meshZL, m_meshXM;

		double X_min, X_max, Z_min, Z_max;
		double **V_Old, **V_New;
		double *X, *Z;

		double dx, dz, dt; 
		double **MAu, **MAd, **MAl; 

		double sum_rate;
		double cur_x;
		int s_z_idx, e_z_idx, s_t_idx, e_t_idx;

		int xcen_idx, zcen_idx;
 
		double value;

		double *SrTemp, *SrVTemp, *SrV1Temp; // greek 계산을 위한 변수 당일값, 다음날 값 저장
 
		int outN;	
		int Baseidx;
		 
		double disc_div;
		double Srv;
		double *volskew;
 
		// 단순 초기화
		i = 0;
		j = 0;
		k = 0;
	 
		intraLoop = 0;
	
		for(s_idx=0; s_idx < m_matuN; s_idx++) {
			if(m_matu[s_idx]>=0)
				break;
		}

		if(m_matu[m_matuN-1] > 0) {
			if(m_matu[s_idx] == 0 && m_ctime >= 15) {
				s_idx++;
				m_bval[s_idx] = m_sval;
			}
		}

		 
		//cur_x = log(m_sval/m_bval[s_idx]);
		cur_x = log(m_sval/(m_xalpha*m_bval[s_idx]));

		sum_rate = 0.0;
		for(i=0; i<s_idx; i++) { 
			//if(m_bval[i+1]/m_bval[i] - 1.0 >= 0) 				
			if(m_bval[i+1]/m_bval[i] - m_xalpha >= 0) 
				sum_rate += m_LCap;
			else
				sum_rate += m_LFloor;

			//sum_rate += MAX(MIN(m_bval[i+1]/m_bval[i] - 1.0, m_LCap), m_LFloor);
		}

		X_min = -1.0*m_Xmaxlevel;
		X_max = 1.0*m_Xmaxlevel; 
		Nx = m_meshXM;
		xcen_idx = Nx/2;

		dx = (X_max - X_min)/((double)Nx);
		X = (double *)malloc((size_t)((Nx+1)*sizeof(double)));

		for(i=0; i<=Nx; i++)
			X[i] = X_min + i*dx;

	
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
	
		Z_min = sum_rate + m_LFloor*(m_matuN-s_idx); 
		Z_max = sum_rate + m_LCap*(m_matuN-s_idx);
		Nz = (int)floor((Z_max-Z_min)/m_Zmesh + 0.5); //반올림 
		//Nz = 2*m_ZFactor*(m_matuN-s_idx);
			
		dz = (Z_max-Z_min)/((double)Nz);
	 

		Z = (double *)malloc((size_t)((Nz+1)*sizeof(double)));

		for(i=0; i<=Nz; i++) 
			Z[i] = Z_min + i*dz;

		zcen_idx = (int)floor((sum_rate - Z_min)/dz);
 
		MAu = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		MAd = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		MAl = (double **)malloc((size_t)(Nx+1)*sizeof(double*));
		for(i=0; i<=Nx; i++) {
			MAu[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
			MAd[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
			MAl[i] = (double *)malloc((size_t)(Nz+1)*sizeof(double));
		}

	 

		V_Old = (double **)malloc((size_t)((Nx+1)*sizeof(double *)));
		V_New = (double **)malloc((size_t)((Nx+1)*sizeof(double *)));
		for(i=0; i<=Nx; i++) {
			V_Old[i] = (double *)malloc((size_t)((Nz+1)*sizeof(double)));
			V_New[i] = (double *)malloc((size_t)((Nz+1)*sizeof(double)));
		}
		
		volskew = (double*)malloc((size_t)(m_volSN)*sizeof(double));
 
	  
	 	
 
		double tempvalue;
		int NGP;
		double alpha, alphaz;
		double xtemp;

		// 만기 payoff setting
 
		for(j=0; j<=Nz; j++)
			V_New[xcen_idx][j] = AMOUNT * m_prate * (m_NAQ + MAX(MIN(MAX(Z[j],m_LFloor*m_matuN),m_LCap*m_matuN),m_GFloor));
	 
	/*		
		if(m_matu[m_matuN-1] == 0) {
			double tempvalue, alpha;
			int NGP;
			for(i=0; i<=Nx; i++) {
				tempvalue = sum_rate + MAX(MIN(exp(X[i])-1,m_LCap),m_LFloor);
				NGP = (int)floor((tempvalue-Z_min)/dz);
				if(NGP<0) NGP = 0;
				if(NGP>=Nz) NGP = Nz-1;

				alpha = (tempvalue-Z[NGP])/(Z[NGP+1]-Z[NGP]);
				V_Old[i][zcen_idx] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
			}

			for(i=0; i<=Nx; i++)
				V_New[i][zcen_idx] = V_Old[i][zcen_idx];

		}


		if(m_curtgreek == 1 && m_matu[m_matuN-1] == 0) { // 기본 그릭	Theta를 위한 +1일 저장	 		
		
			for(i=0; i<outN + 4; i++) {		
				xtemp = log(SrTemp[i]/m_bval[s_idx]);					 
			
				NGP = (int)floor((xtemp-X_min)/dx);			
				if(NGP<0) NGP = 0;			
				if(NGP>=Nx) NGP = Nx-1;
			
				alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);			
				SrV1Temp[i] = (1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx];
			} //for(i)
		}//if(m_curgreek)
	*/  

		double r, q, v, rd;
		double *dailyIFrf, *dailyIFdc;	
		// intra Day에서는 입력 파라미터 고정
		// 하루마다 해당 파라미터 계산
		q = m_divrate;

		if(m_matu[m_matuN-1]>0) { // 현재 curtN 는 m_matuN <-- 만기해당..
			dailyIFrf = (double *)malloc((size_t)(m_matu[m_matuN-1])*sizeof(double));
			dailyIFdc = (double *)malloc((size_t)(m_matu[m_matuN-1])*sizeof(double));
		
			for(i=0; i<m_matu[m_matuN-1]; i++) {
				dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
				dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
			}
		} //if(m_matu[curtN-1]>0)
	
	   


		// time.. 계산 시작
	
		for(term=m_matuN-1; term>=s_idx; term--) {
			//s_z_idx = m_ZFactor*(m_matuN-term); ? m_ZFactor*(term - m_matuN) 이상하넹..
			//e_z_idx = Nz-s_z_idx;
			s_z_idx = (int)floor((sum_rate + m_LFloor*(term-s_idx) - Z_min)/dz);
			e_z_idx = (int)floor((sum_rate + m_LCap*(term-s_idx) - Z_min)/dz) +1; // add

			for(j=s_z_idx; j<=e_z_idx; j++) {
				for(i=0; i<=Nx; i++) {
					if(exp(X[i])-1 >=0 )
						tempvalue = Z[j] + m_LCap;
					else
						tempvalue = Z[j] + m_LFloor;
					//tempvalue = Z[j] + MAX(MIN(exp(X[i])-1,m_LCap),m_LFloor);
					NGP = (int)floor((tempvalue-Z_min)/dz);
					if(NGP<0) NGP = 0;
					if(NGP>=Nz) NGP = Nz-1;

					alpha = (tempvalue-Z[NGP])/(Z[NGP+1]-Z[NGP]);		
					V_Old[i][j] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
				} // for(i)
			} // for(j)
	
		///////////////////////////////////////////--------------Debug	Mode e	
			for(i=0; i<=Nx; i++) 
				for(j=s_z_idx; j<=e_z_idx; j++)
					V_New[i][j] = V_Old[i][j]; 
	 	
 
			   
			//if(m_curtgreek == 1 && (m_matu[m_matuN-1] == 0 || m_matu[m_matuN-1] == 1)) { // 기본 그릭	Theta를 위한 +1일 저장			 	
			if(m_curtgreek == 1 && (m_matu[m_matuN-1] == 0 || m_matu[term] == 1)) { // 기본 그릭	Theta를 위한 +1일 저장	 
				for(i=0; i<outN + 4; i++) {		
					xtemp = log(SrTemp[i]/(m_xalpha*m_bval[s_idx]));
					NGP = (int)floor((xtemp - X_min)/dx);
					if(NGP<0) NGP = 0;
					if(NGP>=Nx) NGP = Nx-1;

					alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
					alphaz = (sum_rate - Z[zcen_idx])/dz;
					SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
				}
			}

			if(term>s_idx) {
				e_t_idx = m_matu[term-1];
			}
			else {
				e_t_idx = 0;
			}

			for(s_t_idx = m_matu[term]-1; s_t_idx >= e_t_idx; s_t_idx--) {

				if(s_t_idx == m_matu[term]-1)
					intraTN = m_meshTN1;
				else
					intraTN = m_meshTN2;
				dt = (1.0/YFactor)/intraTN;

				r = dailyIFrf[s_t_idx];
				rd = dailyIFdc[s_t_idx];

				if(m_voltype == 0) 
					volskew[0] = m_vol[0][0];
				else
					curttimeVol(volskew, s_t_idx);
			
				//A*V_new = V_old 생성 Z[j]마다.. 
				for(j=s_z_idx; j<=e_z_idx; j++) {
					for(i=1; i<=Nx-1; i++) {
						if(m_voltype == 2) {
							Srv = realS(X[i],Z[j],m_bval[s_idx],sum_rate)/m_volbval;
							v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0);
						}
						else {
							v = volskew[0];
						}

						MAl[i][j] = (r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
						MAd[i][j] = 1 + rd*dt + v*v*dt/(dx*dx);
						MAu[i][j] = -(r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
					} // for(i)
									
					//B.C
					//MAd[1][j] = MAd[1][j] + 2*MAl[1][j];				
					//MAu[1][j] = MAu[1][j] - MAl[1][j];
					//log변환 B.C
					MAd[1][j] = MAd[1][j] + 4*(1+dx)/(2+3*dx) * MAl[1][j];
					MAu[1][j] = MAu[1][j] - (2+dx)/(2+3*dx) * MAl[1][j];
				
					//MAl[Nx-1][j] = MAl[Nx-1][j] - MAu[Nx-1][j];				
					//MAd[Nx-1][j] = MAd[Nx-1][j] + 2*MAu[Nx-1][j];
					MAl[Nx-1][j] = MAl[Nx-1][j] - (2-dx)/(2-3*dx) * MAu[Nx-1][j];				
					MAd[Nx-1][j] = MAd[Nx-1][j] + 4*(1-dx)/(2-3*dx) *MAu[Nx-1][j];

					//LU
					for(i=2; i<=Nx-1; i++) {
						MAl[i][j] = MAl[i][j]/MAd[i-1][j];
						MAd[i][j] = MAd[i][j] - MAl[i][j]*MAu[i-1][j];
					} 
				}//for(j)  // 행렬 완료
				 

				for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {

					//Vold에 Vnew 복사	
					for(i=0; i<=Nx; i++) 			
						for(j=s_z_idx; j<=e_z_idx; j++)				
							V_Old[i][j] = V_New[i][j]; 

					for(j=s_z_idx; j<=e_z_idx; j++) {

						//L*y = b forward 풀기
						V_New[1][j] = V_Old[1][j]; // Vold 복사 했으므로.. 꼭 필요한건 아니지만.. 로직상..
						for(i=2; i<=Nx-1; i++) 
							V_New[i][j] = V_Old[i][j] - MAl[i][j]*V_New[i-1][j];

						//U*x = y back풀기
						V_New[Nx-1][j] = V_New[Nx-1][j]/MAd[Nx-1][j];
						for(i=Nx-2; i>=1; i--)
							V_New[i][j] = (V_New[i][j] - MAu[i][j]*V_New[i+1][j])/MAd[i][j];

						//V_New[0][j] = 2*V_New[1][j] - V_New[2][j];
						//V_New[Nx][j] = 2*V_New[Nx-1][j] - V_New[Nx-2][j];
						//log 변환 B.C 적용
						V_New[0][j] = 4*(1+dx)/(2+3*dx) * V_New[1][j] - (2+dx)/(2+3*dx) * V_New[2][j];
						V_New[Nx][j] = 4*(1-dx)/(2-3*dx) * V_New[Nx-1][j] - (2-dx)/(2-3*dx) * V_New[Nx-2][j];
					
					} //for(j) z풀기
				}//for(intraLoop)					 

				//이산배당 체크
				disc_div = 0.0;
				for(i=0; i<m_divN; i++) {
					if(s_t_idx == (m_divBDay[i]-1)) {
						disc_div = m_divCash[i];
						break;
					}
				}

				// 적용
				if(disc_div > 0.0) {
					double divS[2];
					double pole[] = {0, 1};
					double alpha, divalpha;
					double Srtemp;
					int NGP;

					//Vold에 Vnew 복사	
					for(i=0; i<=Nx; i++) 			
						for(j=s_z_idx; j<=e_z_idx; j++)				
							V_Old[i][j] = V_New[i][j]; 

					divS[0] = (m_divbval*m_divApy - disc_div);
					divS[1] = (m_divbval*m_divApy + disc_div);

					for(j=s_z_idx; j<=e_z_idx; j++) {
						for(i=0; i<=Nx; i++) {
							Srtemp = realS(X[i],Z[j],m_bval[s_idx],sum_rate);
							divalpha = interp1(Srtemp,divS,pole,2,1,0);
							xtemp = log((Srtemp - divalpha*disc_div)/realS(0,Z[j],m_bval[s_idx],sum_rate));

							NGP = (int)floor((xtemp - X_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nx) NGP = Nx-1;

							alpha = (xtemp - X[NGP])/(X[NGP+1]-X[NGP]);

							V_New[i][j] = (1.0-alpha)*V_Old[NGP][j] + alpha*V_Old[NGP+1][j];
						}//for(i)
					}//for(j)

				}//if(disc_div)

				if(m_curtgreek == 1 && s_t_idx == 1) { // 기본 그릭 Theta를 위한 +1일 저장							
			
					if(m_matu[s_idx] == 0) {	
						double ztemp;

						for(i=0; i<outN + 4; i++) {	
							if(SrTemp[i]/m_bval[s_idx]-m_xalpha >= 0)
								ztemp = Z[zcen_idx] + m_LCap;
							else
								ztemp = Z[zcen_idx] + m_LFloor;

							//ztemp = Z[zcen_idx] + MAX(MIN(SrTemp[i]/m_bval[s_idx]-1,m_LCap),m_LFloor);
							NGP = (int)floor((ztemp - Z_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nz) NGP = Nz-1;

							alpha = (ztemp-Z[NGP])/(Z[NGP+1]-Z[NGP]);
							SrV1Temp[i] = (1.0-alpha)*V_New[xcen_idx][NGP] + alpha*V_New[xcen_idx][NGP+1];
						}// for(i)
					}
					else {			 
						for(i=0; i<outN + 4; i++) {		
							xtemp = log(SrTemp[i]/(m_xalpha*m_bval[s_idx]));
							NGP = (int)floor((xtemp - X_min)/dx);
							if(NGP<0) NGP = 0;
							if(NGP>=Nx) NGP = Nx-1;

							alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
							alphaz = (sum_rate - Z[zcen_idx])/dz;				 
							//SrV1Temp[i] = (1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx];
							SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
						}// for(i)
					}

				}// if(m_curtgreek)

			}//for(s_t_idx)
		}//for(term)

	 

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
			logsgnuplot.open("c:/logland/debug_clqagnufdprice.txt",ios::out);	 
			for(i=0; i<=Nx; i++) {
				alphaz = (sum_rate - Z[zcen_idx])/dz;			
				logsgnuplot << exp(X[i])*m_xalpha*m_bval[s_idx] <<"     "<< (1-alphaz)*V_New[i][zcen_idx] + alphaz*V_New[i][zcen_idx+1] << endl;
			}		 
			logsgnuplot.close();
		 
		}
		///////////////////////////////////////////--------------Debug	Mode e	
	 
	 
		if(m_curtgreek == 1 || m_curtgreek == 2) {
				
	 		
			 
			for(i=0; i<outN + 4; i++) {		
				xtemp = log(SrTemp[i]/(m_xalpha*m_bval[s_idx]));
				NGP = (int)floor((xtemp - X_min)/dx);
				if(NGP<0) NGP = 0;
				if(NGP>=Nx) NGP = Nx-1;

				alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
				alphaz = (sum_rate-Z[zcen_idx])/dz;
				SrVTemp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + alphaz*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
			}// for(i)

		 
		
			Baseidx = SoutRange + 2; // m_sval에 해당하는 idx 위 아래 2개씩 더 계산해서.. +2 필요
		}

		if(m_curtgreek == 1) { // 기본 그릭  
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);
			alphaz = (sum_rate - Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + (alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]);
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
		
				if(m_matu[s_idx] == 0 && m_ctime < 15) {		
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

						
					//A*V_new = V_old 생성 Z[j]마다.. 
					for(j=s_z_idx; j<=e_z_idx; j++) {
						for(i=1; i<=Nx-1; i++) {
							if(m_voltype == 2) {
								Srv = realS(X[i],Z[j],m_bval[s_idx],sum_rate)/m_volbval;
								v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0);
							}
							else {
								v = volskew[0];
							}

							MAl[i][j] = (r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
							MAd[i][j] = 1 + rd*dt + v*v*dt/(dx*dx);
							MAu[i][j] = -(r-q-v*v/2.0)*dt/(2*dx) - (v*v/2.0)*(dt/(dx*dx));
						} // for(i)
									
						//B.C
						//MAd[1][j] = MAd[1][j] + 2*MAl[1][j];				
						//MAu[1][j] = MAu[1][j] - MAl[1][j];
								
						//MAl[Nx-1][j] = MAl[Nx-1][j] - MAu[Nx-1][j];				
						//MAd[Nx-1][j] = MAd[Nx-1][j] + 2*MAu[Nx-1][j];

						MAd[1][j] = MAd[1][j] + 4*(1+dx)/(2+3*dx) * MAl[1][j];				
						MAu[1][j] = MAu[1][j] - (2+dx)/(2+3*dx) * MAl[1][j];
								
						MAl[Nx-1][j] = MAl[Nx-1][j] - (2-dx)/(2-3*dx) * MAu[Nx-1][j];								
						MAd[Nx-1][j] = MAd[Nx-1][j] + 4*(1-dx)/(2-3*dx) *MAu[Nx-1][j];

						//LU
						for(i=2; i<=Nx-1; i++) {
							MAl[i][j] = MAl[i][j]/MAd[i-1][j];
							MAd[i][j] = MAd[i][j] - MAl[i][j]*MAu[i-1][j];
						} 
					}//for(j)  // 행렬 완료
 
					///////////////

  
		
					for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 확인		
 
						//Vold에 Vnew 복사	
						for(i=0; i<=Nx; i++) 			
							for(j=s_z_idx; j<=e_z_idx; j++)				
								V_Old[i][j] = V_New[i][j]; 

						for(j=s_z_idx; j<=e_z_idx; j++) {

							//L*y = b forward 풀기
							V_New[1][j] = V_Old[1][j]; // Vold 복사 했으므로.. 꼭 필요한건 아니지만.. 로직상..
							for(i=2; i<=Nx-1; i++) 
								V_New[i][j] = V_Old[i][j] - MAl[i][j]*V_New[i-1][j];

							//U*x = y back풀기
							V_New[Nx-1][j] = V_New[Nx-1][j]/MAd[Nx-1][j];
							for(i=Nx-2; i>=1; i--)
								V_New[i][j] = (V_New[i][j] - MAu[i][j]*V_New[i+1][j])/MAd[i][j];

							//V_New[0][j] = 2*V_New[1][j] - V_New[2][j];
							//V_New[Nx][j] = 2*V_New[Nx-1][j] - V_New[Nx-2][j];							
							V_New[0][j] = 4*(1+dx)/(2+3*dx) * V_New[1][j] - (2+dx)/(2+3*dx) * V_New[2][j];					
							V_New[Nx][j] = 4*(1-dx)/(2-3*dx) * V_New[Nx-1][j] - (2-dx)/(2-3*dx) * V_New[Nx-2][j];
					
						} //for(j) z풀기


						if(m_curtgreek == 1) { // 만기에
								 
							for(i=0; i<outN + 4; i++) {		
								xtemp = log(SrTemp[i]/(m_xalpha*m_bval[s_idx]));
								NGP = (int)floor((xtemp - X_min)/dx);
								if(NGP<0) NGP = 0;
								if(NGP>=Nx) NGP = Nx-1;

								alpha = (xtemp-X[NGP])/(X[NGP+1]-X[NGP]);
								alphaz = (sum_rate-Z[zcen_idx])/dz;
								SrV1Temp[i] = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) +(alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
							}// for(i)
		 
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
				
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);
			alphaz = (sum_rate - Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) + (alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;
			//value = interp1(slevel,S,newV,SM+1,SmeshType,1);	

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
				
			NGP = (int)floor((cur_x-X_min)/dx);
			if(NGP<0) NGP = 0;
			if(NGP>=Nx) NGP = Nx-1;
			alpha = (cur_x - X[NGP])/(X[NGP+1]-X[NGP]);

			alphaz = (sum_rate-Z[zcen_idx])/dz;
			value = (1-alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx] + alpha*V_New[NGP+1][zcen_idx]) +(alphaz)*((1.0-alpha)*V_New[NGP][zcen_idx+1] + alpha*V_New[NGP+1][zcen_idx+1]) ;

			//value = interp1(slevel,S,newV,SM+1,SmeshType,1); // value - m_Ogreeks[0] 을 이용하기 위하여
		}
	 

		free(X);
		free(Z);

		for(i=0;i<=Nx; i++) {
			free(V_New[i]);
			free(V_Old[i]);
			free(MAu[i]);
			free(MAd[i]);
			free(MAl[i]);
		}	 	 
		free(V_New);
		free(V_Old);
		free(MAu);
		free(MAd);
		free(MAl);
	 
 
		free(volskew);

		return value;
	}
		

} // fd_price();



void CliquetA::curttimeVol(double *vol, int day)
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




double CliquetA::realS(double x, double z, double s0, double z0)
{
	
	double centers;
	double tempz;
	double realsval;


	centers = s0;
	tempz = z-z0;

	if(z>=0) {
		while(tempz>=m_LCap) {
			tempz = tempz - m_LCap;
			centers = centers * (1+m_LCap);
		}
	}
	else {
		while(tempz<=m_LFloor) {
			tempz = tempz - m_LFloor;
			centers = centers * (1+m_LFloor);
		}
	}

	centers = centers * (1+tempz);

	realsval = m_xalpha*centers * exp(x);

	return realsval;
}

  
 