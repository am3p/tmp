#include "StepDownCDC.h"
#include <iostream>
#include <fstream>

extern odysseyRandom *opr;

 

StepDownCDC::StepDownCDC(int kiFlag, int koFlag, double *sval, double *bval, double *xval, double *kihval, double *kohval, double *Cpn,
	      int *matuN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  double *CDlegInfo,
		  int irateN, int *irateBDay,   double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int *corrtype, int *corrN, int *corrBDay, double *corr, long SimulNo)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	
			
	 
	
	
	int i, j;	 
	 
	//입력 데이터 내부 전역 변수로 저장   
	  
	m_matuN = matuN[0];
	m_CDN = matuN[1];
	m_matu = new int[m_matuN];
	m_ObBDay = new int[m_CDN]; //cd
	m_PayBDay = new int[m_CDN]; //cd
	m_CDadjval = bval[AsNUM3];
	m_FlegFactor = new double[m_CDN];
	m_ObCDrate = new double[m_CDN];

	
	for(i=0; i<AsNUM3; i++) {
		m_sval[i] = sval[i];
		m_bval[i] = bval[i];
		m_kihval[i] = kihval[i];
	}

	m_xval = new double *[AsNUM3];
	for(i=0; i<AsNUM3; i++)
		m_xval[i] = new double[m_matuN];
	
	for(i=0; i<AsNUM3; i++) {
		for(j=0; j<m_matuN; j++) {
			m_xval[i][j] = xval[m_matuN*i+j];
		}
	}

	
	m_kiFlag = kiFlag;
	 

	m_Cpn = new double[m_matuN+2];
 

	for(i=0; i<m_matuN; i++) {	 
		m_Cpn[i] = Cpn[i];
		m_matu[i] = matu[i];		 		
	}
	m_Cpn[m_matuN] = Cpn[m_matuN]; // dummy cpn
	m_Cpn[m_matuN+1] = Cpn[m_matuN+1]; // 원금보장 조정 cpn
		 
	/* 
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
	*/

	//cd
	for(i=0; i<m_CDN; i++) {
		m_ObBDay[i] = matu[m_matuN + i];
		m_PayBDay[i] = matu[m_matuN + m_CDN + i];
		m_FlegFactor[i] = CDlegInfo[i];
		m_ObCDrate[i] = CDlegInfo[m_CDN + i];
	}



	
	m_irateN = irateN;
	m_irateBDay = new int[m_irateN];

	m_iRFrate = new double*[AsNUM3];
	for(i=0; i<AsNUM3; i++)
		m_iRFrate[i] = new double[m_irateN];
 
	m_iDCrate = new double[m_irateN];
	m_iIRSrate = new double[m_irateN];

	for(i=0; i<AsNUM3; i++)
		for(j=0; j<m_irateN; j++)
			m_iRFrate[i][j] = iRFrate[irateN*i+j];
	
	for(i=0; i<m_irateN; i++) {
		m_irateBDay[i] = irateBDay[i];
	  
		m_iDCrate[i] = iDCrate[i];
		m_iIRSrate[i] = iIRSrate[i];
	}

	for(i=0; i<AsNUM3; i++) {
		m_divrate[i] = divrate[i];
		m_divbval[i] = divbval[i];
		m_divN[i] = divN[i];
	}

	m_MAXdivN = m_divN[0];
	for(i=1; i<AsNUM3; i++) {
		if(m_MAXdivN < m_divN[i])
			m_MAXdivN = m_divN[i];
	}

	m_divBDay = new int*[AsNUM3];
	m_divCash = new double*[AsNUM3];
	for(i=0; i<AsNUM3; i++) {
		m_divBDay[i] = new int[m_MAXdivN];
		m_divCash[i] = new double[m_MAXdivN];
	}

	int sumdivN;
	sumdivN = 0;
	for(i=0; i<AsNUM3; i++) {		
		for(j=0; j<m_divN[i]; j++) {
			m_divBDay[i][j] = divBDay[sumdivN + j];
			m_divCash[i][j] = divCash[sumdivN + j];
		}
		sumdivN = sumdivN + m_divN[i];
	}
	
	m_divApy = divApy;

	m_voltype = voltype[0];
	for(i=0; i<AsNUM3; i++) {
		m_volbval[i] = volbval[i];
	}
	 
	m_volTN = volTN[0];
	m_volSN = volSN[0];
	m_volBDay = new int[m_volTN];
	for(i=0; i<m_volTN; i++)
		m_volBDay[i] = volBDay[i];
	m_volSmness = new double[m_volSN];
	for(i=0; i<m_volSN; i++)
		m_volSmness[i] = volSmness[i];

	m_vol = new double**[AsNUM3];
	for(i=0; i<AsNUM3; i++) {
		m_vol[i] = new double*[m_volSN];
		for(j=0; j<m_volSN; j++) {
			m_vol[i][j] = new double[m_volTN];
		}
	}

	int k;

	for(k=0; k<AsNUM3; k++) {
		for(i=0; i<m_volSN; i++) {
			for(j=0; j<m_volTN; j++) {
				m_vol[k][i][j] = vol[m_volSN*m_volTN*k + m_volTN*i + j];
			}
		}
	}

	 

	m_corrtype = corrtype[0];
	m_corrN = corrN[0];
	m_corrBDay = new int[m_corrN];
	m_corr = new double**[m_corrN];
	for(i=0; i<m_corrN; i++) {
		m_corr[i] = new double*[AsNUM3];
		for(j=0; j<AsNUM3; j++) {		
			m_corr[i][j] = new double[AsNUM3];
		}
	}
	
	for(i=0; i<m_corrN; i++) {
		m_corrBDay[i] = corrBDay[i];
		for(j=0; j<AsNUM3; j++) {
			for(k=0; k<AsNUM3; k++) {
				m_corr[i][j][k] = corr[AsNUM3*AsNUM3*i + AsNUM3*j + k];
			}
		}
	}
    

	// QT
	
	m_FXvoltype = voltype[1];
	for(i=0; i<AsNUM3; i++)
		m_FXvolbval[i] = volbval[AsNUM3 + i];
	 
	m_FXvolTN = volTN[1];
	m_FXvolSN = volSN[1];
	m_FXvolBDay = new int[m_FXvolTN];
	for(i=0; i<m_FXvolTN; i++)
		m_FXvolBDay[i] = volBDay[m_volTN+i];
	m_FXvolSmness = new double[m_FXvolSN];
	for(i=0; i<m_FXvolSN; i++)
		m_FXvolSmness[i] = volSmness[m_volSN+i];

	m_FXvol = new double**[AsNUM3];
	for(k=0; k<AsNUM3; k++) {
		m_FXvol[k] = new double*[m_FXvolSN];
		for(i=0; i<m_FXvolSN; i++) {
			m_FXvol[k][i] = new double[m_FXvolTN];
		}
	}

	for(k=0; k<AsNUM3; k++) {
		for(i=0; i<m_FXvolSN; i++) {
			for(j=0; j<m_FXvolTN; j++) {
				m_FXvol[k][i][j] = vol[AsNUM3*m_volSN*m_volTN + k*m_FXvolSN*m_FXvolTN + i*m_FXvolTN + j];
			}
		}
	}

	 

	m_FXcorrtype = corrtype[1];
	m_FXcorrN = corrN[1];
	m_FXcorrBDay = new int[m_FXcorrN];
	m_FXcorr = new double*[AsNUM3];

		
	for(i=0; i<m_FXcorrN; i++)  
		m_FXcorrBDay[i] = corrBDay[m_corrN+i];

	for(i=0; i<AsNUM3; i++) 
		m_FXcorr[i] = new double[m_FXcorrN];

	for(i=0; i<AsNUM3; i++) {
		for(j=0; j<m_FXcorrN; j++) {		
			m_FXcorr[i][j] = corr[m_corrN*AsNUM3*AsNUM3 + m_FXcorrN*i + j];
		}
	}

	 
     
	 
	//m_endS = m_xval[m_matuN-2]; // 마지막 조기상환일에 필요 조기상환일에 세팅!! 

	 

	// ki 체크
	for(i=0; i<AsNUM3; i++) {
		if(m_sval[i] <= m_kihval[i]) {
			m_kiFlag = 1; // 하나라도 치면..
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

StepDownCDC::~StepDownCDC()
{	
	
	int i,k;	 
	//delete m_evalN;
	delete m_matu;
	delete m_ObBDay;//cd
	delete m_PayBDay; //cd
	delete m_FlegFactor;
	delete m_ObCDrate;
	delete m_Cpn;
	delete m_irateBDay;	
	delete m_iDCrate;  
	delete m_iIRSrate; 	 			
	delete m_volBDay;	 
	delete m_volSmness;	 	 
	delete m_corrBDay;	
	delete m_FXvolBDay;
	delete m_FXvolSmness;
	delete m_FXcorrBDay;
	 
 

	//
	for(k=0; k<AsNUM3; k++) {
		for(i=0; i<m_volSN; i++) {
			delete m_vol[k][i];
		}
		for(i=0; i<m_FXvolSN; i++) {
			delete m_FXvol[k][i];
		}
	}

	for(k=0; k<AsNUM3; k++) {
		delete m_vol[k];
		delete m_FXvol[k];

		delete m_xval[k];
		delete m_iRFrate[k];
		delete m_divBDay[k];
		delete m_divCash[k];
		delete m_FXcorr[k];
	}
	
	 
	
	delete m_vol;
	delete m_FXvol;

	delete m_xval;
	delete m_iRFrate;
	delete m_divBDay;
	delete m_divCash;
	delete m_FXcorr;
	//
	
	 
	

	for(k=0; k<m_corrN; k++) {
		for(i=0; i<AsNUM3; i++) {
			delete m_corr[k][i];
		}		
	}	
	for(k=0; k<m_corrN; k++) {
		delete m_corr[k];
	}
	delete m_corr;
	
	rnd->randseq(2,m_noiters,MAXDAYS,NULL);
	rnd->randnum(2,m_noiters,MAXAST,NULL,7);
	delete rnd;
}

double StepDownCDC::getValue() { 
	
	double value;	 
	
	value = mc_price(); 
		
	m_Ogreeks[S3Priceidx] = value;
	 
  

	return value;

}



 

double StepDownCDC::getDeltaGammaTheta() {
 
	double tempvalue[AsNUM3][2];
	double value;
	int i,k;

	double tempsval[AsNUM3];

	for(i=0; i<AsNUM3; i++)
		tempsval[i] = m_sval[i];

	for(i=0; i<AsNUM3; i++) {
		m_sval[i] = tempsval[i]*(1+0.01);
		value = mc_price();
		tempvalue[i][0] = value;

		m_sval[i] = tempsval[i]*(1-0.01);
		value = mc_price();
		tempvalue[i][1] = value;

		m_sval[i] = tempsval[i];
	}
	
	m_Ogreeks[S3Delta1idx] = (tempvalue[0][0] - tempvalue[0][1])/2;
	m_Ogreeks[S3Delta2idx] = (tempvalue[1][0] - tempvalue[1][1])/2;
	m_Ogreeks[S3Delta3idx] = (tempvalue[2][0] - tempvalue[2][1])/2;
	
	m_Ogreeks[S3Gamma1idx] = tempvalue[0][0] - 2*m_Ogreeks[S3Priceidx] + tempvalue[0][1];
	m_Ogreeks[S3Gamma2idx] = tempvalue[1][0] - 2*m_Ogreeks[S3Priceidx] + tempvalue[1][1];
	m_Ogreeks[S3Gamma3idx] = tempvalue[2][0] - 2*m_Ogreeks[S3Priceidx] + tempvalue[2][1];

	 
	if(m_matu[m_matuN-1] <= 0) {		
		value = m_Ogreeks[S3Priceidx];
	}
	else {
	
		for(i=0; i<m_matuN; i++)		
			m_matu[i] = m_matu[i] - 1;
		for(k=0; k<AsNUM3; k++)
			for(i=1;i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] - 1;
			
		for(i=0; i<m_CDN; i++) { //cd
			m_ObBDay[i] = m_ObBDay[i] - 1;				
			m_PayBDay[i] = m_PayBDay[i] - 1;
		}

		value = mc_price();

		for(i=0; i<m_matuN; i++)
			m_matu[i] = m_matu[i] +1;
		for(k=0; k<AsNUM3; k++)
			for(i=1; i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] + 1;
				
		for(i=0; i<m_CDN; i++) {
			m_ObBDay[i] = m_ObBDay[i] + 1;				
			m_PayBDay[i] = m_PayBDay[i] + 1;
		}
	}
	
	for(i=0; i<m_matuN-1; i++) {
		if(m_matu[i] == 0) {
			int CntAs;
			CntAs = AsNUM3;
			for(k=0; k<AsNUM3; k++) {
				if(m_sval[k] < m_xval[k][i]) {
					CntAs--;
					break;
				}
			}
			if(CntAs == AsNUM3) 
				value = m_Ogreeks[S3Priceidx];

			break;
		}
	}

	m_Ogreeks[S3Thetaidx] = value - m_Ogreeks[S3Priceidx];	 
 

	return value;

}

double StepDownCDC::getVega(int vegaN, int *vegaidx) {
		 
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
	

	m_Ogreeks[S3Vega1idx] = 0;
	 

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
	 

		m_Ogreeks[S3Vega1idx] += m_Otvegas[i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[0][i][j] = tempv[i][j]; // vol 복원

	/////////////////////v1 end

	
	// v2 시작
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[1][i][j]; // vol 저장
	

	m_Ogreeks[S3Vega2idx] = 0;
	 

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
						m_vol[1][k][j] = tempv[k][j] + curtpertubation;
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
							m_vol[1][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[1][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[1][k][j] = tempv[k][j] + curtpertubation;
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
		m_Otvegas[vegaN+i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
	 

		m_Ogreeks[S3Vega2idx] += m_Otvegas[vegaN+i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[1][i][j] = tempv[i][j]; // vol 복원

	/////////////////////v2 end

	
	// v3 시작
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[2][i][j]; // vol 저장
	

	m_Ogreeks[S3Vega3idx] = 0;
	 

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
						m_vol[2][k][j] = tempv[k][j] + curtpertubation;
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
							m_vol[2][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[2][k][j] = tempv[k][j] + curtpertubation;
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
								m_vol[2][k][j] = tempv[k][j] + curtpertubation;
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
		m_Otvegas[2*vegaN+i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
	 

		m_Ogreeks[S3Vega3idx] += m_Otvegas[2*vegaN+i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[2][i][j] = tempv[i][j]; // vol 복원

	/////////////////////v3 end

	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	 
	return m_Ogreeks[twoVegaidx1];

}


double StepDownCDC::getRho(int rhoN,int *rhoidx) {
	
	 
	double value;

	double *tempirate;  
	int i,j,k,subk;

	int pertubationX[3];
	double pertubationY[3];		  

	int same_q[AsNUM3+1];

	tempirate = new double[m_irateN];
	 
	for(i=0; i<=AsNUM3; i++)
		same_q[i] = -1;
	//same_q[AsNUM3] = 1; //이줄 주석 처리하면 CD Sawp 처리 가능



	m_Ogreeks[S3Rhorf1idx] = 0; //rf rho sum
	m_Ogreeks[S3Rhorf2idx] = 0; //rf rho sum
	m_Ogreeks[S3Rhorf3idx] = 0; //rf rho sum
	m_Ogreeks[S3Rhodcidx] = 0; //rf rho sum
	m_Ogreeks[S3RhoCDidx] = 0; //rf rho sum

	for(i=0; i<(AsNUM3+1)*MAXRHOS; i++)
		m_Otrfrhos[i] = 0;

	//rf1 rho 시작

	for(k=0; k<AsNUM3; k++) {
		if(same_q[k] != -1)
			continue;

		for(i=0; i<m_irateN; i++) 
			tempirate[i] = m_iRFrate[0][i];

		for(subk=k+1; subk<AsNUM3; subk++)
			if(m_iRateSeq[k] == m_iRateSeq[subk])
				same_q[subk] = 0;

		if(m_iRateSeq[k] == m_iRateSeq[AsNUM3] && same_q[AsNUM3] == -1)
			same_q[AsNUM3] = 0;
		 
		for(i=0; i<rhoN; i++) {
			if(rhoN == 1) { // 평행이동 rho
				for(j=0; j<m_irateN; j++) {
					m_iRFrate[k][j] = tempirate[j] + RhoPertubation;	
					for(subk=k+1; subk<AsNUM3; subk++)
						if(same_q[subk] == 0) 					
							m_iRFrate[subk][j] = tempirate[j] + RhoPertubation;	
					if(same_q[AsNUM3] == 0)
						m_iIRSrate[j] = tempirate[j] + RhoPertubation;
				}				
			}
			else {
				if(i == 0) {		
					pertubationX[0] = m_irateBDay[rhoidx[i]];
					pertubationX[1] = m_irateBDay[rhoidx[i+1]];

					pertubationY[0] = RhoPertubation;
					pertubationY[1] = 0.0;

					for(j=0; j<m_irateN; j++) {			
						m_iRFrate[k][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
						for(subk=k+1; subk<AsNUM3; subk++)
							if(same_q[subk] == 0)
								m_iRFrate[subk][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
						if(same_q[AsNUM3] == 0)
							m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
					}
				}
				else {
					if(i == rhoN-1) {					
						pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
						pertubationX[1] = m_irateBDay[rhoidx[i]];			

						pertubationY[0] = 0.0;			
						pertubationY[1] = RhoPertubation;			
			
						for(j=0; j<m_irateN; j++) {							
							m_iRFrate[k][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
							for(subk=k+1; subk<AsNUM3; subk++)
								if(same_q[subk] == 0)
									m_iRFrate[subk][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
							if(same_q[AsNUM3] == 0)
								m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
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
							m_iRFrate[k][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat
							for(subk=k+1; subk<AsNUM3; subk++)
								if(same_q[subk] == 0)
									m_iRFrate[subk][j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat
							if(same_q[AsNUM3] == 0)
								m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat							
						}	

					}
				}
			}// if(rhoN == 0) else

			value = mc_price();

			m_Otrfrhos[rhoN*k+i] = value - m_Ogreeks[S3Priceidx]; // 10bp term rf rho 
			m_Ogreeks[S3Rhorf1idx + k] += m_Otrfrhos[rhoN*k+i];			 
		}
 

		for(i=0; i<m_irateN; i++) {
			m_iRFrate[k][i] = tempirate[i]; // RF rate 복원						 
		}

		for(subk=k+1;subk<AsNUM3; subk++)
			if(same_q[subk] == 0) {
				for(i=0; i<m_irateN; i++)
					m_iRFrate[k][i] = tempirate[i];
				same_q[subk] = 1;
			}
		
		if(same_q[AsNUM3] == 0) {
			for(i=0; i<m_irateN; i++)
				m_iIRSrate[i] = tempirate[i];
			same_q[AsNUM3] = 1;
		}
	}
 

//////////////////////////////
	//dc rho 시작
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iDCrate[i]; // DC rate 저장
	}
	 
	m_Ogreeks[S3Rhodcidx] = 0; // dc rho sum
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

		m_Otdcrhos[i] = value - m_Ogreeks[S3Priceidx]; // 10bp term rf rho 
		m_Ogreeks[S3Rhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate 복원
	}
 	
//////////////////////////////
	 
	//cd rho 시작 // CD swap 만
	if(same_q[AsNUM3] == -1) {
		same_q[AsNUM3] = 0;

		for(i=0; i<m_irateN; i++) {
			tempirate[i] = m_iIRSrate[i]; // DC rate 저장
		}
	 
		m_Ogreeks[S3RhoCDidx] = 0; // dc rho sum
		for(i=0; i<rhoN; i++) {
			if(rhoN == 1) { // 평행이동 rho
				for(j=0; j<m_irateN; j++) 
					m_iIRSrate[j] = tempirate[j] + RhoPertubation;
			}
			else {
				if(i == 0) {		
					pertubationX[0] = m_irateBDay[rhoidx[i]];
					pertubationX[1] = m_irateBDay[rhoidx[i+1]];

					pertubationY[0] = RhoPertubation;
					pertubationY[1] = 0.0;

					for(j=0; j<m_irateN; j++) {			
						m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat
					}
				}
				else {
					if(i == rhoN-1) {					
						pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
						pertubationX[1] = m_irateBDay[rhoidx[i]];			

						pertubationY[0] = 0.0;			
						pertubationY[1] = RhoPertubation;			
			
						for(j=0; j<m_irateN; j++) {							
							m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2개, 1:비균등, 0:양쪽 flat				
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
							m_iIRSrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3개, 1:비균등, 0:양쪽 flat				
						}	

					}
				}
			} // if(rhoN == 0) else

			value = mc_price();

			m_Otrfrhos[rhoN*AsNUM3+i] = value - m_Ogreeks[S3Priceidx]; // 10bp term rf rho 
			m_Ogreeks[S3RhoCDidx] += m_Otrfrhos[rhoN*AsNUM3+i];
		}
 

			
		for(i=0; i<m_irateN; i++) {
			m_iIRSrate[i] = tempirate[i]; // DC rate 복원
		}
		same_q[AsNUM3] = 1;
	}
	


	delete tempirate;
 

		
	return m_Ogreeks[S3Rhorf1idx] + m_Ogreeks[S3Rhorf2idx] + m_Ogreeks[S3Rhorf3idx] + m_Ogreeks[S3Rhodcidx] + m_Ogreeks[S3RhoCDidx];  // 의미 없음 

}


 
 

double StepDownCDC::mc_price()
{
	long *rsequns;
	double **randoms;
	long nostep;

	double simsval[AsNUM3];
	double ***Cholcor;
	
	int simkiFlag;
	double simprice, simcdprice;
	double price, cdprice;
	double dT, sqdT;
		
	int i,j,k;

	dT = 1/YFactor;
	sqdT = sqrt(dT);
	double mudt[AsNUM3];
	double vldt[AsNUM3];
	//m_modeltype = 0;

	double disc_div[AsNUM3];
	double Srv, FXSrv;
	double *volskew, *FXvolskew;

	double r[AsNUM3], q[AsNUM3], v[AsNUM3];
	double FXv[AsNUM3], FXcor[AsNUM3];

	double **dailyIFrf, *dailyZdc;

	//cd
	double *fdObCDrate, *FlegCpn, *CpnCDleg;
 
	 
	

	nostep = MAX(m_matu[m_matuN-1],1);
	 

	rsequns = (long*)malloc((size_t)((nostep)*sizeof(long)));
	randoms =(double **)malloc((size_t)((AsNUM3)*sizeof(double *)));	
	if(!rsequns || !randoms) return -100;

	for ( i = 0; i < AsNUM3; i++ ) {
		randoms[i] = (double *)malloc((size_t)((m_noiters)*sizeof(double)));		 
		if(!randoms[i]) {
			for(j=0; j<i; j++)
				free(randoms[j]);
			return -200;
		}
	}

	if ( rnd->randnum(1,m_noiters,AsNUM3,randoms,0) != 0 ) {
		free(rsequns);
		for(i=0; i<AsNUM3; i++)
			free(randoms[i]);
		free(randoms);
		return -600;
	}
	   
	dailyIFrf = new double*[AsNUM3];
	for(k=0; k<AsNUM3; k++) {
		dailyIFrf[k] = new double[nostep];		
	}
	dailyZdc = new double[nostep+1];

	for(k=0; k<AsNUM3; k++) {
		for(i=0; i<nostep; i++) { //사용할 때 i=1:nostep i 대응 [i-1]로
			dailyIFrf[k][i] = interp1(i+1,m_irateBDay,m_iRFrate[k],m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate[k],m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
		}
	}

	for(i=0; i<=nostep; i++) { //사용할때 i 대응 [i]로 사용
		dailyZdc[i] = interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0);
	}

	
	Cholcor = new double**[m_corrN];
	for(i=0; i<m_corrN; i++) {
		Cholcor[i] = new double*[AsNUM3];
		for(j=0; j<AsNUM3; j++) {
			Cholcor[i][j] = new double[AsNUM3];
			for(k=0; k<AsNUM3; k++)
				Cholcor[i][j][k] = 0.0;		
		}
	}

	volskew  = new double[m_volSN];
	FXvolskew  = new double[m_FXvolSN];
	 
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
 

	int coi, coj, cok;
	double tempsum;
	for(i=0; i<m_corrN; i++) {
		for(cok=0; cok<AsNUM3; cok++) {
			tempsum = 0;
			for(coi=0; coi<cok; coi++)
				tempsum = tempsum + Cholcor[i][cok][coi]*Cholcor[i][cok][coi];

			Cholcor[i][cok][cok] = sqrt(m_corr[i][cok][cok]-tempsum);
						
			for(coi=cok+1; coi<AsNUM3; coi++) {
				tempsum = 0;
				for(coj=0; coj<cok; coj++)
					tempsum = tempsum + Cholcor[i][coi][coj]*Cholcor[i][cok][coj];

				Cholcor[i][coi][cok] = (m_corr[i][coi][cok] - tempsum)/Cholcor[i][cok][cok];
			}
		}
	}//Chol 분해!!

	for(i=0; i<AsNUM3; i++)
		q[i] = m_divrate[i];

	int si,l,sk,disci,cddayi;
	int CntAs;
	double crand;
	int curtCi, crandi;
	double Wps;
	price = 0;
	cdprice = 0;

	if(m_modeltype == 10 || m_xval[0][m_matuN-1] == m_kihval[0])
		m_kiFlag = 1;
 
	
	for(si=0; si<m_noiters; si++) { //MC SIM
		simkiFlag = m_kiFlag;
		simprice = 0;
		simcdprice = 0;

		if(rnd->randseq(1,si,nostep,rsequns) != 0) {
			for(i=0; i<AsNUM3; i++) {
				free(randoms[i]);
				delete dailyIFrf[i];				 
			}
			free(randoms);
			delete dailyIFrf;
			free(rsequns);
			delete dailyZdc;
			delete volskew;
			delete FXvolskew;

			for(i=0; i<m_corrN; i++) {
				for(j=0; j<AsNUM3; j++) {
					delete Cholcor[i][j];
				}
				delete Cholcor[i];
			}
			delete Cholcor;

	
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
					for(k=0; k<AsNUM3; k++)					
						simsval[k] = m_sval[k];
				}
				else { //if l
					for(k=0; k<AsNUM3; k++) {
						r[k] = dailyIFrf[k][l-1];

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

						//FXvol
						if(m_FXvoltype == 0) {
							FXvolskew[0] = m_FXvol[k][0][0];
						}
						else {
							curttimeFXVol(FXvolskew, l, k);
						}
						if(m_FXvoltype == 2) {
							FXSrv = m_FXvolbval[k]/m_FXvolbval[k];  // ATM만.. 사용
							FXv[k] = interp1(FXSrv, m_FXvolSmness, FXvolskew, m_FXvolSN,1,0);
						}
						else {
							FXv[k] = FXvolskew[0];
						}
						 
						//FXcor 
						FXcor[k] = curttimeFXCorr(l, k);

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

						mudt[k] = (r[k]-q[k]-0.5*v[k]*v[k] - FXcor[k]*FXv[k]*v[k])*dT;
						vldt[k] = v[k]*sqdT;
					}

					//cor
					curtCi = curttimeCorridx(l);

					for(k=0; k<AsNUM3; k++) {		
						crand = 0;
						for(crandi=0; crandi<AsNUM3; crandi++) 
							crand = crand + Cholcor[curtCi][k][crandi]*randoms[crandi][rsequns[l-1]];

						simsval[k] = simsval[k] * exp(mudt[k] + vldt[k]*crand);

						if(simkiFlag == 0)
							if(simsval[k] <= m_kihval[k])
							simkiFlag = 1;
					}

				}//if l==0 else

				for(cddayi=0; cddayi<m_CDN-1; cddayi++) { // 만기는 따로 처리 하기 위해서
					if(l == m_PayBDay[cddayi]) {
						simcdprice = simcdprice + AMOUNT * (FlegCpn[cddayi]) * exp(-dailyZdc[l]*l/YFactor);
						break;
					}
				}


			}//for l

			CntAs = AsNUM3;
			for(k=0; k<AsNUM3; k++) {
				if(simsval[k] < m_xval[k][sk]) {
					CntAs--;
					break;
				}				
			}

			if(CntAs == AsNUM3) {
				simprice = AMOUNT * (m_Cpn[sk])*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				simcdprice = simcdprice + AMOUNT * CpnCDleg[sk] * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				break;
			}

			if(sk==m_matuN-1) {
				if(simkiFlag == 0)
					simprice = AMOUNT *(m_Cpn[sk+1])*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				else {
					Wps = simsval[0]/m_bval[0];
					for(k=1; k<AsNUM3; k++) 
						Wps = MIN(Wps,simsval[k]/m_bval[k]);
					simprice = AMOUNT * (-1 + MIN(Wps,1)) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				}

				simcdprice = simcdprice + AMOUNT * CpnCDleg[sk] * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
			}
				
		} // for sk

		price = price + simprice/m_noiters;
		cdprice = cdprice + simcdprice/m_noiters;

	} // MC SIM
 
	
		 
		
	//메모리 해제
	for(i=0; i<AsNUM3; i++) {
		free(randoms[i]);
		delete dailyIFrf[i];	 
	}
	free(randoms);
	delete dailyIFrf;
	free(rsequns);
	delete dailyZdc;
	delete volskew;
	delete FXvolskew;

	for(i=0; i<m_corrN; i++) {
		for(j=0; j<AsNUM3; j++) {
			delete Cholcor[i][j];
		}
		delete Cholcor[i];
	}
	delete Cholcor;
			
	delete fdObCDrate;
	delete FlegCpn;
	delete CpnCDleg; 

	// 메모리 해제
 
			 
		
	return price-cdprice;
	 		



} //mc_price();



void StepDownCDC::curttimeVol(double *vol, int day, int Asidx)
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
 


int StepDownCDC::curttimeCorridx(int day)
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

	return findidx;
	 
	
}
 


void StepDownCDC::curttimeFXVol(double *vol, int day, int Asidx)
{
	int i;
	
	if(m_FXvoltype == 1) { // term imp vol
		double v1, v2;
		int t1, t2;
		double *termvol;

		termvol = new double[m_FXvolTN];
		for(i=0; i<m_FXvolTN; i++) 
			termvol[i] = m_FXvol[Asidx][0][i];
		
		t1 = day;
		v1 = interp1(t1,m_FXvolBDay,termvol,m_FXvolTN,1,0);
		t2 = day+1;
		v2 = interp1(t2,m_FXvolBDay,termvol,m_FXvolTN,1,0);

		vol[0] = sqrt(MAX(FXminVol*FXminVol, (v2*v2*t2 - v1*v1*t1)/(t2-t1))); // implied forward vol

		delete termvol;
	}
	else { // step function 형 local vol
		int findidx;
		
		findidx = m_FXvolTN-1;
		for(i=0; i<m_FXvolTN; i++) {
			if(day <= m_FXvolBDay[i]) {
				findidx = i;
				break;
			}
		}		 

		for(i=0; i<m_FXvolSN; i++) 
			vol[i] = m_FXvol[Asidx][i][findidx];
	}
	
}
 
double StepDownCDC::curttimeFXCorr(int day, int Asidx)
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

	return m_FXcorr[Asidx][findidx];
	 
	
}


/*
double StepDownCDC::curttimeFXCorr2(int day)
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




int StepDownCDC::findinteridx(double x, double *Dx, int DxN)
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

 
 