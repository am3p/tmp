#include "AirbagStepDownB.h"
#include <iostream>
#include <fstream>

extern odysseyRandom *opr;

AirbagStepDownB::AirbagStepDownB(int kiFlag, int koFlag, double *sval, double *bval, double *xval, double *kihval, double *kohval, double *airbagxval, double *Cpn,
	      int matuN, int airbagN, int *matu, double ctime,
		  int maxevalN, int *evalN, double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int corrtype, int corrN, int *corrBDay, double *corr, long SimulNo)
		  // shiftRate는 배리어 sift overhedge로 배리어를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 배리어 spread overhedge로 배리어를 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	int i, j;	 
	 
	//입력 데이터 내부 전역 변수로 저장   
	  
	m_matuN = matuN;
	m_matu = new int[m_matuN];

	for(i=0; i<AsNUM; i++) {
		m_sval[i] = sval[i];
		m_bval[i] = bval[i];
		m_kihval[i] = kihval[i];
		m_airbagxval[i] = airbagxval[i];
	}

	m_xval = new double*[AsNUM];
	for(i=0; i<AsNUM; i++)
		m_xval[i] = new double[m_matuN];
	
	for(i=0; i<AsNUM; i++) {
		for(j=0; j<m_matuN; j++) {
			m_xval[i][j] = xval[m_matuN*i+j];
		}
	}	 
	 
	m_kiFlag = kiFlag;	 

	m_airbagN = airbagN;

	m_Cpn = new double[m_matuN+2];	 

	for(i=0; i<m_matuN; i++) {	 
		m_Cpn[i] = Cpn[i];
		m_matu[i] = matu[i];
	 	
	}
	m_Cpn[m_matuN] = Cpn[m_matuN]; // dummy cpn
	m_Cpn[m_matuN+1] = Cpn[m_matuN+1]; // 원금보장 조정 cpn
		 
	m_ctime = ctime;

	m_maxevalN = maxevalN;
	m_evalN = new int[m_matuN];
	for(i=0; i<m_matuN; i++) 
		m_evalN[i] = evalN[i];

	m_psval = new double*[AsNUM];
	for(i=0; i<AsNUM; i++)
		m_psval[i] = new double[m_maxevalN];
	for(i=0; i<AsNUM; i++) {
		for(j=0; j<m_maxevalN; j++) {
			m_psval[i][j] = psval[m_maxevalN*i + j];
		}
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


	for(i=0; i<AsNUM; i++) {
		m_divrate[i] = divrate[i];
		m_divbval[i] = divbval[i];
		m_divN[i] = divN[i];
	}

	m_MAXdivN = m_divN[0];
	for(i=1; i<AsNUM; i++) {
		if(m_MAXdivN < m_divN[i])
			m_MAXdivN = m_divN[i];
	}

	m_divBDay = new int*[AsNUM];
	m_divCash = new double*[AsNUM];
	for(i=0; i<AsNUM; i++) {
		m_divBDay[i] = new int[m_MAXdivN];
		m_divCash[i] = new double[m_MAXdivN];
	}

	int sumdivN;
	sumdivN = 0;
	for(i=0; i<AsNUM; i++) {		
		for(j=0; j<m_divN[i]; j++) {
			m_divBDay[i][j] = divBDay[sumdivN + j];
			m_divCash[i][j] = divCash[sumdivN + j];
		}
		sumdivN = sumdivN + m_divN[i];
	}
	
	m_divApy = divApy;

	m_voltype = voltype;
	for(i=0; i<AsNUM; i++) {
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

	m_vol = new double**[AsNUM];
	for(i=0; i<AsNUM; i++) {
		m_vol[i] = new double*[m_volSN];
		for(j=0; j<m_volSN; j++) {
			m_vol[i][j] = new double[m_volTN];
		}
	}

	int k;

	for(k=0; k<AsNUM; k++) {
		for(i=0; i<m_volSN; i++) {
			for(j=0; j<m_volTN; j++) {
				m_vol[k][i][j] = vol[m_volSN*m_volTN*k + m_volTN*i + j];
			}
		}
	}

	m_corrtype = corrtype;
	m_corrN = corrN;
	m_corrBDay = new int[m_corrN];
	m_corr = new double[m_corrN];
	for(i=0; i<m_corrN; i++) {
		m_corrBDay[i] = corrBDay[i];
		m_corr[i] = corr[i];
	}
     

	 
	//m_endS = m_xval[m_matuN-2]; // 마지막 조기상환일에 필요 조기상환일에 세팅!! 

	 

	// ki 체크
	for(i=0; i<AsNUM; i++) {
		if(m_sval[i] <= m_kihval[i]) {
			m_kiFlag = 1;
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

AirbagStepDownB::~AirbagStepDownB()
{	
	int i,k;	 
	delete m_evalN;
 
	delete m_irateBDay;
 
	delete m_iRFrate; 
	delete m_iDCrate;  
	delete m_iIRSrate; 	 
			 
	delete m_volBDay;	 
	delete m_volSmness;	 
 	 
	delete m_corrBDay;
	delete m_corr;
   
	delete m_matu;
 
	delete m_Cpn;
 
	
	for(k=0; k<AsNUM; k++) {
		for(i=0; i<m_volSN; i++) {
			delete m_vol[k][i];
		}
	}

	for(i=0; i<AsNUM; i++) {
		delete m_vol[i];

		delete m_xval[i];
		delete m_psval[i];
		delete m_divBDay[i];
		delete m_divCash[i];
	}

	delete m_vol;

	delete m_xval;
	delete m_psval;
	delete m_divBDay;
	delete m_divCash;

	rnd->randseq(2,m_noiters,MAXDAYS,NULL);
	rnd->randnum(2,m_noiters,MAXAST,NULL,7);
	delete rnd;
}

double AirbagStepDownB::getValue() {
	 
	double value;
	 
	m_curtgreek = 1;
	value = mc_price();
			
	m_Ogreeks[twoPriceidx] = value; 

	return value;

}

double AirbagStepDownB::getDeltaGammaTheta() {
 
	double tempvalue[AsNUM][2];
	double value;
	int i,k;
	double temppsval[AsNUM][10]; // 최대 10영업일 평균으로 ..

	double tempsval[AsNUM];

	for(i=0; i<AsNUM; i++)
		tempsval[i] = m_sval[i];

	for(i=0; i<AsNUM; i++) {
		m_sval[i] = tempsval[i]*(1+0.01);
		value = mc_price();
		tempvalue[i][0] = value;

		m_sval[i] = tempsval[i]*(1-0.01);
		value = mc_price();
		tempvalue[i][1] = value;

		m_sval[i] = tempsval[i];
	}
	
	m_Ogreeks[twoUpDeltaCidx1] = tempvalue[0][0] - m_Ogreeks[twoPriceidx];
	m_Ogreeks[twoDownDeltaCidx1] = m_Ogreeks[twoPriceidx] - tempvalue[0][1];

	m_Ogreeks[twoUpDeltaCidx2] = tempvalue[1][0] - m_Ogreeks[twoPriceidx];
	m_Ogreeks[twoDownDeltaCidx2] = m_Ogreeks[twoPriceidx] - tempvalue[1][1];
	 
	
	m_Ogreeks[twoUpGammaCidx1] = tempvalue[0][0] - 2*m_Ogreeks[twoPriceidx] + tempvalue[0][1];
	m_Ogreeks[twoDownGammaCidx1] = tempvalue[0][0] - 2*m_Ogreeks[twoPriceidx] + tempvalue[0][1];

	m_Ogreeks[twoUpGammaCidx2] = tempvalue[1][0] - 2*m_Ogreeks[twoPriceidx] + tempvalue[1][1];
	m_Ogreeks[twoDownGammaCidx2] = tempvalue[1][0] - 2*m_Ogreeks[twoPriceidx] + tempvalue[1][1];
 
	int tempairbagN;
	tempairbagN = m_airbagN;
	 
	if(m_matu[m_matuN-1] <= 0) {		
		value = m_Ogreeks[twoPriceidx];
	}
	else {

		
		for(i=0; i<m_matuN-1; i++) {
			if(m_matu[i] == 0) {
				int CntAs;
				double tempave;

				CntAs = AsNUM;
			

				for(k=0; k<AsNUM; k++) {

					tempave = 0;
					for(int j=0; j<m_evalN[i]; j++) {
						tempave = tempave + m_psval[k][j];
					}
					tempave = tempave/m_evalN[i];

					if(tempave < m_airbagxval[k]) {
						CntAs--;
						break;
					}
				}
				if(CntAs == AsNUM) 
					m_airbagN++;

				break;
			}
		}


		for(k=0; k<AsNUM; k++) {
			for(i=0; i<m_maxevalN; i++) {
				temppsval[k][i] = m_psval[k][i];
			}
		}
		for(k=0; k<AsNUM; k++) {
			m_psval[k][0] = m_sval[k];
			for(i=1; i<m_maxevalN; i++) {
				m_psval[k][i] = temppsval[k][i-1];
			}
		}
	
		for(i=0; i<m_matuN; i++)		
			m_matu[i] = m_matu[i] - 1;
		for(k=0; k<AsNUM; k++)
			for(i=1;i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] - 1;

		value = mc_price();
		
		for(k=0; k<AsNUM; k++) {
			for(i=0; i<m_maxevalN; i++) {
				m_psval[k][i] = temppsval[k][i];
			}
		}

		for(i=0; i<m_matuN; i++)
			m_matu[i] = m_matu[i] +1;
		for(k=0; k<AsNUM; k++)
			for(i=1; i<m_divN[k]; i++)
				m_divBDay[k][i] = m_divBDay[k][i] + 1;

	}

	m_airbagN = tempairbagN;
	
	for(i=0; i<m_matuN-1; i++) {
		if(m_matu[i] == 0) {
			int CntAs;
			double tempave;

			CntAs = AsNUM;
			

			for(k=0; k<AsNUM; k++) {

				tempave = 0;
				for(int j=0; j<m_evalN[i]; j++) {
					tempave = tempave + m_psval[k][j];
				}
				tempave = tempave/m_evalN[i];

				if(tempave < m_xval[k][i]) {
					CntAs--;
					break;
				}
			}
			if(CntAs == AsNUM) 
				value = m_Ogreeks[twoPriceidx];

			break;
		}
	}

	m_Ogreeks[twoThetaidx] = value - m_Ogreeks[twoPriceidx];	 
 

	return value;

}


double AirbagStepDownB::getVega(int vegaN, int *vegaidx) {
	 
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
			tempv[i][j] = m_vol[0][i][j]; // vol 저장
	

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
				Price[j][udi] =value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+
		} //for(vi)
		//
		/*
		m_Otvegas1[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas1[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas1[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/
		
		//2012-05-23 1% vanna cash / 1% zomma cash
		m_Otvegas1[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		 

		m_Ogreeks[twoVegaidx1] += m_Otvegas1[i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[0][i][j] = tempv[i][j]; // vol 복원

	/////////////////////v1 end

		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[1][i][j]; // vol 저장
	

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
				Price[j][udi] = value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+
		} //for(vi)
		//
		/*
		m_Otvegas2[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
		m_Otvannas2[i] = (Price[2][0] - Price[2][1] - Price[0][0] + Price[0][1])/(4*0.01); // term vanna cash
		m_Otzommas2[i] = (Price[2][0] - 2*Price[1][0] + Price[0][0] - Price[2][1] + 2*Price[1][1] - Price[0][1])/(2*0.01); // term zomma cash
		*/

		//2012-05-23 1% vanna cash / 1% zomma cash
		m_Otvegas2[i] = (Price[1][0] - Price[1][1])/2.0;	 // term vega
	 

		m_Ogreeks[twoVegaidx2] += m_Otvegas2[i];
	 
	}//for(i)
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[1][i][j] = tempv[i][j]; // vol 복원
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	 
	return m_Ogreeks[twoVegaidx1];

}


double AirbagStepDownB::getRho(int rhoN,int *rhoidx) {
	
 
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

		value = mc_price();

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

		value = mc_price();

		m_Otdcrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[twoRhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate 복원
	}
 	


	delete tempirate;

		
	return m_Ogreeks[twoRhorfidx] + m_Ogreeks[twoRhodcidx];

}


/*
double AirbagStepDownB::getCorrDelta(int corrdeltaN,int *corrdeltaidx) {
	
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
*/

 

double AirbagStepDownB::mc_price()
{

 	long *rsequns;
	double **randoms;
	long nostep;

	double simsval[AsNUM];
	double tempave[AsNUM];
	double temppsval[AsNUM][10]; // 최대 10 영업일 가정

	double ***Cholcor;
	double ***tCholcor;
	
	int simkiFlag;
	int simairbagN;
	double simprice;
	double price;
	double dT, sqdT;
		
	int i,j,k;

	dT = 1.0/YFactor;
	sqdT = sqrt(dT);
	double mudt[AsNUM];
	double vldt[AsNUM];
	//m_modeltype = 0;

	double disc_div[AsNUM];
	double Srv;
	double *volskew;

	double r[AsNUM], q[AsNUM], v[AsNUM];
	 
	double *dailyIFrf, *dailyZdc;

 
	

	nostep = MAX(m_matu[m_matuN-1],1);
	 

	rsequns = (long*)malloc((size_t)((nostep)*sizeof(long)));
	randoms =(double **)malloc((size_t)((AsNUM)*sizeof(double *)));	
	if(!rsequns || !randoms) return -100;

	for ( i = 0; i < AsNUM; i++ ) {
		randoms[i] = (double *)malloc((size_t)((m_noiters)*sizeof(double)));		 
		if(!randoms[i]) {
			for(j=0; j<i; j++)
				free(randoms[j]);
			return -200;
		}
	}

	if ( rnd->randnum(1,m_noiters,AsNUM,randoms,0) != 0 ) {
		free(rsequns);
		for(i=0; i<AsNUM; i++)
			free(randoms[i]);
		free(randoms);
		return(-600);
	}

	dailyIFrf = new double[nostep];		
	
	dailyZdc = new double[nostep+1];

	for(k=0; k<AsNUM; k++) {
		for(i=0; i<nostep; i++) { //사용할 때 i=1:nostep i 대응 [i-1]로
			dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
		}
	}

	for(i=0; i<=nostep; i++) { //사용할때 i 대응 [i]로 사용
		dailyZdc[i] = interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0);
	}

	
	Cholcor = new double**[m_corrN];
	for(i=0; i<m_corrN; i++) {
		Cholcor[i] = new double*[AsNUM];
		for(j=0; j<AsNUM; j++) {
			Cholcor[i][j] = new double[AsNUM];
			for(k=0; k<AsNUM; k++) 
				Cholcor[i][j][k] = 0.0;				 			
		}
	}
	
	tCholcor = new double**[m_corrN];
	for(i=0; i<m_corrN; i++) {
		tCholcor[i] = new double*[AsNUM];
		for(j=0; j<AsNUM; j++) {
			tCholcor[i][j] = new double[AsNUM];
			for(k=0; k<AsNUM; k++) {				 
				if(j==k) 
					tCholcor[i][j][k] = 1.0;	
				else				
					tCholcor[i][j][k] = m_corr[i];		
			}
		}
	}

	volskew  = new double[m_volSN];
	 
	 

	int coi, coj, cok;
	double tempsum;
	for(i=0; i<m_corrN; i++) {
		for(cok=0; cok<AsNUM; cok++) {
			tempsum = 0;
			for(coi=0; coi<cok; coi++)
				tempsum = tempsum + Cholcor[i][cok][coi]*Cholcor[i][cok][coi];

			Cholcor[i][cok][cok] = sqrt(tCholcor[i][cok][cok]-tempsum);
						
			for(coi=cok+1; coi<AsNUM; coi++) {
				tempsum = 0;
				for(coj=0; coj<cok; coj++)
					tempsum = tempsum + Cholcor[i][coi][coj]*Cholcor[i][cok][coj];

				Cholcor[i][coi][cok] = (tCholcor[i][coi][cok] - tempsum)/Cholcor[i][cok][cok];
			}
		}
	
	}//Chol 분해!!

	
	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/kzhere1.txt",ios::out);	 
	 
			 
			logsgnuplot << "Cholcor[0][0][0] = " << Cholcor[0][0][0] << endl;
			logsgnuplot << "Cholcor[0][0][1] = " << Cholcor[0][0][1] << endl;
			logsgnuplot << "Cholcor[0][1][0] = " << Cholcor[0][1][0] << endl;
			logsgnuplot << "Cholcor[0][1][1] = " << Cholcor[0][1][1] << endl;
				 
			logsgnuplot << "Cholcor[1][0][0] = " << Cholcor[1][0][0] << endl;
			logsgnuplot << "Cholcor[1][0][1] = " << Cholcor[1][0][1] << endl;
			logsgnuplot << "Cholcor[1][1][0] = " << Cholcor[1][1][0] << endl;
			logsgnuplot << "Cholcor[1][1][1] = " << Cholcor[1][1][1] << endl;
 

		 
			logsgnuplot.close();

 
		}




	for(i=0; i<AsNUM; i++)
		q[i] = m_divrate[i];

	int si,l,sk,disci;
	int CntAs;
	double crand;
	int curtCi, crandi;
	double Wps;
	 

	 


 /*
	int evali;
	int evalq;		 
	int curtN;
	 
	for(i=0; i<m_matuN; i++) {
		if(m_matu[i]>=0) {
			curtN = i+1;
			break;
		}
	}
	 

	evalq = m_evalN[curtN-1] - m_matu[curtN-1];
	
	if(evalq <=1 ) {	
		for(i=0; i<AsNUM; i++) {
			tempave[i] = m_sval[i];
		}		 
	}
	else {
		for(i=0; i<AsNUM; i++) {
			tempave[i] = m_sval[i]*(m_evalN[curtN-1] - (evalq-1));
			for(evali=1; evali<evalq; evali++)
				tempave[i] = tempave[i] + m_psval[i][evali];
			tempave[i] = tempave[i]/m_evalN[curtN-1];
		}
	}	   
*/		
	/*
	using namespace std;
	using std::ofstream;	 
	ofstream logsgnuplot;	 
	logsgnuplot.open("c:/logland/kzhere1.txt",ios::out);	 
	 */


			 
			 
	
	if(0) {
			using namespace std;
			using std::ofstream;	 
			ofstream logsgnuplot;	 
			logsgnuplot.open("c:/logland/kzhere2.txt",ios::out);	 
	 
			 
			logsgnuplot << "m_kihval[0] = " << m_kihval[0] << endl;	 
			logsgnuplot << "m_kihval[1] = " << m_kihval[1] << endl;
	 
			logsgnuplot << "m_airbagxval[0] = " << m_airbagxval[0] << endl;
			logsgnuplot << "m_airbagxval[1] = " << m_airbagxval[1] << endl;
		 
			logsgnuplot << "m_xval[0][0] = " << m_xval[0][0] << endl;
			logsgnuplot << "m_xval[0][1] = " << m_xval[0][1] << endl;
			logsgnuplot << "m_xval[0][2] = " << m_xval[0][2] << endl;
			logsgnuplot << "m_xval[0][3] = " << m_xval[0][3] << endl;
			logsgnuplot << "m_xval[0][4] = " << m_xval[0][4] << endl;
			logsgnuplot << "m_xval[0][5] = " << m_xval[0][5] << endl;
			logsgnuplot << "m_xval[0][6] = " << m_xval[0][6] << endl;
			logsgnuplot << "m_xval[0][7] = " << m_xval[0][7] << endl;
			logsgnuplot << "m_xval[0][8] = " << m_xval[0][8] << endl;

			logsgnuplot << "m_xval[1][0] = " << m_xval[1][0] << endl;
			logsgnuplot << "m_xval[1][1] = " << m_xval[1][1] << endl;
			logsgnuplot << "m_xval[1][2] = " << m_xval[1][2] << endl;
			logsgnuplot << "m_xval[1][3] = " << m_xval[1][3] << endl;
			logsgnuplot << "m_xval[1][4] = " << m_xval[1][4] << endl;
			logsgnuplot << "m_xval[1][5] = " << m_xval[1][5] << endl;
			logsgnuplot << "m_xval[1][6] = " << m_xval[1][6] << endl;
			logsgnuplot << "m_xval[1][7] = " << m_xval[1][7] << endl;
			logsgnuplot << "m_xval[1][8] = " << m_xval[1][8] << endl;


			logsgnuplot << "m_airbagN = " << m_airbagN << endl;
			logsgnuplot << "m_matuN = " << m_matuN << endl;

			logsgnuplot << "m_Cpn[0] = " << m_Cpn[0] << endl;
			logsgnuplot << "m_Cpn[1] = " << m_Cpn[1] << endl;
			logsgnuplot << "m_Cpn[2] = " << m_Cpn[2] << endl;
			logsgnuplot << "m_Cpn[3] = " << m_Cpn[3] << endl;
			logsgnuplot << "m_Cpn[4] = " << m_Cpn[4] << endl;
			logsgnuplot << "m_Cpn[5] = " << m_Cpn[5] << endl;
			logsgnuplot << "m_Cpn[6] = " << m_Cpn[6] << endl;
			logsgnuplot << "m_Cpn[7] = " << m_Cpn[7] << endl;
			logsgnuplot << "m_Cpn[8] = " << m_Cpn[8] << endl;
			logsgnuplot << "m_Cpn[9] = " << m_Cpn[9] << endl;
			 
			 
		 
		 
			logsgnuplot.close();

 
		}
		  

  

	price = 0;
	for(si=0; si<m_noiters; si++) { //MC SIM
		simkiFlag = m_kiFlag;
		simprice = 0;
		simairbagN = m_airbagN;

		for(k=0; k<AsNUM; k++) {
			for(i=0; i<m_maxevalN; i++) {
				temppsval[k][i] = m_psval[k][i];
			}
		}

		if(rnd->randseq(1,si,nostep,rsequns) != 0) {
			for(i=0; i<AsNUM; i++) {
				free(randoms[i]);				 			 
			}
			free(randoms);
			delete dailyIFrf;
			free(rsequns);
			delete dailyZdc;
			delete volskew;
			 

			for(i=0; i<m_corrN; i++) {
				for(j=0; j<AsNUM; j++) {
					delete Cholcor[i][j];
					delete tCholcor[i][j];

				}
				delete Cholcor[i];
				delete tCholcor[i];
			}
			delete Cholcor;
			delete tCholcor;

			return -900;
		}

		l = 0;
		for(sk=0; sk<m_matuN; sk++) {
			if(m_matu[sk] < 0) 
				continue;

			for(;l<=m_matu[sk]; l++) {
				if(l==0) {
					for(k=0; k<AsNUM; k++)					
						simsval[k] = m_sval[k];
				}
				else { //if l
					for(k=0; k<AsNUM; k++) {
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

					//cor
					curtCi = curttimeCorridx(l);

					for(k=0; k<AsNUM; k++) {		
						crand = 0;
						for(crandi=0; crandi<AsNUM; crandi++) 
							crand = crand + Cholcor[curtCi][k][crandi]*randoms[crandi][rsequns[l-1]];

						simsval[k] = simsval[k] * exp(mudt[k] + vldt[k]*crand);

						if(simkiFlag == 0)
							if(simsval[k] <= m_kihval[k])
							simkiFlag = 1;
					}

					for(k=0; k<AsNUM; k++) {
						for(i=0; i<m_maxevalN-1; i++) {
							temppsval[k][m_maxevalN-1-i] = temppsval[k][m_maxevalN-1 - (i+1)];
						}
						temppsval[k][0] = simsval[k];
					}



				}//if l==0 else
			}//for l

			for(k=0; k<AsNUM; k++) {
				tempave[k] = 0;
				for(i=0; i<m_evalN[sk]; i++) 
					tempave[k] = tempave[k] + temppsval[k][i];
				tempave[k] = tempave[k]/m_evalN[sk];
			}

			CntAs = AsNUM;
			for(k=0; k<AsNUM; k++) {
				if(tempave[k] < m_xval[k][sk]) {
					CntAs--;
					break;
				}				
			}

			if(CntAs == AsNUM) {
				simprice = AMOUNT * (1+m_Cpn[sk])*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				break;
			}
			

			
			CntAs = AsNUM;
			for(k=0; k<AsNUM; k++) {
				if(tempave[k] < m_airbagxval[k]) {
					CntAs--;
					break;
				}				
			}

			if(CntAs == AsNUM) 
				simairbagN++;
			

			/*
			if(simsval[0] >= m_airbagxval[0] && simsval[1] >= m_airbagxval[1])
				simairbagN++;
			*/

			if(sk==m_matuN-1) {
				if(simkiFlag == 0)
					simprice = AMOUNT *(1+m_Cpn[sk+1])*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
				else {
					
					Wps = tempave[0]/m_bval[0];
					for(k=1; k<AsNUM; k++) 
						Wps = MIN(Wps,tempave[k]/m_bval[k]);					

					//Wps = MIN(simsval[0]/m_bval[0] -1, simsval[1]/m_bval[1]-1);
					//simprice = AMOUNT * (1+(1-(double)simairbagN/(double)m_matuN)*Wps)*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
					simprice = AMOUNT * (1+(1-1.0*simairbagN/m_matuN)*MIN((Wps-1),0))*exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
					//logsgnuplot << "simairbagN = " << simairbagN << endl;
					//logsgnuplot << "Wps = " << Wps << endl;
				}

			}
				
		} // for sk
		 
			 
 

		price = price + simprice/m_noiters;

	} // MC SIM
 
	
		 
			//logsgnuplot.close();

	
		 
		
	//메모리 해제
	for(i=0; i<AsNUM; i++) {
		free(randoms[i]);		 
	}
	free(randoms);
	delete dailyIFrf;
	free(rsequns);
	delete dailyZdc;
	delete volskew;
	 

	for(i=0; i<m_corrN; i++) {
		for(j=0; j<AsNUM; j++) {
			delete Cholcor[i][j];
			delete tCholcor[i][j];
		}
		delete Cholcor[i];
		delete tCholcor[i];
	}
	delete Cholcor;
	delete tCholcor;

	// 메모리 해제
 
		 
			 
		
	return price;
	 

} // fd_price();

 

void AirbagStepDownB::curttimeVol(double *vol, int day, int Asidx)
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
 


int AirbagStepDownB::curttimeCorridx(int day)
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
  

  