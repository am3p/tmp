#include "TwoAssetDigital.h"
#include <iostream>
#include <fstream>

TwoAssetDigital::	TwoAssetDigital(int cpFlag, double *sval, double *bval, double *xval, double Cpn,
	      int matu, double ctime,
		  int evalN, double *psval, 
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,		  
		  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,		  		  	  
		  int corrtype, int corrN, int *corrBDay, double *corr,
		  int batchFlag, char *kCode,
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
		  double shiftRate, double spreadRate)
		  // shiftRate는 행사가 sift overhedge로 행사가를 bval 기준으로 몇 % 움직일건지 방향까지 반영한 값 예: 0.02 행사가를 오른쪽으로 2% 이동, -0.02 행사가를 왼쪽으로 2%이동
		  // spreadRate는 행사가에서 spread overhedge로 행사가에서 방향 없이 행사가에서 벗어난 정도를 기준가 기준으로 표현 call 일땐 왼쪽으로 적용
{		
	int i, j;
	
	m_shiftRate = shiftRate;
	m_spreadRate = abs(spreadRate);

	//입력 데이터 내부 전역 변수로 저장
	m_cpFlag = 1;//cpFlag; Digital Call 만 구현 해 놓은.

	m_sval1 = sval[0];
	m_sval2 = sval[1];
	m_bval1 = bval[0];
	m_bval2 = bval[1];

	if(m_bval1 == 1.0) {	// m_bval == 1 인경우는 기준가로 나누지 않는 TwoAssetDigital 옵션
		m_xval1 = xval[0]*(1+m_shiftRate);
		m_xval2 = xval[1]*(1+m_shiftRate);
		if(m_cpFlag == 1) { // call 옵션 왼쪽으로 적용		
			m_spreadS1 = m_xval1*(1-m_spreadRate);
			m_spreadS2 = m_xval2*(1-m_spreadRate);
		}
		else {
			m_spreadS1 = m_xval1*(1+m_spreadRate);
			m_spreadS2 = m_xval2*(1+m_spreadRate);
		}
	}
	else {		
		m_xval1 = xval[0] + m_bval1*m_shiftRate;
		m_xval2 = xval[1] + m_bval2*m_shiftRate;
		if(m_cpFlag == 1) {
			m_spreadS1 = m_xval1 - m_bval1*m_spreadRate;
			m_spreadS2 = m_xval2 - m_bval2*m_spreadRate;
		}
		else {
			m_spreadS1 = m_xval1 + m_bval1*m_spreadRate;
			m_spreadS2 = m_xval2 + m_bval2*m_spreadRate;
		}
	}

	m_Cpn = Cpn;
	m_matu = matu;
	m_ctime = ctime;
 
	m_evalN = evalN;
	m_psval1 = new double[m_evalN];
	m_psval2 = new double[m_evalN];
	for(i=0; i<m_evalN; i++) {
		m_psval1[i] = psval[i];
		m_psval2[i] = psval[m_evalN+i];
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
	m_Smaxlevel1 = MAX(SmaxlevelC[0]*m_xval1/m_bval1, SmaxlevelC[1]*m_sval1/m_bval1);
	m_Smaxlevel2 = MAX(SmaxlevelC[0]*m_xval2/m_bval2, SmaxlevelC[1]*m_sval2/m_bval2);
	
	m_meshSM1 = SmeshM[0];
	m_meshSM2 = SmeshM[1];
	m_meshSM3 = SmeshM[2];

	m_meshTN1 = TmeshN[0];
	m_meshTN2 = TmeshN[1];

	 

	int debugFlag;
	debugFlag = 0;
	if(debugFlag) {

		using namespace std;
		using std::ofstream;

		ofstream logs;
		logs.open("c:/logland/debug_twoassetdigitalinput.txt",ios::out);
			
		logs << "m_sval1 = " << m_sval1 << endl;
		logs << "m_bval2 = " << m_bval2 << endl;
		logs << "m_xval1 = " << m_xval1 << endl;
		logs << "m_matu = " << m_matu << endl;  
		logs << "m_ctime = " << m_ctime << endl;   
		logs << "m_irateN = " << m_irateN << endl;
	
		 
	 
		for(i=0; i<m_irateN; i++) {
			logs << "m_irateBDay[" << i <<"]" <<"= "<< m_irateBDay[i] << endl;
		}
		for(i=0; i<m_irateN; i++) {
			logs << "m_irateCDay[" << i <<"]" <<"= "<< m_irateCDay[i] << endl;
		}
		for(i=0; i<m_irateN; i++) {
			logs << "m_iRFrate[" << i <<"]" <<"= "<< m_iRFrate[i] << endl;
		}
		for(i=0; i<m_irateN; i++) {
			logs << "m_iDCrate[" << i <<"]" <<"= "<< m_iDCrate[i] << endl;
		}
		for(i=0; i<m_irateN; i++) {
			logs << "m_iIRSrate[" << i <<"]" <<"= "<< m_iIRSrate[i] << endl;			 
		}

		logs << "m_divrate1 = " << m_divrate1 << endl;
		logs << "m_divbval2 = " << m_divbval2 << endl;
		logs << "m_divN1 = " << m_divN1 << endl;
		 
	 
		for(i=0; i<m_divN1; i++) {
			logs << "m_divBDay[" << i <<"]" <<"= "<< m_divBDay1[i] << endl;
		}
		for(i=0; i<m_divN2; i++) {
			logs << "m_divCash[" << i <<"]" <<"= "<< m_divCash2[i] << endl;			 
		}
		logs << "m_divApy = " << m_divApy << endl;
		logs << "m_voltype = " << m_voltype << endl;
		logs << "m_volbval2 = " << m_volbval2 << endl;
		logs << "m_volTN = " << m_volTN << endl;
		logs << "m_volSN = " << m_volSN << endl;
		 
 
		for(i=0; i<m_volTN; i++)
			logs << "m_volBDay[" << i <<"]" <<"= "<< m_volBDay[i] << endl;
			 
 
		for(i=0; i<m_volSN; i++)
			logs << "m_volSmness[" << i <<"]" <<"= "<< m_volSmness[i] << endl;
			 
	 
		for(i=0; i<m_volSN; i++)
			for(j=0;j<m_volTN; j++)
				logs << "m_vol1[" << i <<"][" << j <<"] = "<< m_vol1[i][j] << endl;

		for(i=0; i<m_volSN; i++)
			for(j=0;j<m_volTN; j++)
				logs << "m_vol2[" << i <<"][" << j <<"] = "<< m_vol2[i][j] << endl;
				 
		logs << "m_batchFlag = " << m_batchFlag << endl;
		 
		logs << "m_Sminlevel1 = " << m_Sminlevel1 << endl;
		logs << "m_Smaxlevel2 = " << m_Smaxlevel2 << endl;
		
		logs << "m_meshSM1 = " << m_meshSM1 << endl;
		logs << "m_meshSM2 = " << m_meshSM2 << endl;
		logs << "m_meshSM3 = " << m_meshSM3 << endl;

		logs << "m_meshTN1 = " << m_meshTN1 << endl;
		logs << "m_meshTN2 = " << m_meshTN2 << endl; 

		logs << "m_batchFlag = " << m_batchFlag << endl;
		logs << "m_kCode = " << m_kCode << endl;

		logs.close();

	}
	
	 
	  
}

TwoAssetDigital::~TwoAssetDigital()
{	
	int i;
	
	 
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

	delete m_corrBDay;
	delete m_corr;

	for(i=0; i<m_volSN; i++) {
		delete m_vol1[i];
		delete m_vol2[i];
	}
	delete m_vol1;
	delete m_vol2;

	delete m_kCode;   


}

double TwoAssetDigital::getValue() {

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


double TwoAssetDigital::getVega(int vegaN, int *vegaidx) {
 
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


double TwoAssetDigital::getRho(int rhoN,int *rhoidx) {
	
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


double TwoAssetDigital::getCorrDelta(int corrdeltaN,int *corrdeltaidx) {
	
	double greeks[MAXGREEKS2];
	double value;

	double *tempcorr;
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		 

	m_curtgreek = 4;
 


	tempcorr = new double[m_corrN];
	//rf rho 시작
	for(i=0; i<m_corrN; i++) {
		tempcorr[i] = m_corr[i]; // RF rate 저장
	}
	 
	m_Ogreeks[twoCorrDeltaidx] = 0; //rf rho sum
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

		m_Otcorrdelta[i] = value - m_Ogreeks[0]; // // 0.1 term correlationdelta 
		m_Ogreeks[twoCorrDeltaidx] += m_Otcorrdelta[i];
	}
 

	for(i=0; i<m_corrN; i++) {
		m_corr[i] = tempcorr[i]; // RF rate 복원
	}
	
 

	delete tempcorr;

		
	return m_Ogreeks[twoCorrDeltaidx];

}


double TwoAssetDigital::cf_price()
{	
	return 7777;
}
 

void TwoAssetDigital::curttimeVol1(double *vol, int day)
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

void TwoAssetDigital::curttimeVol2(double *vol, int day)
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

double TwoAssetDigital::curttimeCorr(int day)
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





double TwoAssetDigital::fd_price(double *greeks)
{
	int i, j, k, matuLoop, intraLoop;
	int SM1, SM2, SmeshType1, SmeshType2, intraTN;

	double dt;
	double dx, dx3, dy ;

	double slevel1, slevel2, xlevel1, xlevel2, spreadSlevel1, spreadSlevel2;
	double *S1, *h1, *S2, *h2;
	double **oldV, **newV, **tU;
	double **XMAu, **XMAd, **XMAl;
	double **YMAu, **YMAd, **YMAl;
	double *vecV1, *vecV2;
	double **coeV1, **coeV2;

	int xidx1, xidx2;

	double value;

	double *SrTemp1, *SrTemp2, **SrVTemp, **SrV1Temp; // greek 계산을 위한 변수 당일값, 다음날 값 저장
	double sleveltemp1, sleveltemp2;		
	int outN;	
	int Baseidx;
	 
	double disc_div1, disc_div2;
	double Srv1, Srv2;
	double *volskew1, *volskew2;

	// 단순 초기화
	i = 0;
	j = 0;
	k = 0;
	matuLoop = 0;
	intraLoop = 0;
	
	 
	int evalq, evali;
	evalq = m_evalN - m_matu;
	
	if(evalq <=1 ) 	{
		slevel1 = m_sval1/m_bval1;	
		slevel2 = m_sval2/m_bval2;
	}
	else {
		slevel1 = m_sval1/m_bval1 * (m_evalN - (evalq-1));
		slevel2 = m_sval2/m_bval2 * (m_evalN - (evalq-1));
		for(evali=1; evali<evalq; evali++) {		
			slevel1 = slevel1 + m_psval1[evali]/m_bval1;
			slevel2 = slevel2 + m_psval2[evali]/m_bval2;
		}
		slevel1 = slevel1/m_evalN;
		slevel2 = slevel2/m_evalN;
	}	 

	xlevel1 = m_xval1/m_bval1;
	xlevel2 = m_xval2/m_bval2;
	spreadSlevel1 = m_spreadS1/m_bval1;
	spreadSlevel2 = m_spreadS2/m_bval2;
	
	 
	double meshxlevel1, meshxlevel2;		
	meshxlevel1 = m_xval1/m_bval1;			
	meshxlevel2 = m_xval2/m_bval2;
	 
			
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
			while (temps <= m_Smaxlevel1 + dx/8) {
				if(temps <= meshxlevel1 - Admesh)
					temph = dx;				 
				if(temps >= meshxlevel1 - Admesh)
					temph = AdmeshR;
				if(temps >= meshxlevel1 + Admesh)
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

		coeV1[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
		coeV2[i] = (double *)malloc((size_t)(SM2+1)*sizeof(double));
	}
	 
	vecV1 = (double *)malloc((size_t)(SM1+1)*sizeof(double));
	 
	vecV2 = (double *)malloc((size_t)(SM2+1)*sizeof(double));
 
	volskew1 = (double*)malloc((size_t)(m_volSN)*sizeof(double));
	volskew2 = (double*)malloc((size_t)(m_volSN)*sizeof(double));

	//만기 payoff	
	if(m_cpFlag == 1) { // call option
		for(i=0; i<=SM1; i++)	
			for(j=0; j<=SM2; j++) {
				if (S1[i]>=xlevel1 && S2[j]>=xlevel2)				
					newV[i][j] = m_Cpn; // TwoAssetDigital Call payoff		
				else
					newV[i][j] = 0.0;
			}
 
		if(m_spreadRate > 0) {
			double *spX;
			double *spY;
			double **spZ;
			double spreadvalue;			

			spX = new double[2];
			spY = new double[2];
			spZ = new double*[2];
			for (i=0; i<2; i++)
				spZ[i] = new double[2];

			for(i=0; i<=SM1; i++) {
				for(j=0; j<=SM2; j++) {
					if(S1[i]>=spreadSlevel1 && S1[i] <=xlevel1) {
						spX[0] = spreadSlevel1;
						spX[1] = xlevel1;
						if(S2[j]>=xlevel2){							 
							spY[0] = xlevel2;
							spY[1] = m_Smaxlevel2;

							spZ[0][0] = 0.0;
							spZ[0][1] = 0.0;
							spZ[1][0] = m_Cpn;
							spZ[1][1] = m_Cpn;
														 
							spreadvalue = interp2(S1[i],S2[j],spX,spY,spZ,1,1);
							newV[i][j] = MAX(newV[i][j],spreadvalue);
						}
						if(S2[j]>=spreadSlevel2 && S2[j]<xlevel2) {							 
							spY[0] = spreadSlevel2;
							spY[1] = xlevel2;

							spZ[0][0] = 0;
							spZ[0][1] = 0;
							spZ[1][0] = 0;
							spZ[1][1] = m_Cpn;
														 
							spreadvalue = interp2(S1[i],S2[j],spX,spY,spZ,1,1);
							newV[i][j] = MAX(newV[i][j],spreadvalue);
						}
					}
					if(S2[j]>=spreadSlevel2 && S2[j]<xlevel2 && S1[i]>xlevel1) {
						spX[0] = xlevel1;
						spX[1] = m_Smaxlevel1;
						spY[0] = spreadSlevel2;
						spY[1] = xlevel2;

						spZ[0][0]=0;
						spZ[0][1]=m_Cpn;
						spZ[1][0]=0;
						spZ[1][1]=m_Cpn;

						spreadvalue = interp2(S1[i],S2[j],spX,spY,spZ,1,1);
						newV[i][j] = MAX(newV[i][j],spreadvalue);
					}
				}
			}
			delete spX;
			delete spY;
			for(i=0; i<2; i++)
				delete spZ[i];
			delete spZ;

		}

	}
	else { // -1일 때 put option 체크 필요..?
	 int elsek;	
	 elsek=3;
	}

	 
			
	if(m_curtgreek == 1 && (m_matu == 1 || m_matu == 0)) { // 기본 그릭	Theta를 위한 +1일 저장
		for(i=0; i<outN + 4; i++) {	
			for(j=0; j<outN + 4; j++) {
				if(m_matu == 1) {
					if(m_ctime>15) { // 장 종료 후 는.. 오늘 종가가 내일 psval의 [0]에 들어 가는 경우
						sleveltemp1 = SrTemp1[i]/m_bval1; // *(m_evalN[curtN-1] - evalq) = 1 = m_matu[curtN-1] 해당
						for(evali=0; evali<evalq; evali++)
							sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
						sleveltemp1 = sleveltemp1 / m_evalN;

						sleveltemp2 = SrTemp2[j]/m_bval2;										
						for(evali=0; evali<evalq; evali++)
							sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
						sleveltemp2 = sleveltemp2 / m_evalN;
					}
					else { // 장중
						sleveltemp1 = (1+MIN(evalq,1))*SrTemp1[i]/m_bval1;					
						for(evali=1; evali<evalq; evali++)
							sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
						sleveltemp1 = sleveltemp1 / m_evalN;
										
						sleveltemp2 = (1+MIN(evalq,1))*SrTemp2[j]/m_bval2;					
						for(evali=2; evali<evalq; evali++)
							sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
						sleveltemp2 = sleveltemp2 / m_evalN;
					}
				}
				else {
					sleveltemp1 = SrTemp1[i]/m_bval1;
					for(evali=1; evali<evalq; evali++)
						sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
					sleveltemp1 = sleveltemp1 / m_evalN;
				
					sleveltemp2 = SrTemp2[j]/m_bval2;
					for(evali=1; evali<evalq; evali++)
						sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
					sleveltemp2 = sleveltemp2 / m_evalN;
				}			
				SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1,SM2);  
			}
		}		 
	}

	 

	double r, q1, q2, v1, v2, cor, rd;
	double *dailyIFrf, *dailyIFdc;	
	// intra Day에서는 입력 파라미터 고정
	// 하루마다 해당 파라미터 계산
	q1 = m_divrate1;
	q2 = m_divrate2;
	 
	if(m_matu>0) {
		dailyIFrf = (double *)malloc((size_t)(m_matu)*sizeof(double));
		dailyIFdc = (double *)malloc((size_t)(m_matu)*sizeof(double));
		
		for(i=0; i<m_matu; i++) {
			dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat
			dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
		}
	} //if(m_matu>0)


	// time.. 계산 시작
	for(matuLoop=m_matu-1; matuLoop>=0; matuLoop--) {
		//time Adaptive mesh
		if(matuLoop == m_matu-1) 
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
					XMAl[i][j] = -0.5*(v1*v1)*(i*i)*dt + (r-q1)*i*dt/2.0 + 0.5*cor*v1*v2*i*j*dt/2.0;
					XMAd[i][j] = v1*v1*i*i*dt + 0.5*rd*dt + 1;
					XMAu[i][j] = -0.5*(v1*v1)*(i*i)*dt - (r-q1)*i*dt/2.0 - 0.5*cor*v1*v2*i*j*dt/2.0;

					coeV1[i][j] = 0.5*cor*v1*v2*i*j*dt/2.0;
				}
				else {						
					XMAl[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*(h1[i]+h1[i-1])) + (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]) + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
					XMAd[i][j] = (v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*h1[i]) + 0.5*rd*dt + 1;
					XMAu[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i]*(h1[i]+h1[i-1])) - (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]) - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);								
					
					coeV1[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);	
				}
			} // 행렬 만들기					
			XMAl[0][j] = 0.0;				
			XMAd[0][j] = 1.0;				
			XMAu[0][j] = 0.0;		

	 


			//경계 선형 조건
			if(SmeshType1 == 0) {
				XMAl[SM1-1][j] = XMAl[SM1-1][j] - XMAu[SM1-1][j];
				XMAd[SM1-1][j] = XMAd[SM1-1][j] + 2*XMAu[SM1-1][j];
			}
			else {			 
				XMAl[SM1-1][j] = XMAl[SM1-1][j] - XMAu[SM1-1][j]*h1[SM1-1]/h1[SM1-2];
				XMAd[SM1-1][j] = XMAd[SM1-1][j] + XMAu[SM1-1][j]*(h1[SM1-1]+h1[SM1-2])/h1[SM1-2];			
			}

			//LU decomposition
			for(i=1; i<=SM1-1; i++) {
				XMAl[i][j] = XMAl[i][j]/XMAd[i-1][j];
				XMAd[i][j] = XMAd[i][j] - XMAl[i][j]*XMAu[i-1][j];
			}				
		} //X 행렬 끝

		// Y행령
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
					YMAl[i][j] = -0.5*(v2*v2)*(j*j)*dt + (r-q2)*j*dt/2.0 + 0.5*cor*v1*v2*i*j*dt/2.0;
					YMAd[i][j] = v2*v2*j*j*dt + 0.5*rd*dt + 1;
					YMAu[i][j] = -0.5*(v2*v2)*(j*j)*dt - (r-q2)*j*dt/2.0 - 0.5*cor*v1*v2*i*j*dt/2.0;						 

					coeV2[i][j] = 0.5*cor*v1*v2*i*j*dt/2.0;
				}
				else {						
					YMAl[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*(h2[j]+h2[j-1])) + (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]) + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);						
					YMAd[i][j] = (v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*h2[j]) + 0.5*rd*dt + 1;
					YMAu[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j]*(h2[j]+h2[j-1])) - (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]) - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);					

					coeV2[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);
				}
			} // 행렬 만들기					
			YMAl[i][0] = 0.0;				
			YMAd[i][0] = 1.0;				
			YMAu[i][0] = 0.0;				
	 
			//경계 선형 조건
			if(SmeshType2 == 0) {
				YMAl[i][SM2-1] = YMAl[i][SM2-1] - YMAu[i][SM2-1];
				YMAd[i][SM2-1] = YMAd[i][SM2-1] + 2*YMAu[i][SM2-1];
			}
			else {			 
				YMAl[i][SM2-1] = YMAl[i][SM2-1] - YMAu[i][SM2-1]*h2[SM2-1]/h2[SM2-2];
				YMAd[i][SM2-1] = YMAd[i][SM2-1] + YMAu[i][SM2-1]*(h2[SM2-1]+h2[SM2-2])/h2[SM2-2];			
			}

			//LU decomposition
			for(j=1; j<=SM2-1; j++) {
				YMAl[i][j] = YMAl[i][j]/YMAd[i][j-1];
				YMAd[i][j] = YMAd[i][j] - YMAl[i][j]*YMAu[i][j-1];
			}			
		} // Y 행렬 끝
 
	 
		
		for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {
		 
			//memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
			for(i=0; i<=SM1; i++) 
				for(j=0; j<=SM2; j++)
					oldV[i][j] = newV[i][j];

			for(i=0; i<=SM1; i++)
				tU[i][0]=0.0; // D.C 적용
		 
			for(j=1; j<SM2; j++) {  //y축 고정시켜서.. 			 

				for(i=1; i<SM1; i++) { //x축 풀기 위하여 행렬 만들기				 						
					vecV1[i] =oldV[i][j] + coeV1[i][j] *(-tU[i+1][j-1]+tU[i-1][j-1]);					 
				} // 행렬 만들기					 		

				tU[0][j] = 0;
				for(i=1; i<SM1; i++) {
					tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
				}
				tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
				for(i=SM1-2; i>=0; i--) {
					tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
				}
				
				if(SmeshType1 == 0)				
					tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
				else
					tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
				
			} // for(j y축.. 고정 다 품
			for(i=0; i<=SM1; i++) {
				if(SmeshType2 == 0) 
					tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];
				else
					tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];
			} // 일차 tU 완결
	 

			//////////////////////////////////////////////////////// newV 구하기!!
			for(j=0; j<=SM2; j++)
				newV[0][j]=0.0; // D.C 적용
		 
			for(i=1; i<SM1; i++) {  //y축 고정시켜서.. 				 

				for(j=1; j<SM2; j++) { //x축 풀기 위하여 행렬 만들기					 						
					vecV2[j] =tU[i][j] + coeV2[i][j] *(-newV[i-1][j+1]+newV[i-1][j-1]);					   
				} // 행렬 만들기						 		

				newV[i][0] = 0;
				for(j=1; j<SM2; j++) {
					newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];
				}
				newV[i][SM2-1] = newV[i][SM2-1]/YMAd[i][SM2-1];
				for(j=SM2-2; j>=0; j--) {
					newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];
				}
				
				if(SmeshType2 == 0)				
					newV[i][SM2] = 2*newV[i][SM2-1] - newV[i][SM2-2];
				else
					newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * newV[i][SM2-2];
				
			} // for(j y축.. 고정 다 품


			for(j=0; j<=SM2; j++) {
				if(SmeshType1 == 0) 
					newV[SM1][j] = 2*newV[SM1-1][j] - newV[SM1-2][j];
				else
					newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * newV[SM1-2][j];
			} // 일차 tU 완결
 


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
			for(i=0; i<=SM1; i++)
				for(j=0; j<=SM2; j++)
					oldV[i][j] = newV[i][j];

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

		} //if(disc_div > 0.0)
		//이산배당 적용 end
 							
	 
			 
		if(m_curtgreek == 1 && matuLoop == 1) { // 기본 그릭	Theta를 위한 +1일 저장
			for(i=0; i<outN + 4; i++) {
				for(j=0; j<outN + 4; j++) {
					if(evalq<=0) {
						sleveltemp1 = SrTemp1[i]/m_bval1;
						sleveltemp2 = SrTemp2[j]/m_bval2;
					}
					else { // n영업일 평균 처리
						if(m_ctime>15) { // 장종료
							sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN-evalq); //  = m_matu
							for(evali=0; evali<evalq; evali++)
								sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
							sleveltemp1 = sleveltemp1 / m_evalN;
														
							sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN-evalq); //  = m_matu
							for(evali=0; evali<evalq; evali++)
								sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
							sleveltemp2 = sleveltemp2 / m_evalN;
						}
						else { // 장중
							sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN-evalq+1);
							for(evali=1; evali<evalq; evali++)
								sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
							sleveltemp1 = sleveltemp1 / m_evalN;
																					
							sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN-evalq+1);
							for(evali=1; evali<evalq; evali++)
								sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
							sleveltemp2 = sleveltemp2 / m_evalN;
						}
					}				

					SrV1Temp[i][j] = interp2(sleveltemp1, sleveltemp2, S1, S2, newV, SM1, SM2); // (0: flat, 1: 선형)
				}// for(j
			} //for(i
		} //(m_curtgreek == 1 && matuLoop == 1)
						
	}//for(matuLoop)	

	if(m_matu>0) {
		free(dailyIFrf);
		free(dailyIFdc);

	}
		
	 
///////////////////////////////////////////--------------Debug	Mode s
 
	int debugFlag;
	debugFlag = 1;
	if(debugFlag) {

		using namespace std;
		using std::ofstream;

		ofstream logsf;
		ofstream logsgnuplot;
		//logsf.open("c:/logland/debug_fdprice.txt",ios::out);
		logsgnuplot.open("c:/logland/debug_twodgnufdprice.txt",ios::out);
	 
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
			 
		//logsf.close();
		logsgnuplot.close();
		 
	  


		logsgnuplot.open("c:/logland/debug_avetwodnufdprice.txt",ios::out);	 
		double stemp1,stemp2, valuetemp;
				
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

		for(i=0; i<=SM1; i++) {
			for(j=0; j<=SM2; j++) {
				if(evalq <=1) {
					stemp1 = S1[i]*m_bval1;
					stemp2 = S2[j]*m_bval2;
				}
				else {
					stemp1 = S1[i]*m_bval1 * (m_evalN-(evalq-1));			
					for(evali=1;evali<evalq; evali++)				
						stemp1 = stemp1 + m_psval1[evali];
					stemp1 = stemp1 / m_evalN;
										
					stemp2 = S2[i]*m_bval2 * (m_evalN-(evalq-1));			
					for(evali=1;evali<evalq; evali++)				
						stemp2 = stemp2 + m_psval2[evali];
					stemp2 = stemp2 / m_evalN;
				}

				valuetemp = interp2(stemp1/m_bval1, stemp2/m_bval2, S1,S2, newV, SM1,SM2);
												
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

	 
	if(m_curtgreek == 1 || m_curtgreek == 2) {
		for(i=0; i<outN + 4; i++) {		
			for(j=0; j<outN + 4; j++) {
				if(evalq <= 1) {
					sleveltemp1 = SrTemp1[i]/m_bval1;
					sleveltemp2 = SrTemp2[j]/m_bval2;
				}
				else {
					sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN - (evalq-1));
					for(evali=1; evali<evalq; evali++)
						sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
					sleveltemp1 = sleveltemp1/m_evalN;
										
					sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN - (evalq-1));
					for(evali=1; evali<evalq; evali++)
						sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
					sleveltemp2 = sleveltemp2/m_evalN;
				}
			
				SrVTemp[i][j] = interp2(sleveltemp1,sleveltemp2, S1,S2, newV, SM1, SM2); // (0: flat, 1: 선형)
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

		greeks[twoCrossGammaCidx] = (SrVTemp[Baseidx+1][Baseidx+1] - SrVTemp[Baseidx+1][Baseidx-1] - SrVTemp[Baseidx-1][Baseidx+1] + SrVTemp[Baseidx-1][Baseidx-1]) /(4* 0.01); // [9] = cross gamma cash

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
		
		//2012-05-23 1% gamma cash로 수정
		greeks[twoUpGammaCidx1] = (SrVTemp[Baseidx+2][Baseidx] - 2*SrVTemp[Baseidx+1][Baseidx] + SrVTemp[Baseidx][Baseidx]); // [5] = up gamma cash = gamma * (0.01s)^2 
		greeks[twoDownGammaCidx1] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx-1][Baseidx] + SrVTemp[Baseidx-2][Baseidx]); // [6] = down gamma cash		
		greeks[twoUpGammaCidx2] = (SrVTemp[Baseidx][Baseidx+2] - 2*SrVTemp[Baseidx][Baseidx+1] + SrVTemp[Baseidx][Baseidx]); // [7] = up gamma cash = gamma * (0.01s)^2 
		greeks[twoDownGammaCidx2] = (SrVTemp[Baseidx][Baseidx] - 2*SrVTemp[Baseidx][Baseidx-1] + SrVTemp[Baseidx][Baseidx-2]); // [8] = down gamma cash
		
		//2012-05-23 1% cross gamma cash로 수정
		greeks[twoCrossGammaCidx] = (SrVTemp[Baseidx+1][Baseidx+1] - SrVTemp[Baseidx+1][Baseidx-1] - SrVTemp[Baseidx-1][Baseidx+1] + SrVTemp[Baseidx-1][Baseidx-1]) /4; // [9] = 1% cross gamma cash 

		greeks[twoThetaidx] = (SrV1Temp[Baseidx][Baseidx] - SrVTemp[Baseidx][Baseidx]); // [10] = 1 day theta

		//2012-05-23 1% charm cash로 수정
		greeks[twoCharmCidx1] = (SrV1Temp[Baseidx+1][Baseidx]-SrV1Temp[Baseidx-1][Baseidx]-SrVTemp[Baseidx+1][Baseidx]+SrVTemp[Baseidx-1][Baseidx])/2; // [6] = 1 day 1% charm cash
		greeks[twoCharmCidx2] = (SrV1Temp[Baseidx][Baseidx+1]-SrV1Temp[Baseidx][Baseidx-1]-SrVTemp[Baseidx][Baseidx+1]+SrVTemp[Baseidx][Baseidx-1])/2; // [6] = 1 day 1% charm cash



		/*
		if(m_matu == 0)
		{
			for(i=1; i<BasicGRsN1; i++) 
				greeks[i] = 0;  // 만기에는 그릭 0 !!
		}
		*/
		 

		// free 전에.. file이나 memory로 출력 필요!!


 

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
		
			if(m_matu == 0 && m_ctime < 15) {		
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
				cor=m_corr[0]; 
	 
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
								XMAl[i][j] =-0.5*(v1*v1)*(i*i)*dt + (r-q1)*i*dt/2.0 + 0.5*cor*v1*v2*i*j*dt/2.0;
								XMAd[i][j] =v1*v1*i*i*dt + 0.5*rd*dt + 1;
								XMAu[i][j] =-0.5*(v1*v1)*(i*i)*dt - (r-q1)*i*dt/2.0 - 0.5*cor*v1*v2*i*j*dt/2.0;

								coeV1[i][j] =0.5*cor*v1*v2*i*j*dt/2.0;
							}
							else {
								XMAl[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*(h1[i]+h1[i-1])) + (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]) + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
								XMAd[i][j] = (v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*h1[i]) + 0.5*rd*dt + 1;
								XMAu[i][j] = -(v1*v1)*(S1[i]*S1[i])*dt/(h1[i-1]*(h1[i]+h1[i-1])) - (r-q1)*S1[i]*dt/(h1[i]+h1[i-1]) - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);

								coeV1[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h1[i]+h1[i-1])*h2[j-1]);
							}
						} // 행렬 만들기					
						XMAl[0][j] = 0.0;				
						XMAd[0][j] = 1.0;				
						XMAu[0][j] = 0.0;				

						//경계 선형 조건
						if(SmeshType1 == 0) {
							XMAl[SM1-1][j] = XMAl[SM1-1][j] - XMAu[SM1-1][j];
							XMAd[SM1-1][j] = XMAd[SM1-1][j] + 2*XMAu[SM1-1][j];
						}
						else {			 
							XMAl[SM1-1][j] = XMAl[SM1-1][j] - XMAu[SM1-1][j]*h1[SM1-1]/h1[SM1-2];
							XMAd[SM1-1][j] = XMAd[SM1-1][j] + XMAu[SM1-1][j]*(h1[SM1-1]+h1[SM1-2])/h1[SM1-2];			
						}

						//LU decomposition
						for(i=1; i<=SM1-1; i++) {
							XMAl[i][j] = XMAl[i][j]/XMAd[i-1][j];
							XMAd[i][j] = XMAd[i][j] - XMAl[i][j]*XMAu[i-1][j];
						}						
					} // X 행렬 만들기
 		 
					// Y 행렬 만들기
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
								YMAl[i][j] = -0.5*(v2*v2)*(j*j)*dt + (r-q2)*j*dt/2.0 + 0.5*cor*v1*v2*i*j*dt/2.0;
								YMAd[i][j] = v2*v2*j*j*dt + 0.5*rd*dt + 1;
								YMAu[i][j] = -0.5*(v2*v2)*(j*j)*dt - (r-q2)*j*dt/2.0 - 0.5*cor*v1*v2*i*j*dt/2.0;

								coeV2[i][j] = 0.5*cor*v1*v2*i*j*dt/2.0;
							}
							else {
								YMAl[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*(h2[j]+h2[j-1])) + (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]) + 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);						
								YMAd[i][j] = (v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*h2[j]) + 0.5*rd*dt + 1;
								YMAu[i][j] = -(v2*v2)*(S2[j]*S2[j])*dt/(h2[j-1]*(h2[j]+h2[j-1])) - (r-q2)*S2[j]*dt/(h2[j]+h2[j-1]) - 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);

								coeV2[i][j] = 0.5*cor*v1*v2*S1[i]*S2[j]*dt/((h2[j]+h2[j-1])*h1[i-1]);
							}
						} // 행렬 만들기					
						YMAl[i][0] = 0.0;				
						YMAd[i][0] = 1.0;				
						YMAu[i][0] = 0.0;				

						//경계 선형 조건
						if(SmeshType2 == 0) {
							YMAl[i][SM2-1] = YMAl[i][SM2-1] - YMAu[i][SM2-1];
							YMAd[i][SM2-1] = YMAd[i][SM2-1] + 2*YMAu[i][SM2-1];
						}
						else {			 
							YMAl[i][SM2-1] = YMAl[i][SM2-1] - YMAu[i][SM2-1]*h2[SM2-1]/h2[SM2-2];
							YMAd[i][SM2-1] = YMAd[i][SM2-1] + YMAu[i][SM2-1]*(h2[SM2-1]+h2[SM2-2])/h2[SM2-2];			
						}

						//LU decomposition
						for(j=1; j<=SM2-1; j++) {
							YMAl[i][j] = YMAl[i][j]/YMAd[i][j-1];
							YMAd[i][j] = YMAd[i][j] - YMAl[i][j]*YMAu[i][j-1];
						}

					} // Y행렬 만들기
			 


				for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 확인
					//memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
					for(i=0; i<=SM1; i++) 
						for(j=0; j<=SM2; j++)
							oldV[i][j] = newV[i][j];

					for(i=0; i<=SM1; i++)
						tU[i][0]=0.0; // D.C 적용
		 
					for(j=1; j<SM2; j++) {  //y축 고정시켜서.. 	
						for(i=1; i<SM1; i++) { //x축 풀기 위하여 행렬 만들기						
							vecV1[i] =oldV[i][j] + coeV1[i][j]*(-tU[i+1][j-1]+tU[i-1][j-1]);				 
						} // 행렬 만들기		

						tU[0][j] = 0;
						for(i=1; i<SM1; i++) {
							tU[i][j] = vecV1[i] - XMAl[i][j]*tU[i-1][j];
						}
						tU[SM1-1][j] = tU[SM1-1][j]/XMAd[SM1-1][j];
						for(i=SM1-2; i>=0; i--) {
							tU[i][j] = (tU[i][j] - XMAu[i][j]*tU[i+1][j])/XMAd[i][j];
						}
				
						if(SmeshType1 == 0)				
							tU[SM1][j] = 2*tU[SM1-1][j] - tU[SM1-2][j];
						else
							tU[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * tU[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * tU[SM1-2][j];
				
					} // for(j y축.. 고정 다 품

					for(i=0; i<=SM1; i++) {
						if(SmeshType2 == 0) 
							tU[i][SM2] = 2*tU[i][SM2-1] - tU[i][SM2-2];
						else
							tU[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * tU[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * tU[i][SM2-2];
					} // 일차 tU 완결


					//////////////////////////////////////////////////////// newV 구하기!!
					for(j=0; j<=SM2; j++)
						newV[0][j]=0.0; // D.C 적용
		 
					for(i=1; i<SM1; i++) {  //y축 고정시켜서.. 

						for(j=1; j<SM2; j++) { //x축 풀기 위하여 행렬 만들기								
							vecV2[j] =tU[i][j] + coeV2[i][j] *(-newV[i-1][j+1]+newV[i-1][j-1]);							 
						} // 행렬 만들기						  
 
						newV[i][0] = 0;
						for(j=1; j<SM2; j++) {
							newV[i][j] = vecV2[j] - YMAl[i][j]*newV[i][j-1];
						}
						newV[i][SM2-1] = newV[i][SM2-1]/YMAd[i][SM2-1];
						for(j=SM2-2; j>=0; j--) {
							newV[i][j] = (newV[i][j] - YMAu[i][j]*newV[i][j+1])/YMAd[i][j];
						}
				
						if(SmeshType2 == 0)				
							newV[i][SM2] = 2*newV[i][SM2-1] - newV[i][SM2-2];
						else
							newV[i][SM2] = (h2[SM2-1]+h2[SM2-2])/h2[SM2-2] * newV[i][SM2-1] - h2[SM2-1]/h2[SM2-2] * newV[i][SM2-2];
				
					} // for(j y축.. 고정 다 품
					for(j=0; j<=SM2; j++) {
						if(SmeshType1 == 0) 
							newV[SM1][j] = 2*newV[SM1-1][j] - newV[SM1-2][j];
						else
							newV[SM1][j] = (h1[SM1-1]+h1[SM1-2])/h1[SM1-2] * newV[SM1-1][j] - h1[SM1-1]/h1[SM1-2] * newV[SM1-2][j];
					} // new 완결
				 					
					if(m_curtgreek == 1) {  
						for(i=0; i<outN + 4; i++) {	
							for(j=0; j<outN + 4; j++) {
								if(evalq<=1) { 							
									sleveltemp1 = SrTemp1[i]/m_bval1;
									sleveltemp2 = SrTemp2[j]/m_bval2;
								}
								else {
									sleveltemp1 = SrTemp1[i]/m_bval1 * (m_evalN-(evalq-1));
									for(evali=1; evali<evalq; evali++)
										sleveltemp1 = sleveltemp1 + m_psval1[evali]/m_bval1;
									sleveltemp1 = sleveltemp1/m_evalN;

									sleveltemp2 = SrTemp2[j]/m_bval2 * (m_evalN-(evalq-1));
									for(evali=1; evali<evalq; evali++)
										sleveltemp2 = sleveltemp2 + m_psval2[evali]/m_bval2;
									sleveltemp2 = sleveltemp2/m_evalN;
								}

								SrV1Temp[i][j] = interp2(sleveltemp1,sleveltemp2, S1,S2, newV,SM1,SM2); // (0: flat, 1: 선형)
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
			
	 


/////

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
 
	} //if(m_curtgreek == 2) 



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

		free(XMAl[i]);
		free(XMAd[i]);
		free(XMAu[i]);

		free(YMAl[i]);
		free(YMAd[i]);
		free(YMAu[i]);

		free(coeV1[i]);
		free(coeV2[i]);
	}
	free(oldV);
	free(newV);
	free(tU);

	
	free(XMAu);
	free(XMAd);
	free(XMAl);
	free(coeV1);
	free(vecV1);	

	free(YMAu);
	free(YMAd);
	free(YMAl);
	free(coeV2);
	free(vecV2);

	free(volskew1);
	free(volskew2);

	return value;

}

