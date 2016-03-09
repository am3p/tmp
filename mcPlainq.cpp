#include "mcPlain.h"
#include <iostream>
#include <fstream>

extern odysseyRandom *opr;

 

mcPlain::mcPlain(int cpFlag, double sval,  double xval, 		  
	      int matuN, int *matu, 	   
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

	
	for(i=0; i<AsNUM1; i++) {
		m_sval[i] = sval;	
		m_xval[i] = xval;
	}
	 
	
	m_cpFlag = cpFlag;	

	for(i=0; i<m_matuN; i++) {	 	
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

mcPlain::~mcPlain()
{	
	
	int i,k;	 
	//delete m_evalN;
	delete m_matu;

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

double mcPlain::getValue() { 
	
	double value;	 
	
	value = mc_price(); 
	

	return value;

}


 
 

double mcPlain::mc_price()
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
	
	 
	  
	price = 0;
 
 
	
		for(si=0; si<m_noiters; si++) { //MC SIM
			 
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


					}//if l==0 else

				}//for l
 				 

				if(sk==m_matuN-1) {
					
					if(m_cpFlag == 1) {
						simprice =  AMOUNT * MAX(0,simsval[0]-m_xval[0]) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
					 
					}
					else {
						simprice = AMOUNT * MAX(0,m_xval[0] - simsval[0]) * exp(-dailyZdc[m_matu[sk]]*m_matu[sk]/YFactor);
					}
				}
				
			} // for sk

			price = price + simprice/m_noiters;

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
	 

	// 메모리 해제
 
		 
			
	return price;
	 		



} //mc_price();



void mcPlain::curttimeVol(double *vol, int day, int Asidx)
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
double mcPlain::curttimeFXCorr2(int day)
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




int mcPlain::findinteridx(double x, double *Dx, int DxN)
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

 
 