#include "Plain.h"
#include <iostream>
#include <fstream>

Plain::Plain(int cpFlag, double prate, double sval, double bval, double xval, 
	      int matu, double ctime,
		  int evalN,  double *psval,
		  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
		  double divrate, double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,	
		  int batchFlag, char *kCode,
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN)
{		
	int i, j;
	
	//�Է� ������ ���� ���� ������ ����
	m_cpFlag = cpFlag;
	m_prate = prate;

	m_sval = sval;
	m_bval = bval;
	m_xval = xval;
	m_matu = matu;
	m_ctime = ctime;

	 
	m_evalN = evalN;
	m_psval = new double[m_evalN];
	for(i=0; i<m_evalN; i++)
		m_psval[i] = psval[i];	 
	
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
	m_vol = new double*[m_volSN]; // m_volSN��.. 1�̻�
	for(i=0; i<m_volSN; i++)
		m_vol[i] = new double[m_volTN];
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[i][j] = vol[i*m_volTN+j];

	m_batchFlag = batchFlag;
	m_kCode = new char[kCodeN+1];
	strcpy_s(m_kCode,kCodeN,kCode);


	m_Sminlevel = Sminlevel;
	m_Smaxlevel = MAX(SmaxlevelC[0]*m_xval/m_bval, SmaxlevelC[1]*m_sval/m_bval);
	
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
		logs.open("c:/logland/debug_plaininput.txt",ios::out);
		
		logs << "m_cpFlag = " << m_cpFlag << endl;
		logs << "m_prate = " << m_prate << endl;
		logs << "m_evalN = " << m_evalN << endl;

		for(i=0; i<m_evalN; i++)		
			logs << "m_psval[" << i << "] = " << m_psval[i] << endl;

		logs << "m_sval = " << m_sval << endl;
		logs << "m_bval = " << m_bval << endl;
		logs << "m_xval = " << m_xval << endl;
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

		logs << "m_divrate = " << m_divrate << endl;
		logs << "m_divbval = " << m_divbval << endl;
		logs << "m_divN = " << m_divN << endl;
		 
	 
		for(i=0; i<m_divN; i++) {
			logs << "m_divBDay[" << i <<"]" <<"= "<< m_divBDay[i] << endl;
		}
		for(i=0; i<m_divN; i++) {
			logs << "m_divCash[" << i <<"]" <<"= "<< m_divCash[i] << endl;			 
		}
		logs << "m_divApy = " << m_divApy << endl;
		logs << "m_voltype = " << m_voltype << endl;
		logs << "m_volbval = " << m_volbval << endl;
		logs << "m_volTN = " << m_volTN << endl;
		logs << "m_volSN = " << m_volSN << endl;
		 
 
		for(i=0; i<m_volTN; i++)
			logs << "m_volBDay[" << i <<"]" <<"= "<< m_volBDay[i] << endl;
			 
 
		for(i=0; i<m_volSN; i++)
			logs << "m_volSmness[" << i <<"]" <<"= "<< m_volSmness[i] << endl;
			 
	 
		for(i=0; i<m_volSN; i++)
			for(j=0;j<m_volTN; j++)
				logs << "m_vol[" << i <<"][" << j <<"] = "<< m_vol[i][j] << endl;
	 
		 


		logs << "m_batchFlag = " << m_batchFlag << endl;
		logs << "m_kCode = " << m_kCode << endl; 

		logs << "m_Sminlevel = " << m_Sminlevel << endl;
		logs << "m_Smaxlevel = " << m_Smaxlevel << endl;
		
		logs << "m_meshSM1 = " << m_meshSM1 << endl;
		logs << "m_meshSM2 = " << m_meshSM2 << endl;
		logs << "m_meshSM3 = " << m_meshSM3 << endl;

		logs << "m_meshTN1 = " << m_meshTN1 << endl;
		logs << "m_meshTN2 = " << m_meshTN2 << endl; 

	 


		logs.close();

	}

	  
}

Plain::~Plain()
{	
	int i;
	
	 
	delete m_psval;

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

	delete m_kCode;

}

double Plain::getValue() {

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


double Plain::getVega(int *vegaN, int *vegaidx) {

	double greeks[MAXGREEKS1];
	double value;
	int i,j   ;

	double **tempv;
	
	  
	double Price[2];
	 
	
	m_curtgreek = 2;
		
	tempv = new double*[m_volSN]; // m_volSN��.. 1�̻�
	for(i=0; i<m_volSN; i++)
		tempv[i] = new double[m_volTN];
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[i][j]; // vol ����
	

	m_Ogreeks[Vegaidx] = 0;
	m_Ogreeks[VannaCidx] = 0;
	m_Ogreeks[ZommaCidx] = 0;

	int PerSN, PerTN;
	double *PerSArr;
	int *PerTArr;
	double **PerVArr;

	int vegaSN, vegaTN;
	int *vegaSidx, *vegaTidx;
	int vi, vj;
	
	PerSArr = new double[3];
	PerTArr = new int[3];
	PerVArr = new double*[3];
	for(i=0; i<3; i++)
		PerVArr[i] = new double[3];


	vegaSN = vegaN[0];
	vegaTN = vegaN[1];
	vegaSidx = new int[vegaSN];
	vegaTidx = new int[vegaTN];



	for(vi=0; vi<vegaSN; vi++)
		vegaSidx[vi] = vegaidx[vi];
	for(vj=0; vj<vegaTN; vj++)
		vegaTidx[vj] = vegaidx[vegaSN+vj];

	if(m_voltype == 0){
		vegaSN=1;
		vegaTN=1;
	}

	for(vi=0; vi<vegaSN;vi++) {
		for(vj=0; vj<vegaTN; vj++) {
			
			for(i=0; i<3; i++) 
				for(j=0; j<3; j++)
					PerVArr[i][j]=0.0;

			if(vegaSN==1) {
				PerSN = 1;
				PerSArr[0]=m_volSmness[vegaSidx[vi]];
				if(vegaTN==1) {
					PerTN=1;
					PerTArr[0]=m_volBDay[vegaTidx[vj]];
					PerVArr[0][0] = VegaPertubation;
				}
				else{
					if(vj==0) {
						PerTN=2;
						PerTArr[0]=m_volBDay[vegaTidx[vj]];
						PerTArr[1]=m_volBDay[vegaTidx[vj+1]];
						PerVArr[0][0]=VegaPertubation;
					}
					else{						 
						if(vj==vegaTN-1) {
							PerTN=2;						 
						}
						else{
							PerTN=3;							 
							PerTArr[2]=m_volBDay[vegaTidx[vj+1]];
						}
						PerTArr[0]=m_volBDay[vegaTidx[vj-1]];
						PerTArr[1]=m_volBDay[vegaTidx[vj]];
						PerVArr[0][1]=VegaPertubation;
					}
				}
			} // if(vegaSN==1)
			else{
				if(vegaTN==1) {
					PerTN=1;
					PerTArr[0] = m_volBDay[vegaTidx[vj]];
					if(vi==0) {
						PerSN=2;
						PerSArr[0] = m_volSmness[vegaSidx[vi]];
						PerSArr[1] = m_volSmness[vegaSidx[vi+1]];
						PerVArr[0][0] = VegaPertubation;
					}
					else {
						if(vi==vegaSN-1) {
							PerSN=2;
						}
						else{
							PerSN=3;
							PerSArr[2] = m_volSmness[vegaSidx[vi+1]];
						}
						PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
						PerSArr[1] = m_volSmness[vegaSidx[vi]];
						PerVArr[1][0] = VegaPertubation;
					}
				} //vegaTN==1
				else{
					if(vi==0) {
						PerSN=2;
						PerSArr[0] = m_volSmness[vegaSidx[vi]];
						PerSArr[1] = m_volSmness[vegaSidx[vi+1]];
						if(vj==0) {
							PerTN=2;
							PerTArr[0] = m_volBDay[vegaTidx[vj]];
							PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
							PerVArr[0][0] = VegaPertubation;
						} // vj==0
						else {
							if(vj==vegaTN-1) {
								PerTN=2;
							}
							else {
								PerTN=3;
								PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
							}
							PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
							PerTArr[1] = m_volBDay[vegaTidx[vj]];
							PerVArr[0][1] = VegaPertubation;
						} // else vj==0
					} // vi==0
					else {
						if(vi==vegaSN-1) {
							PerSN=2;
							PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
							PerSArr[1] = m_volSmness[vegaSidx[vi]];
							if(vj==0) {
								PerTN=2;
								PerTArr[0] = m_volBDay[vegaTidx[vj]];
								PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
								PerVArr[1][0] = VegaPertubation;
							} // vj==0
							else {
								if(vj==vegaTN-1) {
									PerTN=2;
								}
								else {
									PerTN=3;
									PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
								}
								PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
								PerTArr[1] = m_volBDay[vegaTidx[vj]];
								PerVArr[1][1] = VegaPertubation;
							} // else vj==0
						} // vi==vegaSN-1
						else {
							PerSN=3;
							PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
							PerSArr[1] = m_volSmness[vegaSidx[vi]];
							PerSArr[2] = m_volSmness[vegaSidx[vi+1]];
							if(vj==0) {
								PerTN=2;
								PerTArr[0] = m_volBDay[vegaTidx[vj]];
								PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
								PerVArr[1][0] = VegaPertubation;
							} // vj==0
							else {
								if(vj==vegaTN-1) {
									PerTN=2;
								}
								else {
									PerTN=3;
									PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
								}
								PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
								PerTArr[1] = m_volBDay[vegaTidx[vj]];
								PerVArr[1][1] = VegaPertubation;
							} // else vj==0
						} // else vi==vegaSN-1
					} //else vi==0					
				} // else vegaTN==1
			} // else vegaSN==1

			for(i=0; i<m_volSN; i++) 
				for(j=0; j<m_volTN; j++)
					m_vol[i][j] = tempv[i][j] + volinter2d(m_volSmness[i],m_volBDay[j],PerSN,PerTN,PerSArr,PerTArr,PerVArr);

			value = fd_price(greeks);
				
			Price[0] = value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+

///////////////////////////////////////////////
			for(i=0; i<3; i++) 
				for(j=0; j<3; j++)
					PerVArr[i][j]=0.0;

			if(vegaSN==1) {
				PerSN = 1;
				PerSArr[0]=m_volSmness[vegaSidx[vi]];
				if(vegaTN==1) {
					PerTN=1;
					PerTArr[0]=m_volBDay[vegaTidx[vj]];
					PerVArr[0][0] = -VegaPertubation;
				}
				else{
					if(vj==0) {
						PerTN=2;
						PerTArr[0]=m_volBDay[vegaTidx[vj]];
						PerTArr[1]=m_volBDay[vegaTidx[vj+1]];
						PerVArr[0][0]=-VegaPertubation;
					}
					else{						 
						if(vj==vegaTN-1) {
							PerTN=2;						 
						}
						else{
							PerTN=3;							 
							PerTArr[2]=m_volBDay[vegaTidx[vj+1]];
						}
						PerTArr[0]=m_volBDay[vegaTidx[vj-1]];
						PerTArr[1]=m_volBDay[vegaTidx[vj]];
						PerVArr[0][1]=-VegaPertubation;
					}
				}
			} // if(vegaSN==1)
			else{
				if(vegaTN==1) {
					PerTN=1;
					PerTArr[0] = m_volBDay[vegaTidx[vj]];
					if(vi==0) {
						PerSN=2;
						PerSArr[0] = m_volSmness[vegaSidx[vi]];
						PerSArr[1] = m_volSmness[vegaSidx[vi+1]];
						PerVArr[0][0] = -VegaPertubation;
					}
					else {
						if(vi==vegaSN-1) {
							PerSN=2;
						}
						else{
							PerSN=3;
							PerSArr[2] = m_volSmness[vegaSidx[vi+1]];
						}
						PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
						PerSArr[1] = m_volSmness[vegaSidx[vi]];
						PerVArr[1][0] = -VegaPertubation;
					}
				} //vegaTN==1
				else{
					if(vi==0) {
						PerSN=2;
						PerSArr[0] = m_volSmness[vegaSidx[vi]];
						PerSArr[1] = m_volSmness[vegaSidx[vi+1]];
						if(vj==0) {
							PerTN=2;
							PerTArr[0] = m_volBDay[vegaTidx[vj]];
							PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
							PerVArr[0][0] = -VegaPertubation;
						} // vj==0
						else {
							if(vj==vegaTN-1) {
								PerTN=2;
							}
							else {
								PerTN=3;
								PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
							}
							PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
							PerTArr[1] = m_volBDay[vegaTidx[vj]];
							PerVArr[0][1] = -VegaPertubation;
						} // else vj==0
					} // vi==0
					else {
						if(vi==vegaSN-1) {
							PerSN=2;
							PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
							PerSArr[1] = m_volSmness[vegaSidx[vi]];
							if(vj==0) {
								PerTN=2;
								PerTArr[0] = m_volBDay[vegaTidx[vj]];
								PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
								PerVArr[1][0] = -VegaPertubation;
							} // vj==0
							else {
								if(vj==vegaTN-1) {
									PerTN=2;
								}
								else {
									PerTN=3;
									PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
								}
								PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
								PerTArr[1] = m_volBDay[vegaTidx[vj]];
								PerVArr[1][1] = -VegaPertubation;
							} // else vj==0
						} // vi==vegaSN-1
						else {
							PerSN=3;
							PerSArr[0] = m_volSmness[vegaSidx[vi-1]];
							PerSArr[1] = m_volSmness[vegaSidx[vi]];
							PerSArr[2] = m_volSmness[vegaSidx[vi+1]];
							if(vj==0) {
								PerTN=2;
								PerTArr[0] = m_volBDay[vegaTidx[vj]];
								PerTArr[1] = m_volBDay[vegaTidx[vj+1]];
								PerVArr[1][0] = -VegaPertubation;
							} // vj==0
							else {
								if(vj==vegaTN-1) {
									PerTN=2;
								}
								else {
									PerTN=3;
									PerTArr[2] = m_volBDay[vegaTidx[vj+1]];
								}
								PerTArr[0] = m_volBDay[vegaTidx[vj-1]];
								PerTArr[1] = m_volBDay[vegaTidx[vj]];
								PerVArr[1][1] = -VegaPertubation;
							} // else vj==0
						} // else vi==vegaSN-1
					} //else vi==0					
				} // else vegaTN==1
			} // else vegaSN==1

			for(i=0; i<m_volSN; i++) 
				for(j=0; j<m_volTN; j++)
					m_vol[i][j] = tempv[i][j] + volinter2d(m_volSmness[i],m_volBDay[j],PerSN,PerTN,PerSArr,PerTArr,PerVArr);

			value = fd_price(greeks);
				
			Price[1] = value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+



/////////////////////

			m_Otvegas[vi*vegaTN+vj] = (Price[0] - Price[1])/2.0;	 // term vega

			m_Ogreeks[Vegaidx] += m_Otvegas[vi*vegaTN+vj];
		} // for vj
	} // for vi

		

 
 
		 
	
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[i][j] = tempv[i][j]; // vol ����
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	
	for(i=0; i<3; i++)
		delete PerVArr[i];
	delete PerVArr;

	delete PerSArr;
	delete PerTArr;
	
	delete vegaSidx;
	delete vegaTidx;


	return m_Ogreeks[Vegaidx];

}



double Plain::getTVega(int *vegaN, int *vegaidx) {

	double greeks[MAXGREEKS1];
	double value;
	int i,j   ;

	double **tempv;
	
	  
	double Price[2];
	 
	double tempsval;
	 
	tempv = new double*[m_volSN]; // m_volSN��.. 1�̻�
	for(i=0; i<m_volSN; i++)
		tempv[i] = new double[m_volTN];
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			tempv[i][j] = m_vol[i][j]; // vol ����
	

	m_Ogreeks[Vegaidx] = 0;
	m_Ogreeks[VannaCidx] = 0;
	m_Ogreeks[ZommaCidx] = 0;

	 

	int vegaSN, vegaTN;
	int *vegaSidx, *vegaTidx;
	int vi, vj;
	
	tempsval = m_sval;
	 


	vegaSN = vegaN[0];
	vegaTN = vegaN[1];
	vegaSidx = new int[vegaSN];
	vegaTidx = new int[vegaTN];



	for(vi=0; vi<vegaSN; vi++)
		vegaSidx[vi] = vegaidx[vi];
	for(vj=0; vj<vegaTN; vj++)
		vegaTidx[vj] = vegaidx[vegaSN+vj];

	if(m_voltype == 0){
		vegaSN=1;
		vegaTN=1;
	}

	for(vi=0; vi<vegaSN;vi++) {
		m_sval = tempsval*m_volSmness[vegaSidx[vi]];
		for(vj=0; vj<vegaTN; vj++) {
			


			for(i=0; i<m_volSN; i++) 
				for(j=0; j<m_volTN; j++)
					m_vol[i][j] = tempv[i][j];

			for(i=0; i<m_volSN; i++)
				for(j=vj; j<m_volTN; j++)
					m_vol[i][j] = m_vol[i][j] + 0.01;

			value = fd_price(greeks);
				
			Price[0] = value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+

///////////////////////////////////////////////


			for(i=0; i<m_volSN; i++) 
				for(j=0; j<m_volTN; j++)
					m_vol[i][j] = tempv[i][j];

			
			for(i=0; i<m_volSN; i++) 
				for(j=vj; j<m_volTN; j++)
					m_vol[i][j] = m_vol[i][j] - 0.01;

			value = fd_price(greeks);
				
			Price[1] = value; // udi=0: vol+ , udi=1: vol-  //j=0:s-, j=1:s0, j=2:s+



/////////////////////

			m_Otvegas[vi*vegaTN+vj] = (Price[0] - Price[1])/2.0;	 // term vega

			m_Ogreeks[Vegaidx] += m_Otvegas[vi*vegaTN+vj];
		} // for vj
	} // for vi

		 
 
 
		 
	m_sval = tempsval;
		
	for(i=0; i<m_volSN; i++)
		for(j=0;j<m_volTN; j++)
			m_vol[i][j] = tempv[i][j]; // vol ����
	 
	
	for(i=0; i<m_volSN; i++)
		delete tempv[i];
	delete tempv;
	
	  
	delete vegaSidx;
	delete vegaTidx;


	return m_Ogreeks[Vegaidx];

}


double Plain::getRho(int rhoN,int *rhoidx) {

	double greeks[MAXGREEKS1];
	double value;

	double *tempirate;
	int i,j;

	int pertubationX[3];
	double pertubationY[3];		 

	m_curtgreek = 3;
 


	tempirate = new double[m_irateN];
	//rf rho ����
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iRFrate[i]; // RF rate ����
	}
	 
	m_Ogreeks[Rhorfidx] = 0; //rf rho sum
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // �����̵� rho
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
					m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2��, 1:��յ�, 0:���� flat
				}
			}
			else {
				if(i == rhoN-1) {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];			

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;			
			
					for(j=0; j<m_irateN; j++) {							
						m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2��, 1:��յ�, 0:���� flat				
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
						m_iRFrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3��, 1:��յ�, 0:���� flat				
					}	

				}
			}
		}// if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otrfrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[Rhorfidx] += m_Otrfrhos[i];
	}
 

	for(i=0; i<m_irateN; i++) {
		m_iRFrate[i] = tempirate[i]; // RF rate ����
	}
	

//////////////////////////////
	//dc rho ����
	for(i=0; i<m_irateN; i++) {
		tempirate[i] = m_iDCrate[i]; // DC rate ����
	}
	 
	m_Ogreeks[Rhodcidx] = 0; // dc rho sum
	for(i=0; i<rhoN; i++) {
		if(rhoN == 1) { // �����̵� rho
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
					m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2��, 1:��յ�, 0:���� flat
				}
			}
			else {
				if(i == rhoN-1) {					
					pertubationX[0] = m_irateBDay[rhoidx[i-1]];			
					pertubationX[1] = m_irateBDay[rhoidx[i]];			

					pertubationY[0] = 0.0;			
					pertubationY[1] = RhoPertubation;			
			
					for(j=0; j<m_irateN; j++) {							
						m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,2,1,0); // 2:Data 2��, 1:��յ�, 0:���� flat				
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
						m_iDCrate[j] = tempirate[j] + interp1(m_irateBDay[j],pertubationX,pertubationY,3,1,0); // 3:Data 3��, 1:��յ�, 0:���� flat				
					}	

				}
			}
		} // if(rhoN == 0) else

		value = fd_price(greeks);

		m_Otdcrhos[i] = value - m_Ogreeks[0]; // 10bp term rf rho 
		m_Ogreeks[Rhodcidx] += m_Otdcrhos[i];
	}
 

			
	for(i=0; i<m_irateN; i++) {
		m_iDCrate[i] = tempirate[i]; // DC rate ����
	}
 	


	delete tempirate;



	return m_Ogreeks[Rhorfidx] + m_Ogreeks[Rhodcidx];

}

double Plain::cf_price()
{
	m_voltype = 3;
	return 7777;
}

double Plain::fd_price(double *greeks)
{
	int i, j, k, matuLoop, intraLoop;
	int SM, SmeshType, intraTN;

	double ds, ds3, dt;

	double slevel, xlevel;
	double *S, *h;
	double *oldV, *newV;
	double *MAu, *MAd, *MAl;

	int xidx;

	double value;

	double *SrTemp, *SrVTemp, *SrV1Temp; // greek ����� ���� ���� ���ϰ�, ������ �� ����
	double sleveltemp;		
	int outN;	
	int Baseidx;
	 
	double disc_div;
	double Srv;
	//double *volskew;

	// �ܼ� �ʱ�ȭ
	i = 0;
	j = 0;
	k = 0;
	matuLoop = 0;
	intraLoop = 0;	
	 
	int evalq, evali;
	evalq = m_evalN - m_matu;
	
	if(evalq <=1 ) 	
		slevel = m_sval/m_bval;	
	else {
		slevel = m_sval/m_bval * (m_evalN - (evalq-1));
		for(evali=1; evali<evalq; evali++)
			slevel = slevel + m_psval[evali]/m_bval;
		slevel = slevel/m_evalN;
	}	 

	 

	xlevel = m_xval/m_bval;
	 
 
			
	ds = (m_Smaxlevel - m_Sminlevel)/m_meshSM1;
	xidx = (int)floor((xlevel-m_Sminlevel)/ds);
	 
		
	// Adaptive mesh ���� ���� ���� ���� �յ� ���� ����
	if((m_meshSM2 <=0) || (m_meshSM3 <=0) || ds <= AdmeshR ) { // || (xidx1 - m_meshSM2 < 0) || (m_meshSM1 < xlevel1 + m_meshSM2)) { 
		SmeshType = 0; // �յ�
		SM = m_meshSM1;
		 
		S = (double *)malloc((size_t)(SM+1)*sizeof(double));
		//S mesh generate
		for(i=0; i<=SM; i++) 
			S[i] = m_Sminlevel + ds*i;
		h = (double *)malloc((size_t)(1)*sizeof(double)); // free(h)�� ���� ���� ����
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
				if(temps <= xlevel - Admesh + myeps)
					temph = ds; 
				if(temps > xlevel - Admesh + myeps)
					temph = AdmeshR;
				if(temps > xlevel + Admesh + myeps)
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
 
	// greek ����� ���� �۾�
	if(m_curtgreek == 1 || m_curtgreek == 2) { // �⺻ �׸�		 		 

		outN = 2*SoutRange + 1; // +- 20%���
		SrTemp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // up down gamma�� ���� �� �Ʒ� 2���� �� 
		SrVTemp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // ���� ����
		SrV1Temp = (double *)malloc((size_t)(outN + 4)*sizeof(double)); // ������ ����

		for(i=0; i<outN + 4; i++) {
			SrTemp[i] = m_sval*(1-(SoutRange+2)*0.01 + 0.01*i); //m_sval�� 78%~122%���� 1%�� ���� 0~44 45�� �߽� 22			
		}		 
	}




	
	oldV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	newV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAu = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAd = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAl = (double *)malloc((size_t)(SM+1)*sizeof(double));

	//volskew = (double*)malloc((size_t)(m_volSN)*sizeof(double));

	//���� payoff
	if(m_cpFlag == 1) { // call option
		for(i=0; i<=SM; i++)		
			newV[i] = m_prate * MAX(S[i]-xlevel,0); // max(Sr - X,0)/B = max(Sr - X,0) /B * BA / BA = max(Sr/B *BA - X/B *BA) /BA
	}
	else { // -1�� �� put option üũ �ʿ�..?
		for(i=0; i<=SM; i++)
			newV[i] = m_prate * MAX(xlevel-S[i],0);
	} 
		
	if(m_curtgreek == 1 && (m_matu == 1 || m_matu == 0)) { // �⺻ �׸�	Theta�� ���� +1�� ����
		for(i=0; i<outN + 4; i++) {						
			if(m_matu == 1) {
				if(m_ctime>15) { // �� ���� �� ��.. ���� ������ ���� psval�� [0]�� ��� ���� ���
					sleveltemp = SrTemp[i]/m_bval; // *(m_evalN[curtN-1] - evalq) = 1 = m_matu[curtN-1] �ش�
					for(evali=0; evali<evalq; evali++)
						sleveltemp = sleveltemp + m_psval[evali]/m_bval;
					sleveltemp = sleveltemp / m_evalN;
				}
				else { // ����
					sleveltemp = (1+MIN(evalq,1))*SrTemp[i]/m_bval;					
					for(evali=1; evali<evalq; evali++)
						sleveltemp = sleveltemp + m_psval[evali]/m_bval;
					sleveltemp = sleveltemp / m_evalN;
				}
			}
			else {
				sleveltemp = SrTemp[i]/m_bval;
				for(evali=1; evali<evalq; evali++)
					sleveltemp = sleveltemp + m_psval[evali]/m_bval;
				sleveltemp = sleveltemp / m_evalN;
			}
			
			SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: ����)
		}		 
	}




	double r, q, v, rd;
	double *dailyIFrf, *dailyIFdc;	
	// intra Day������ �Է� �Ķ���� ����
	// �Ϸ縶�� �ش� �Ķ���� ���
	q = m_divrate;

	if(m_matu>0) {
		dailyIFrf = (double *)malloc((size_t)(m_matu)*sizeof(double));
		dailyIFdc = (double *)malloc((size_t)(m_matu)*sizeof(double));
		
		for(i=0; i<m_matu; i++) {
			dailyIFrf[i] = interp1(i+1,m_irateBDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iRFrate,m_irateN,1,0)*i; // ��յ� ����, ���� �� flat
			dailyIFdc[i] = interp1(i+1,m_irateBDay,m_iDCrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateBDay,m_iDCrate,m_irateN,1,0)*i;
		}
	} //if(m_matu>0)


	// time.. ��� ����
	for(matuLoop=m_matu-1; matuLoop>=0; matuLoop--) {
		//time Adaptive mesh
		if(matuLoop == m_matu-1) 
			intraTN = m_meshTN1;		
		else  
			intraTN = m_meshTN2;		
		dt = (1.0/YFactor)/intraTN;

		r = dailyIFrf[matuLoop]; // �ݸ�: One Day implied forward rate �� �̿�
		rd = dailyIFdc[matuLoop];

		//if(m_voltype == 0) {
		//	volskew[0] = m_vol[0][0];
		//	//v=m_vol[0][0];//volskew[0] = m_vol[0][0];
		//}
		//else {
		//	curttimeVol(volskew, matuLoop);				
		//}

		// A*newV = oldV A ����
		for(i=1; i<=SM-1; i++) {
			if(m_voltype == 2) { //local vol	
				Srv = S[i]*m_bval / m_volbval;			
				//v = interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
				v = volST(Srv,matuLoop);//interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
			}
			else { // const vol or imp vol
				//v=volskew[0];
				v=m_vol[0][0];//volskew[0] = m_vol[0][0];
			}

			if(SmeshType == 0) { // u mesh
				MAl[i] = (r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;
				MAd[i] = 1 + rd*dt + (v*v)*(i*i)*dt;
				MAu[i] = -(r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;				 
			}
			else { // A mesh
				MAl[i] = (r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i-1]*(h[i]+h[i-1]));
				MAd[i] = 1 + rd*dt + (v*v)*(S[i]*S[i])*dt/(h[i-1]*h[i]);
				MAu[i] = -(r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i]*(h[i]+h[i-1]));
			}
		} //for(i) A�����

		//��� ���� ����
		if(SmeshType == 0) {
			MAd[1] = MAd[1] + 2*MAl[1];
			MAu[1] = MAu[1] - MAl[1];

			MAl[SM-1] = MAl[SM-1] - MAu[SM-1];
			MAd[SM-1] = MAd[SM-1] + 2*MAu[SM-1];
		}
		else {
			MAd[1] = MAd[1] + MAl[1]*(h[0]+h[1])/h[1];
			MAu[1] = MAu[1] - MAl[1]*h[0]/h[1];

			MAl[SM-1] = MAl[SM-1] - MAu[SM-1]*h[SM-1]/h[SM-2];
			MAd[SM-1] = MAd[SM-1] + MAu[SM-1]*(h[SM-1]+h[SM-2])/h[SM-2];			
		}

		//LU decomposition
		for(i=2; i<=SM-1; i++) {
			MAl[i] = MAl[i]/MAd[i-1];
			MAd[i] = MAd[i] - MAl[i]*MAu[i-1];
		}
		
		for(intraLoop=intraTN; intraLoop>=1; intraLoop--) {
			memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //���� payoff�� newV�� ���������Ƿ�..
		
			//������ ����? �ϸ缭 A�� LU���� newV�� tempV ����..
			//���� ���� �������� ��� 0~SM �߿�.. 1���� SM-1 ���� Ǯ�� 0�� SM�� ����������� �ذ�

			// L*tempV = oldV forward�� Ǯ��
			newV[1] = oldV[1]; // <- oldV�� ���������Ƿ� �� �ʿ��� �۾��� �ƴ�
			for(i=2; i<=SM-1; i++)
				newV[i] = oldV[i] - MAl[i]*newV[i-1];

			//U*newV = tempV backward�� Ǯ��
			newV[SM-1] = newV[SM-1]/MAd[SM-1];
			for(i=SM-2; i>=1; i--)
				newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

			//��� ���� ����
			if(SmeshType == 0) { // �յ� ���� 			
				newV[0] = 2*newV[1] - newV[2];
				newV[SM] = 2*newV[SM-1] - newV[SM-2];
			}
			else { // A mesh
				newV[0] = (h[0]+h[1])/h[1] * newV[1] - h[0]/h[1] * newV[2];
				newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
			}
		}//for(intraLoop)

		//�̻��� ����
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
						newV[i] = interp1(S[i],divSlevel,divV,2,1,0); //2���̰� ���� ���� �� 1,0) ������� 
					}
				}
			}//for(i)

		} //if(disc_div > 0.0)
		//�̻��� ���� end
		
	 
		if(m_curtgreek == 1 && matuLoop == 1) { // �⺻ �׸�	Theta�� ���� +1�� ����
			for(i=0; i<outN + 4; i++) {		
				if(evalq<=0) {
					sleveltemp = SrTemp[i]/m_bval;
				}
				else { // n������ ��� ó��
					if(m_ctime>15) { // ������
						sleveltemp = SrTemp[i]/m_bval * (m_evalN-evalq); //  = m_matu
						for(evali=0; evali<evalq; evali++)
							sleveltemp = sleveltemp + m_psval[evali]/m_bval;
						sleveltemp = sleveltemp / m_evalN;
					}
					else { // ����
						sleveltemp = SrTemp[i]/m_bval * (m_evalN-evalq+1);
						for(evali=1; evali<evalq; evali++)
							sleveltemp = sleveltemp + m_psval[evali]/m_bval;
						sleveltemp = sleveltemp / m_evalN;
					}
				}
				
				SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: ����)
			}		 
		}
		
	 
				
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
		logsgnuplot.open("c:/logland/debug_pgnufdprice.txt",ios::out);
	/*	 
		logsf << "SmeshType = " << SmeshType << endl;
		logsf << "m_Smaxlevel = " << m_Smaxlevel << endl;
		logsf << "m_Sminlevel = " << m_Sminlevel << endl;
		logsf << "slevel = " << slevel << endl;
		logsf << "xlevel = " << xlevel << endl;

		logsf << "ds = " << ds << endl;
		logsf << "ds3 = " << ds3 << endl;
		logsf << "dt = " << dt << endl;		
		logsf << "xidx = " << xidx << endl;
		logsf << "SM = " << SM << endl;
	 	 
		for(i=0; i<=SM; i++) {
			logsf << "S[" << i <<"]" <<"= "<< S[i] << endl;
		}		 	 
	 

		for(i=0; i<=SM; i++) {
			logsf << "newV[" << i <<"]" <<"= "<< newV[i] << endl;			 
		}		
	 */
		for(i=0; i<=SM; i++) {
			logsgnuplot << S[i]*m_bval <<"     "<< newV[i] << endl;
		}

		//logsf.close();
		logsgnuplot.close();
		 
	 
		logsgnuplot.open("c:/logland/debug_avepgnufdprice.txt",ios::out);	 
		double stemp, valuetemp;
		for(i=0; i<=SM; i++) {
			if(evalq <=1)
				stemp = S[i]*m_bval;
			else {
				stemp = S[i]*m_bval * (m_evalN-(evalq-1));			
				for(evali=1;evali<evalq; evali++)				
					stemp = stemp + m_psval[evali];
				stemp = stemp / m_evalN;
			}

			valuetemp = interp1(stemp/m_bval,S,newV,SM+1,SmeshType,1);
			logsgnuplot << S[i]*m_bval <<"     "<< valuetemp << endl;
		}		 
		logsgnuplot.close();
	}

///////////////////////////////////////////--------------Debug	Mode e	
		
	if(0) {
		using namespace std;
		using std::ofstream;
		ofstream logs;
		logs.open("c:/logland/debug_kzhere0.txt",ios::out);			

		logs << "evalq = " << evalq << endl;
		logs << "m_evalN = " << m_evalN << endl;
		logs << "m_psval[1] = " << m_psval[1] << endl;

		logs.close();		
	}

	 
	if(m_curtgreek == 1 || m_curtgreek == 2) {
		for(i=0; i<outN + 4; i++) {			
			if(evalq <= 1)
				sleveltemp = SrTemp[i]/m_bval;
			else {
				sleveltemp = SrTemp[i]/m_bval * (m_evalN - (evalq-1));
				for(evali=1; evali<evalq; evali++)
					sleveltemp = sleveltemp + m_psval[evali]/m_bval;
				sleveltemp = sleveltemp/m_evalN;
			}
			
			SrVTemp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: ����)
		}
		
		Baseidx = SoutRange + 2; // m_sval�� �ش��ϴ� idx �� �Ʒ� 2���� �� ����ؼ�.. +2 �ʿ�
	}

		
	if(0) {
		using namespace std;
		using std::ofstream;
		ofstream logs;
		logs.open("c:/logland/debug_kzhere0.txt",ios::out);			

		logs << "slevel = " << slevel << endl;

		logs.close();		
	}



	if(m_curtgreek == 1) { // �⺻ �׸�  
		value = interp1(slevel,S,newV,SM+1,SmeshType,1);
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
		greeks[UpDeltaCidx] = (SrVTemp[Baseidx+1] - SrVTemp[Baseidx]); // [1] = up delta cash = delta * 0.01s  //2012-05-23 ���� 1% delta cash�� ����
		//greeks[DownDeltaCidx] = (SrVTemp[Baseidx] - SrVTemp[Baseidx-1]) / 0.01; // [2] = down delta cash
		greeks[DownDeltaCidx] = (SrVTemp[Baseidx] - SrVTemp[Baseidx-1]); // [2] = down delta cash   // 2012-05-23 ����
		//greeks[UpGammaCidx] = (SrVTemp[Baseidx+2] - 2*SrVTemp[Baseidx+1] + SrVTemp[Baseidx]) / 0.01; // [3] = up gamma cash = gamma * s^2 * 0.01
		greeks[UpGammaCidx] = (SrVTemp[Baseidx+2] - 2*SrVTemp[Baseidx+1] + SrVTemp[Baseidx]); // [3] = up gamma cash = gamma * (0.01S)^2 // 2012-05-23 ����
		//greeks[DownGammaCidx] = (SrVTemp[Baseidx] - 2*SrVTemp[Baseidx-1] + SrVTemp[Baseidx-2]) / 0.01; // [4] = down gamma cash
		greeks[DownGammaCidx] = (SrVTemp[Baseidx] - 2*SrVTemp[Baseidx-1] + SrVTemp[Baseidx-2]); // [4] = down gamma cash // 2012-05-23 ����
		greeks[Thetaidx] = (SrV1Temp[Baseidx] - SrVTemp[Baseidx]); // [5] = 1 day theta
		//greeks[CharmCidx] = (SrV1Temp[Baseidx+1]-SrV1Temp[Baseidx-1]-SrVTemp[Baseidx+1]+SrVTemp[Baseidx-1])/(2*0.01); // [6] = 1 day charm cash
		greeks[CharmCidx] = (SrV1Temp[Baseidx+1]-SrV1Temp[Baseidx-1]-SrVTemp[Baseidx+1]+SrVTemp[Baseidx-1])/2; // [6] = 1 day charm cash //2012-05-23 1 day 1% charm cash ����
		
		/*
		if(m_matu == 0)
		{
			for(i=1; i<BasicGRsN1; i++) 
				greeks[i] = 0;  // ���⿡�� �׸� 0 !!
		}
		*/
		 

		// free ����.. file�̳� memory�� ��� �ʿ�!!


///////////////////////////////////////////--------------Debug	Mode output			
		debugFlag = 0;	
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf;
		 
			logsf.open("c:/logland/debug_fdout.txt",ios::out);
		 
		 
	  
	 
			for(i=0; i<BasicGRsN1; i++) {
				logsf << "greeks[" << i << "] = "<< greeks[i] << endl;
			}
			logsf << "value = "<< value << endl;
			logsf << "slevel = "<< slevel << endl;
			logsf << "xlevel = "<< xlevel << endl;
			logsf << "m_sval = "<< m_sval << endl;
			logsf << "m_bval = "<< m_bval << endl;

			logsf.close();
	 
			logsf.open("c:/logland/debug_afdout.txt",ios::out);
		 
		 
	  
	 
			for(i=0; i<outN + 4; i++) {
				logsf << SrTemp[i] <<"     "<< SrVTemp[i] << endl;
			}

			logsf.close();
		 
			logsf.open("c:/logland/debug_afdout1.txt",ios::out);
		  
		
			for(i=0; i<outN + 4; i++) {
				logsf << SrTemp[i] <<"     "<< SrV1Temp[i] << endl;
		
			}

			logsf.close();
		} //if(debugFlag)
///////////////////////////////////////////--------------Debug	Mode e	

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
		
			if(m_matu == 0 && m_ctime < 15) {		
				char *tempstring;

				tempstring = new char[kCodeN+1];
				intraTN = 24;//m_meshTN1;						 
				dt = (1.0/YFactor)/intraTN;

				r = m_iRFrate[0]; 
				rd = m_iDCrate[0];

				//if(m_voltype == 0) {
				//	volskew[0] = m_vol[0][0];
				//}
				//else {
				//	curttimeVol(volskew, 0);				
				//}

				// A*newV = oldV A ����
				for(i=1; i<=SM-1; i++) {
					if(m_voltype == 2) { //local vol	
						Srv = S[i]*m_bval / m_volbval;			
						//v =  interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
						v = volST(Srv,matuLoop);//interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
					}
					else { // const vol or imp vol
						//v = volskew[0];
						v = m_vol[0][0];//volskew[0];
					}

					if(SmeshType == 0) { // u mesh
						MAl[i] = (r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;
						MAd[i] = 1 + rd*dt + (v*v)*(i*i)*dt;
						MAu[i] = -(r-q)*i*dt/2 - (v*v)*(i*i)*dt/2;				 
					}
					else { // A mesh
						MAl[i] = (r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i-1]*(h[i]+h[i-1]));
						MAd[i] = 1 + rd*dt + (v*v)*(S[i]*S[i])*dt/(h[i-1]*h[i]);
						MAu[i] = -(r-q)*S[i]*dt/(h[i]+h[i-1]) - (v*v)*(S[i]*S[i])*dt/(h[i]*(h[i]+h[i-1]));
					}
				} //for(i) A�����

				//��� ���� ����
				if(SmeshType == 0) {
					MAd[1] = MAd[1] + 2*MAl[1];
					MAu[1] = MAu[1] - MAl[1];

					MAl[SM-1] = MAl[SM-1] - MAu[SM-1];
					MAd[SM-1] = MAd[SM-1] + 2*MAu[SM-1];
				}
				else {
					MAd[1] = MAd[1] + MAl[1]*(h[0]+h[1])/h[1];
					MAu[1] = MAu[1] - MAl[1]*h[0]/h[1];

					MAl[SM-1] = MAl[SM-1] - MAu[SM-1]*h[SM-1]/h[SM-2];
					MAd[SM-1] = MAd[SM-1] + MAu[SM-1]*(h[SM-1]+h[SM-2])/h[SM-2];			
				}

				//LU decomposition
				for(i=2; i<=SM-1; i++) {
					MAl[i] = MAl[i]/MAd[i-1];
					MAd[i] = MAd[i] - MAl[i]*MAu[i-1];
				}
		
				for(intraLoop=intraTN; intraLoop>=19; intraLoop--) {  // here 19 Ȯ��
					memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //���� payoff�� newV�� ���������Ƿ�..
		
					//������ ����? �ϸ缭 A�� LU���� newV�� tempV ����..
					//���� ���� �������� ��� 0~SM �߿�.. 1���� SM-1 ���� Ǯ�� 0�� SM�� ����������� �ذ�

					// L*tempV = oldV forward�� Ǯ��
					newV[1] = oldV[1]; // <- oldV�� ���������Ƿ� �� �ʿ��� �۾��� �ƴ�
					for(i=2; i<=SM-1; i++)
						newV[i] = oldV[i] - MAl[i]*newV[i-1];

					//U*newV = tempV backward�� Ǯ��
					newV[SM-1] = newV[SM-1]/MAd[SM-1];
					for(i=SM-2; i>=1; i--)
						newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

					//��� ���� ����
					if(SmeshType == 0) { // �յ� ���� 			
						newV[0] = 2*newV[1] - newV[2];
						newV[SM] = 2*newV[SM-1] - newV[SM-2];
					}
					else { // A mesh
						newV[0] = (h[0]+h[1])/h[1] * newV[1] - h[0]/h[1] * newV[2];
						newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
					}
						
			 
									
					if(m_curtgreek == 1) { // ���⿡
						for(i=0; i<outN + 4; i++) {	
							if(evalq<=1) 							
								sleveltemp = SrTemp[i]/m_bval;
							else {
								sleveltemp = SrTemp[i]/m_bval * (m_evalN-(evalq-1));
								for(evali=1; evali<evalq; evali++)
									sleveltemp = sleveltemp + m_psval[evali]/m_bval;
								sleveltemp = sleveltemp/m_evalN;
							}

							SrV1Temp[i] = interp1(sleveltemp,S,newV,SM+1,SmeshType,1); // (0: flat, 1: ����)
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
		for(i=0; i<3; i++) {
			greeks[i] = SrVTemp[Baseidx + (i-1)]; // greeks[0]: s-, greeks[1]: s0, greeks[2]: s+
		}
		
		
		// free ����.. file�̳� memory�� ��� �ʿ�!!
		// ��, file �����  v+ ���� v-����, i������ �ʿ�


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

		value = interp1(slevel,S,newV,SM+1,SmeshType,1); // value - m_Ogreeks[0] �� �̿��ϱ� ���Ͽ�
	}
	 

 
	
	free(S);
	free(h);
	free(oldV);
	free(newV);
	free(MAu);
	free(MAd);
	free(MAl);

	//free(volskew);

	return value;

}

 
 
double Plain:: volST(double S, int T)
{
	 

	int i,j;
	int sidx, tidx;

	double tempv;
	sidx = -1;
	for(i=0; i<m_volSN; i++)
		if(S>m_volSmness[i])
			sidx=i;

	tidx = -1;
	for(j=0; j<m_volTN; j++)
		if(T>m_volBDay[j])
			tidx = j;

	if(sidx==-1) {
		if(tidx==-1){
			tempv = m_vol[0][0];
		}
		else{
			if(tidx==m_volTN-1) {
				tempv = m_vol[0][m_volTN-1];
			}
			else {
				tempv = m_vol[0][tidx+1]*(T-m_volBDay[tidx])/(m_volBDay[tidx+1]-m_volBDay[tidx]);
				tempv = tempv + m_vol[0][tidx]*(T-m_volBDay[tidx+1])/(m_volBDay[tidx]-m_volBDay[tidx+1]);			
			}
		}
	}
	else{
		if(sidx==m_volSN-1){
			if(tidx==-1){
				tempv = m_vol[m_volSN-1][0];
			}
			else{
				if(tidx==m_volTN-1) {
					tempv = m_vol[m_volSN-1][m_volTN-1];
				}
				else {
					tempv = m_vol[m_volSN-1][tidx+1]*(T-m_volBDay[tidx])/(m_volBDay[tidx+1]-m_volBDay[tidx]);
					tempv = tempv + m_vol[m_volSN-1][tidx]*(T-m_volBDay[tidx+1])/(m_volBDay[tidx]-m_volBDay[tidx+1]);			
				}
			}
		}
		else{						
			if(tidx==-1){				
				tempv = m_vol[sidx+1][0]*(S-m_volSmness[sidx])/(m_volSmness[sidx+1]-m_volSmness[sidx]);					
				tempv = tempv + m_vol[sidx][0]*(S-m_volSmness[sidx+1])/(m_volSmness[sidx]-m_volSmness[sidx+1]);
			}
			else{
				if(tidx==m_volTN-1) {									
					tempv = m_vol[sidx+1][m_volTN-1]*(S-m_volSmness[sidx])/(m_volSmness[sidx+1]-m_volSmness[sidx]);									
					tempv = tempv + m_vol[sidx][m_volTN-1]*(S-m_volSmness[sidx+1])/(m_volSmness[sidx]-m_volSmness[sidx+1]);
				}
				else {
					tempv = m_vol[sidx+1][tidx+1]*(S-m_volSmness[sidx])*(T-m_volBDay[tidx])/((m_volSmness[sidx+1]-m_volSmness[sidx])*(m_volBDay[tidx+1]-m_volBDay[tidx]));
					tempv = tempv +m_vol[sidx][tidx+1]*(S-m_volSmness[sidx+1])*(T-m_volBDay[tidx])/((m_volSmness[sidx]-m_volSmness[sidx+1])*(m_volBDay[tidx+1]-m_volBDay[tidx]));	
					tempv = tempv +m_vol[sidx+1][tidx]*(S-m_volSmness[sidx])*(T-m_volBDay[tidx+1])/((m_volSmness[sidx+1]-m_volSmness[sidx])*(m_volBDay[tidx]-m_volBDay[tidx+1]));
					tempv = tempv +m_vol[sidx][tidx]*(S-m_volSmness[sidx+1])*(T-m_volBDay[tidx+1])/((m_volSmness[sidx]-m_volSmness[sidx+1])*(m_volBDay[tidx]-m_volBDay[tidx+1]));		
				}
			}
		}
	}

	return tempv;
	
}

double Plain:: volinter2d(double S, int T, int SN, int TN, double *SArr, int *TArr, double **VArr)
{
	 

	int i,j;
	int sidx, tidx;

	double tempv;
	sidx = -1;
	for(i=0; i<SN; i++)
		if(S>SArr[i])
			sidx=i;

	tidx = -1;
	for(j=0; j<TN; j++)
		if(T>TArr[j])
			tidx = j;

	if(sidx==-1) {
		if(tidx==-1){
			tempv = VArr[0][0];
		}
		else{
			if(tidx==TN-1) {
				tempv = VArr[0][TN-1];
			}
			else {
				tempv = VArr[0][tidx+1]*(T-TArr[tidx])/(TArr[tidx+1]-TArr[tidx]);
				tempv = tempv + VArr[0][tidx]*(T-TArr[tidx+1])/(TArr[tidx]-TArr[tidx+1]);			
			}
		}
	}
	else{
		if(sidx==SN-1){
			if(tidx==-1){
				tempv = VArr[SN-1][0];
			}
			else{
				if(tidx==TN-1) {
					tempv = VArr[SN-1][TN-1];
				}
				else {
					tempv = VArr[SN-1][tidx+1]*(T-TArr[tidx])/(TArr[tidx+1]-TArr[tidx]);
					tempv = tempv + VArr[SN-1][tidx]*(T-TArr[tidx+1])/(TArr[tidx]-TArr[tidx+1]);			
				}
			}
		}
		else{						
			if(tidx==-1){				
				tempv = VArr[sidx+1][0]*(S-SArr[sidx])/(SArr[sidx+1]-SArr[sidx]);					
				tempv = tempv + VArr[sidx][0]*(S-SArr[sidx+1])/(SArr[sidx]-SArr[sidx+1]);
			}
			else{
				if(tidx==TN-1) {									
					tempv = VArr[sidx+1][TN-1]*(S-SArr[sidx])/(SArr[sidx+1]-SArr[sidx]);									
					tempv = tempv + VArr[sidx][TN-1]*(S-SArr[sidx+1])/(SArr[sidx]-SArr[sidx+1]);
				}
				else {
					tempv = VArr[sidx+1][tidx+1]*(S-SArr[sidx])*(T-TArr[tidx])/((SArr[sidx+1]-SArr[sidx])*(TArr[tidx+1]-TArr[tidx]));
					tempv = tempv + VArr[sidx][tidx+1]*(S-SArr[sidx+1])*(T-TArr[tidx])/((SArr[sidx]-SArr[sidx+1])*(TArr[tidx+1]-TArr[tidx]));	
					tempv = tempv + VArr[sidx+1][tidx]*(S-SArr[sidx])*(T-TArr[tidx+1])/((SArr[sidx+1]-SArr[sidx])*(TArr[tidx]-TArr[tidx+1]));
					tempv = tempv + VArr[sidx][tidx]*(S-SArr[sidx+1])*(T-TArr[tidx+1])/((SArr[sidx]-SArr[sidx+1])*(TArr[tidx]-TArr[tidx+1]));		
				}
			}
		}
	}


	return tempv;
}

//
//void Plain::curttimeVol(double *vol, int day)
//{
//	int i;
//	
//	if(m_voltype == 1) { // term imp vol
//		double v1, v2;
//		int t1, t2;
//		double *termvol;
//
//		termvol = new double[m_volTN];
//		for(i=0; i<m_volTN; i++) 
//			termvol[i] = m_vol[0][i];
//		
//		t1 = day;
//		v1 = interp1(t1,m_volBDay,termvol,m_volTN,1,0);
//		t2 = day+1;
//		v2 = interp1(t2,m_volBDay,termvol,m_volTN,1,0);
//
//		vol[0] = sqrt(MAX(minVol*minVol, (v2*v2*t2 - v1*v1*t1)/(t2-t1))); // implied forward vol
//		delete termvol;
//	}
//	else { // step function �� local vol
//		int findidx;
//		
//		findidx = m_volTN-1;
//		for(i=0; i<m_volTN; i++) {
//			if(day <= m_volBDay[i]) {
//				findidx = i;
//				break;
//			}
//		}		 
//
//		for(i=0; i<m_volSN; i++) 
//			vol[i] = m_vol[i][findidx];
//	}
//	
//}
