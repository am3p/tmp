#include "Plainq.h"
#include <iostream>
#include <fstream>

Plainq::Plainq(int cpFlag,  double sval,  double xval, 
	      int matu, 		   
		  int irateN,   int *irateCDay, double *iRFrate,  
		  int divrateN, int *divrateCDay, double *divrate, 
		  double divbval, int divN, int *divBDay, double *divCash, double divApy,
		  int voltype, double volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,			  
		  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN)
{		
	int i, j;
	
	//입력 데이터 내부 전역 변수로 저장
	m_cpFlag = cpFlag;
	 

	m_sval = sval;
 
	m_xval = xval;
	m_matu = matu;
	 
	 
	
	m_irateN = irateN;
	 
	m_irateCDay = new int[m_irateN];
	m_iRFrate = new double[m_irateN];
	 
	 
	for(i=0; i<m_irateN; i++) {		 
		m_irateCDay[i] = irateCDay[i];
		m_iRFrate[i] = iRFrate[i];		 
	}
	
	m_divrateN = divrateN;
	 
	m_divrateCDay = new int[m_divrateN];
	m_divrate = new double[m_divrateN];
	 
	 
	for(i=0; i<m_divrateN; i++) {		 
		m_divrateCDay[i] = divrateCDay[i];
		m_divrate[i] = divrate[i];		  
	}

	 
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
	 


	m_Sminlevel = Sminlevel;
	m_Smaxlevel = MAX(SmaxlevelC[0]*m_xval , SmaxlevelC[1]*m_sval );
	
	m_meshSM1 = SmeshM[0];
	m_meshSM2 = SmeshM[1];
	m_meshSM3 = SmeshM[2];

	m_meshTN1 = TmeshN[0];
	m_meshTN2 = TmeshN[1];
 
 
	  
}

Plainq::~Plainq()
{	
	int i;
	
	 
 
	 
	delete m_irateCDay;  
	delete m_iRFrate; 
	 
	delete m_divrateCDay;  
	delete m_divrate; 
			
	delete m_divBDay;
	delete m_divCash;
			
	delete m_volBDay;	 
	delete m_volSmness;	 

	for(i=0; i<m_volSN; i++)
		delete m_vol[i];
	delete m_vol;
 
}

double Plainq::getValue() {

	double greeks[MAXGREEKS1];
	double value;
	int i;

	 
	value = fd_price(greeks);
	

	for(i=0;i<1; i++) {
		m_Ogreeks[i] = greeks[i];
	}
 

	return value;

}
 

double Plainq::getcfValue() {

	double value;
	double r, q;
	double v;
	double d1, d2;

	r = interp1(m_matu,m_irateCDay,m_iRFrate,m_irateN,1,0);
	q = interp1(m_matu,m_divrateCDay,m_divrate,m_divrateN,1,0);
	v = m_vol[0][0];

	d1 = (log(m_sval/m_xval) + (r-q+0.5*v*v)*m_matu/365.0)/(v*sqrt(m_matu/365.0));
	d2 = d1 - v*sqrt(m_matu/365.0);
	
	value = 0.0;
	if(m_cpFlag==1 || m_cpFlag ==0) 
		value = m_sval*exp(-q*m_matu/365.0)*cnorm(d1) - m_xval*exp(-r*m_matu/365.0)*cnorm(d2);
	
	if(m_cpFlag==-1 || m_cpFlag ==0) 
		value = value + m_xval*exp(-r*m_matu/365.0)*cnorm(-d2) - m_sval*exp(-q*m_matu/365.0)*cnorm(-d1) ;
	

	 

	return value;

} 


double Plainq::fd_price(double *greeks)
{
	int i, j, k, matuLoop, intraLoop;
	int SM, SmeshType, intraTN;

	double ds, ds3, dt;

	double slevel, xlevel;
	double *S, *h;
	double *oldV, *newV;
	double *MAu, *MAd, *MAl;
	double m_bval = 1.0;

	int xidx;

	double value;

	 
	//double sleveltemp;		
	//int outN;	
	//int Baseidx;
	 
	double disc_div;
	double Srv;
	double *volskew;

	// 단순 초기화
	i = 0;
	j = 0;
	k = 0;
	matuLoop = 0;
	intraLoop = 0;	
	slevel = m_sval/m_bval;	
 

	xlevel = m_xval/m_bval;
	 
 
			
	ds = (m_Smaxlevel - m_Sminlevel)/m_meshSM1;
	xidx = (int)floor((xlevel-m_Sminlevel)/ds);



	 	if(0) {

		using namespace std;
		using std::ofstream;

		ofstream logs;
		logs.open("c:/logland/debug_pq.txt",ios::out);
		
		 

		 
			logs << "m_vol[0][0] = " <<   m_vol[0][0] << endl;
			logs << "m_sval = " <<   m_sval << endl;
			logs << "m_cpFlag = " <<   m_cpFlag << endl;
			logs << "m_matu = " <<   m_matu << endl;
			logs << "m_xval = " <<   m_xval << endl;
			logs << "m_iRFrate[0] = " <<   m_iRFrate[0] << endl;
			logs << "m_divrate[0] = " <<   m_divrate[0] << endl;
			logs << "m_Smaxlevel = " <<   m_Smaxlevel << endl;
			logs << "m_meshSM1 = " <<   m_meshSM1 << endl;
	 

		logs.close();

}
		
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
  



	
	oldV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	newV = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAu = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAd = (double *)malloc((size_t)(SM+1)*sizeof(double));
	MAl = (double *)malloc((size_t)(SM+1)*sizeof(double));

	volskew = (double*)malloc((size_t)(m_volSN)*sizeof(double));

	//만기 payoff
	if(m_cpFlag == 1) { // call option
		for(i=0; i<=SM; i++)		
			newV[i] =   MAX(S[i]-xlevel,0); // max(Sr - X,0)/B = max(Sr - X,0) /B * BA / BA = max(Sr/B *BA - X/B *BA) /BA
	}
	else { // -1일 때 put option 체크 필요..?
		if(m_cpFlag == -1) {
			for(i=0; i<=SM; i++)			
				newV[i] =   MAX(xlevel-S[i],0);
		}
		else 
			for(i=0; i<=SM; i++)			
				newV[i] =   MAX(S[i]-xlevel,0) + MAX(xlevel-S[i],0);
	} 
	 

	double r, q, v, rd;
	double *dailyIFrf ;
	double *dailyIFq;
	// intra Day에서는 입력 파라미터 고정
	// 하루마다 해당 파라미터 계산
	//q = m_divrate;

	if(m_matu>0) {
		dailyIFrf = (double *)malloc((size_t)(m_matu)*sizeof(double));		 
		dailyIFq = (double *)malloc((size_t)(m_matu)*sizeof(double));
		
		for(i=0; i<m_matu; i++) {
			dailyIFrf[i] = interp1(i+1,m_irateCDay,m_iRFrate,m_irateN,1,0)*(i+1) - interp1(i,m_irateCDay,m_iRFrate,m_irateN,1,0)*i; // 비균등 분할, 양쪽 끝 flat			 
			dailyIFq[i] = interp1(i+1,m_divrateCDay,m_divrate,m_divrateN,1,0)*(i+1) - interp1(i,m_divrateCDay,m_divrate,m_divrateN,1,0)*i;
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
		rd = r;
		q = dailyIFq[matuLoop];
 
		
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
				//v = volST(Srv,matuLoop);//interp1(Srv,m_volSmness,volskew,m_volSN,1,0); 
			}
			else { // const vol  
				v=volskew[0];
				//v=m_vol[0][0];//volskew[0] = m_vol[0][0];
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
		} //for(i) A만들기

		//경계 선형 조건
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
			memcpy(oldV,newV,(size_t)((SM+1)*sizeof(double))); //만기 payoff를 newV에 저장했으므로..
		
			//변수명 절약? 하며서 A에 LU저장 newV가 tempV 역할..
			//끝점 선형 조건으로 계산 0~SM 중에.. 1부터 SM-1 까지 풀고 0과 SM은 경계조건으로 해결

			// L*tempV = oldV forward로 풀기
			newV[1] = oldV[1]; // <- oldV에 복사했으므로 꼭 필요한 작업은 아님
			for(i=2; i<=SM-1; i++)
				newV[i] = oldV[i] - MAl[i]*newV[i-1];

			//U*newV = tempV backward로 풀기
			newV[SM-1] = newV[SM-1]/MAd[SM-1];
			for(i=SM-2; i>=1; i--)
				newV[i] = (newV[i] - MAu[i]*newV[i+1])/MAd[i];

			//경계 조건 적용
			if(SmeshType == 0) { // 균등 분할 			
				newV[0] = 2*newV[1] - newV[2];
				newV[SM] = 2*newV[SM-1] - newV[SM-2];
			}
			else { // A mesh
				newV[0] = (h[0]+h[1])/h[1] * newV[1] - h[0]/h[1] * newV[2];
				newV[SM] = (h[SM-1]+h[SM-2])/h[SM-2] * newV[SM-1] - h[SM-1]/h[SM-2] * newV[SM-2];
			}
		}//for(intraLoop)

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

		} //if(disc_div > 0.0)
		//이산배당 적용 end
		
	  
	 
				
	}//for(matuLoop)	

	if(m_matu>0) {
		free(dailyIFrf);		 
		free(dailyIFq);

	}

  
 
		value = interp1(slevel,S,newV,SM+1,SmeshType,1);
 	
		greeks[Priceidx] = value;  
		 
 
  
 
	
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

 
 
double Plainq:: volST(double S, int T)
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

double Plainq:: volinter2d(double S, int T, int SN, int TN, double *SArr, int *TArr, double **VArr)
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


void Plainq::curttimeVol(double *vol, int day)
{
	int i;
	
	 // step function 형 local vol
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
