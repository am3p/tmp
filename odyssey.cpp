#include "odyssey.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "Thread.h"

#include "Basic.h"
#include "Plain.h"  
#include "Plainq.h"  
#include "mcPlain.h"  
#include "Digital.h"
#include "Barrier.h"
#include "BinaryBarrier.h"

#include "CliquetA.h" // prate 양수 입력: 보통 cl, 음수 입력하면 digi cl 내부적 1, -1 로 분기하고 abs(prate)로 적용됨!
#include "ZeroBond.h"
#include "FRN.h" 

#include "StepDownA.h"  // 0, 1, 2, 10, 11, 25(주의 25는 Smax 5 입력->SmeshS1도 +), 14(Super SD) -> 잔존 : 0, 1, 2, 10, 11
#include "StepDownCDA.h"  // 0, 1, 2, 10, 11, 25(주의 25는 Smax 5 입력->SmeshS1도 +), 14(Super SD)  -> 잔존 : 0, 1, 2, 10, 11

#include "TwoAssetDigital.h"  

#include "StepDownB.h"  // 0, 1, 2, 10, 11, 201, 14(Super SD)   -> 잔존 : 0, 1, 2, 10, 11, 201
#include "StepDownCDB.h"  // 0, 1, 2, 10, 11, 201, 14(Super SD)  -> 잔존 : 0, 1, 2, 10, 11, 

#include "RStepUpA.h"  // 0, 1, 10, 11 ** 발행은 0 -> 잔존 : 0
#include "RStepUpCDA.h"  // 0, 1, 10, 11 ** 발행은 0 -> 잔존 : 0

#include "HStepDownB.h" // 0

#include "StepDownCPA.h" // 0, 1, 10, 11 -> 잔존 : 0, 10
#include "StepDownCPCDA.h"  // 0, 1, 10, 11 -> 잔존 : 0, 10

#include "StepDownCPB.h"  //0, 1, 10, 11, 121(121:KO 즉시지급형) -> 잔존: 0, 10
#include "StepDownCPCDB.h" //0, 1, 10, 11, 121(121:KO 즉시지급형) -> 잔존: 0, 10

#include "StepDownCTB.h"  // 0, 1, 10, 11, 3 -> 잔존: 0, 3
#include "StepDownCTCDB.h"  // 0, 1, 10, 11, 3 -> 잔존: 0, 3

#include "StepDownBCB.h"  // 0, 1, 2, 10, 11, 201  ** 현재발행된 건 Q10
#include "StepDownBCCDB.h"  // 0, 1, 2, 10, 11, 201  ** 현재발행된 건 Q10


#include "StabilityA.h" // 0 //MC
#include "AirbagStepDownB.h" // 0:원금비보장 다이나믹 참여율 //MC 
#include "LizardStepDownB.h" // 0:원금비보장 //MC  //FD가능한 상품

#include "BarrierB.h"
#include "PlainB.h"


#include "TomA.h" // 0 //MC
#include "TomCDA.h" // 0 //MC

//3S

#include "StepDownC.h"  // 0, 10
#include "StepDownCDC.h" // 0, 10

#include "StepDownCPC.h"  // 0, 10
#include "StepDownCPCDC.h" // 0, 10

#include "VBStepDownA.h"  // 0, 1, 14(Super SD)
#include "VBStepDownCDA.h" 

#include "VBStepDownB.h"  // 0, 1, 14(Super SD)
#include "VBStepDownCDB.h"  // 0, 1, 14(Super SD)
#include "mcPlain.h" // 0 //MC


#include <WinSock.h>




//For MC
#include "odysseyRandom.h" 
odysseyRandom *opr;
 
int RandUp = 0;

// ip check 
WSADATA WSAData;
char szHostName[128] = "";
struct sockaddr_in SocketAddress;
struct hostent     *pHost        = 0;
	
char aszIPAddresses[10][16]; // maximum of ten IP addresses
char temps[16];
char veriip[10][16]={"1.29.","1.11.","172.17.","10.0."};  //otc, 리스크, IT, virtual
int ipfindi, ipi, iptargeti, ipCnt;
int veriipq;


int veriipqFlag = 1; // 1이면 check 안함 , 0 이면 체크함
// ip check

 
 double ODYSSEY_API  _stdcall PlainFn(int uAssetN, int cpFlag, double prate,
	                                  double *sval, double *bval, double *xval, int xvaltype,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int *vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN)
/*
                                     (기초자산수, cpFlag(call:1, put:-1), 현재가, 기준주가(1이면 보통의 plain), 행사가(real 또는 0.9등), 행사가(1:real, 100: 기준가 기준)
									 만기개수, 영업일기준 조기상환 포함 만기일정, 계산호출시간(15시 전이냐 후냐만 중요),
									 금리 기간구조 개수, 영업일 기준의 금리 만기들, 달력일 기준의 금리 만기들, 무위험 금리 제로 기간구조, 할인 금리 제로 기간구조, IRS 제로,
									 연속배당율(보통은 0, 필요할때만), 이산배당의 기준주가, 이산배당 개수(더미 포함한 개수), 이산배당 일정(더미 포함), 이산배당금(더미 포함), 배당적용기준,
									 변동성 타입(0:상수, 1:imp, 2:local), 변동성 기준주가, 변동성 기간구조 개수, 변동성 머니니스 개수, 변동성 기간구조 만기, 변동성 주가 머니니스, 변동성,
									 그릭산출유형(0:total, 1:value, ud, dd, ug,dg, t, c, 2: 1 + vega, va, zo, 3: 1 + rf rho, rd rho)
									 term별 vega 산출 개수, 해당 index, term별 rho 산출 개수, 해당 indx,
									 value & 그릭 출력 배열,
									 term별 베가 출력, term별 바나 출력, term별 좀마 출력,
									 term별 rf 로우 출력, term별 dc 로우 출력,
									 batch file 생성 유무 flag,
									 FDM Sminlevel, Smaxlevel 계수(max(c1 * x/b, c2 * sc/b)의 c1, c2), Smesh 개수(균등, 특이점 전후 개수, 세분), Tmesh 개수(만기, 평일) )
									 
									 예:
									PlainFn(1, 1, {276}, {276}, {276}, 1,
									  1, {252}, 15.2,
									  20, {63, 126, 189, 252, 315, 378, 441, 504, 567, 630, 693, 756, 819, 882, 945, 1008, 1071, 1134, 1197, 1260},
									  {91, 183, 274, 365, 456, 548, 639, 730, 821, 913, 1004, 1095, 1186, 1278, 1369, 1460, 1551, 1643, 1734, 1825},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0}, {276}, {2}, {-10, 126}, {0, 3}, 0.75,
									  2, {276}, 7, 5, {63, 84, 105, 126, 252, 504, 756}, {0.85, 0.9, 0.95, 1.0, 1.05},
									  {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
									   1,
									   7, {0, 1, 2, 3, 4, 5, 6}, 5, {1, 3, 5, 7, 11}
									   *Ogrkkes,
									   *Otvegas, *Otvannas, *Otzommas,
									   *Otrfrhos, *Otdcrhos,
									   0, {5, 1.5}, {10000, 100, 10}, {24, 6} )
*/

{
	if(matu[0] < 0) 
		return -999999;

	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}



		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}


	double value;	
	double inxval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	if(xvaltype == 100) {
		inxval = xval[0]*bval[0];

	}
	else {
		inxval = xval[0];
	}
  
	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 


	Plain *plain = new Plain(cpFlag, prate, sval[0], bval[0], inxval, 
	      matu[0], ctime,
		  evalN[0], psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN);
			  
	/* value & 기본 그릭 산출 */    
	value = plain->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = plain->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = plain->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = plain->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN[0]*vegaN[1]; i++) {
			Otvegas[i] = plain->m_Otvegas[i];
			//Otvannas[i] = plain->m_Otvannas[i];
			//Otzommas[i] = plain->m_Otzommas[i];
		}
	}
	

	/* tVega 관련 그릭 산출 */
	if(greektype == 22) {		 
		value = plain->getTVega(vegaN, vegaidx);
		 

		/* term별 vegas 할당 */
		for(i=0; i<vegaN[0]*vegaN[1]; i++) {
			Otvegas[i] = plain->m_Otvegas[i];
			//Otvannas[i] = plain->m_Otvannas[i];
			//Otzommas[i] = plain->m_Otzommas[i];
		}
	}

	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = plain->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = plain->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = plain->m_Otrfrhos[i];
			Otdcrhos[i] = plain->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
			
	delete psval;
	delete plain;

	return value;
		
} 
  

 double ODYSSEY_API  _stdcall PlainqFn( int cpFlag,  double *sval,  double *xval, int *matu, 
									  int irateN,  int *irateCDay, double *iRFrate,  
									  int divrateN, int *divrateCDay, double *divrate, 
									  double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,									  
									  double *Ogreeks,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN)
 
{
	 

	double value;	
  
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	 

	Plainq *plainq = new Plainq(cpFlag,  sval[0], xval[0], matu[0], 		 
		  irateN, irateCDay, iRFrate,  
		  divrateN, divrateCDay, divrate,
		  divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,			  
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN);
			  
	/* value & 기본 그릭 산출 */    
	value = plainq->getValue();	 
	for(i=0; i<1; i++) {	
		Ogreeks[i] = plainq->m_Ogreeks[i];
	}

	 
	 
	 
 
	 
	delete plainq;

	return value;
		
} 
  



 double ODYSSEY_API  _stdcall DigitalFn(int uAssetN, int cpFlag, double *sval, double *bval, double *xval, int xvaltype, 
									  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
									  double shiftRate, double spreadRate)
/*
                                     (기초자산수, cpFlag(call:1, put:-1), 현재가, 기준주가(1이면 보통의 digital), 행사가(real 또는 0.9등), 행사가기준(1:real, 100: 기준가 기준), 
									 Cpn (액면 1원 기준의 행사 영역에서 받는 금액 예: 기준가의 130% 이상이면 액면의 3% 준다고 하면 0.03 입력)
									 만기개수, 영업일기준 조기상환 포함 만기일정, 계산호출시간(15시 전이냐 후냐만 중요),
									 금리 기간구조 개수, 영업일 기준의 금리 만기들, 달력일 기준의 금리 만기들, 무위험 금리 제로 기간구조, 할인 금리 제로 기간구조, IRS 제로,
									 연속배당율(보통은 0, 필요할때만), 이산배당의 기준주가, 이산배당 개수(더미 포함한 개수), 이산배당 일정(더미 포함), 이산배당금(더미 포함), 배당적용기준,
									 변동성 타입(0:상수, 1:imp, 2:local), 변동성 기준주가, 변동성 기간구조 개수, 변동성 머니니스 개수, 변동성 기간구조 만기, 변동성 주가 머니니스, 변동성,
									 그릭산출유형(0:total, 1:value, ud, dd, ug,dg, t, c, 2: 1 + vega, va, zo, 3: 1 + rf rho, rd rho)
									 term별 vega 산출 개수, 해당 index, term별 rho 산출 개수, 해당 indx,
									 value & 그릭 출력 배열,
									 term별 베가 출력, term별 바나 출력, term별 좀마 출력,
									 term별 rf 로우 출력, term별 dc 로우 출력,
									 batch file 생성 유무 flag,
									 FDM Sminlevel, Smaxlevel 계수(max(c1 * x/b, c2 * sc/b)의 c1, c2), Smesh 개수(균등, 특이점 전후 개수, 세분), Tmesh 개수(만기, 평일),
									 shiftRAte(방향 염두), spreadRate(방향없이, call은 왼쪽으로 풋은 오른쪽으로 적용))
									 
									 예:
									DigitalFn(1, 1, {276}, {276}, {276}, 1,
									  0.03,
									  1, {252}, 15.2,
									  20, {63, 126, 189, 252, 315, 378, 441, 504, 567, 630, 693, 756, 819, 882, 945, 1008, 1071, 1134, 1197, 1260},
									  {91, 183, 274, 365, 456, 548, 639, 730, 821, 913, 1004, 1095, 1186, 1278, 1369, 1460, 1551, 1643, 1734, 1825},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0}, {276}, {2}, {-10, 126}, {0, 3}, 0.75,
									  2, {276}, 7, 5, {63, 84, 105, 126, 252, 504, 756}, {0.85, 0.9, 0.95, 1.0, 1.05},
									  {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
									   1,
									   7, {0, 1, 2, 3, 4, 5, 6}, 5, {1, 3, 5, 7, 11}
									   *Ogrkkes,
									   *Otvegas, *Otvannas, *Otzommas,
									   *Otrfrhos, *Otdcrhos,
									   0, {5, 1.5}, {10000, 100, 10}, {24, 6},
									   0.02, 0.01)
*/

{
	if(matu[0] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double inxval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	if(xvaltype == 100) {
		inxval = xval[0]*bval[0];

	}
	else {
		inxval = xval[0];
	}
  
	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 


	Digital *digital = new Digital(cpFlag, sval[0], bval[0], inxval, Cpn[0], 
	      matu[0], ctime,
		  evalN[0], psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN,
		  shiftRate, spreadRate);
			  
	/* value & 기본 그릭 산출 */    
	value = digital->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = digital->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = digital->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = digital->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = digital->m_Otvegas[i];
			Otvannas[i] = digital->m_Otvannas[i];
			Otzommas[i] = digital->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = digital->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = digital->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = digital->m_Otrfrhos[i];
			Otdcrhos[i] = digital->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
			
	delete psval;
	delete digital;

	return value;
		
} 
 






 double ODYSSEY_API  _stdcall BarrierFn(int uAssetN, int duFlag, int ioFlag, int cpFlag, int knockFlag, double prate,
	                                  double *sval, double *bval, double *xval, double *hval, int xvaltype, 									  
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
									  double shiftRate, double spreadRate)
/*
                                     (기초자산수, cpFlag(call:1, put:-1), 현재가, 기준주가(1이면 보통의 barrier), 행사가(real 또는 0.9등), 행사가기준(1:real, 100: 기준가 기준), 
									 Cpn (액면 1원 기준의 행사 영역에서 받는 금액 예: 기준가의 130% 이상이면 액면의 3% 준다고 하면 0.03 입력)
									 만기개수, 영업일기준 조기상환 포함 만기일정, 계산호출시간(15시 전이냐 후냐만 중요),
									 금리 기간구조 개수, 영업일 기준의 금리 만기들, 달력일 기준의 금리 만기들, 무위험 금리 제로 기간구조, 할인 금리 제로 기간구조, IRS 제로,
									 연속배당율(보통은 0, 필요할때만), 이산배당의 기준주가, 이산배당 개수(더미 포함한 개수), 이산배당 일정(더미 포함), 이산배당금(더미 포함), 배당적용기준,
									 변동성 타입(0:상수, 1:imp, 2:local), 변동성 기준주가, 변동성 기간구조 개수, 변동성 머니니스 개수, 변동성 기간구조 만기, 변동성 주가 머니니스, 변동성,
									 그릭산출유형(0:total, 1:value, ud, dd, ug,dg, t, c, 2: 1 + vega, va, zo, 3: 1 + rf rho, rd rho)
									 term별 vega 산출 개수, 해당 index, term별 rho 산출 개수, 해당 indx,
									 value & 그릭 출력 배열,
									 term별 베가 출력, term별 바나 출력, term별 좀마 출력,
									 term별 rf 로우 출력, term별 dc 로우 출력,
									 batch file 생성 유무 flag,
									 FDM Sminlevel, Smaxlevel 계수(max(c1 * x/b, c2 * sc/b)의 c1, c2), Smesh 개수(균등, 특이점 전후 개수, 세분), Tmesh 개수(만기, 평일),
									 shiftRAte(방향 염두), spreadRate(방향없이, call은 왼쪽으로 풋은 오른쪽으로 적용))
									 
									 예:
									BarrierFn(1, 1, {276}, {276}, {276}, 1,
									  0.03,
									  1, {252}, 15.2,
									  20, {63, 126, 189, 252, 315, 378, 441, 504, 567, 630, 693, 756, 819, 882, 945, 1008, 1071, 1134, 1197, 1260},
									  {91, 183, 274, 365, 456, 548, 639, 730, 821, 913, 1004, 1095, 1186, 1278, 1369, 1460, 1551, 1643, 1734, 1825},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0}, {276}, {2}, {-10, 126}, {0, 3}, 0.75,
									  2, {276}, 7, 5, {63, 84, 105, 126, 252, 504, 756}, {0.85, 0.9, 0.95, 1.0, 1.05},
									  {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
									   1,
									   7, {0, 1, 2, 3, 4, 5, 6}, 5, {1, 3, 5, 7, 11}
									   *Ogrkkes,
									   *Otvegas, *Otvannas, *Otzommas,
									   *Otrfrhos, *Otdcrhos,
									   0, {5, 1.5}, {10000, 100, 10}, {24, 6},
									   0.02, 0.01)
*/

{
	if(matu[0] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double inxval;
	double inhval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	if(xvaltype == 100) {
		inxval = xval[0]*bval[0];
		inhval = hval[0]*bval[0];

	}
	else {
		inxval = xval[0];
		inhval = hval[0];
	}
  
	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 


	Barrier *barrier = new Barrier(duFlag, ioFlag, cpFlag, knockFlag, prate, 
		  sval[0], bval[0], inxval, inhval,
	      matu[0], ctime,	
		  evalN[0], psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN,
		  shiftRate, spreadRate);
			  
	/* value & 기본 그릭 산출 */    
	value = barrier->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = barrier->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = barrier->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = barrier->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = barrier->m_Otvegas[i];
			Otvannas[i] = barrier->m_Otvannas[i];
			Otzommas[i] = barrier->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = barrier->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = barrier->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = barrier->m_Otrfrhos[i];
			Otdcrhos[i] = barrier->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 		
	delete psval;
	delete barrier;

	return value;
		
} 
  





 double ODYSSEY_API  _stdcall BinaryBarrierFn(int uAssetN, int duFlag, int ioFlag, int cpFlag, int knockFlag, 
	                                  double *sval, double *bval, double *xval, double *hval, int xvaltype, 		
									  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
									  double shiftHRate, double shiftXRate, double spreadHRate, double spreadXRate)
{
	if(matu[0] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double inxval;
	double inhval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	if(xvaltype == 100) {
		inxval = xval[0]*bval[0];
		inhval = hval[0]*bval[0];

	}
	else {
		inxval = xval[0];
		inhval = hval[0];
	}
  
	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 


	BinaryBarrier *binarybarrier = new BinaryBarrier(duFlag, ioFlag, cpFlag, knockFlag, sval[0], bval[0], inxval, inhval, Cpn[0],
	      matu[0], ctime,
		  evalN[0], psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN,
		  shiftHRate, shiftXRate, spreadHRate, spreadXRate);
			  
	/* value & 기본 그릭 산출 */    
	value = binarybarrier->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = binarybarrier->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = binarybarrier->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = binarybarrier->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = binarybarrier->m_Otvegas[i];
			Otvannas[i] = binarybarrier->m_Otvannas[i];
			Otzommas[i] = binarybarrier->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = binarybarrier->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = binarybarrier->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = binarybarrier->m_Otrfrhos[i];
			Otdcrhos[i] = binarybarrier->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
		
	delete psval;
	delete binarybarrier;

	return value;
		
} 
 



 double ODYSSEY_API  _stdcall StepDownAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag,  char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double inkihval;
	double *psval;
	int i;

	double inkohval;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];
		inkohval = hval[1]*bval[0];

	}
	else {
		for(i=0; i<matuN; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
		inkohval = hval[1];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 
	StepDownA *stepdowna = new StepDownA(knockFlag[0],knockFlag[1], sval[0], bval[0], inxval, inkihval, inkohval,
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	stepdowna->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowna->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = stepdowna->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowna->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowna->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = stepdowna->m_Otvegas[i];
			Otvannas[i] = stepdowna->m_Otvannas[i];
			Otzommas[i] = stepdowna->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowna->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = stepdowna->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowna->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowna->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete stepdowna;

	return value;
		
} 
  


 double ODYSSEY_API  _stdcall TwoAssetDigitalFn(int uAssetN, int cpFlag, double *sval, double *bval, double *xval, int xvaltype, 
									  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
									  double shiftRate, double spreadRate)
/*
                                     (기초자산수, cpFlag(call:1, put:-1), 현재가, 기준주가(1이면 보통의 twoassetdigital), 행사가(real 또는 0.9등), 행사가기준(1:real, 100: 기준가 기준), 
									 Cpn (액면 1원 기준의 행사 영역에서 받는 금액 예: 기준가의 130% 이상이면 액면의 3% 준다고 하면 0.03 입력)
									 만기개수, 영업일기준 조기상환 포함 만기일정, 계산호출시간(15시 전이냐 후냐만 중요),
									 금리 기간구조 개수, 영업일 기준의 금리 만기들, 달력일 기준의 금리 만기들, 무위험 금리 제로 기간구조, 할인 금리 제로 기간구조, IRS 제로,
									 연속배당율(보통은 0, 필요할때만), 이산배당의 기준주가, 이산배당 개수(더미 포함한 개수), 이산배당 일정(더미 포함), 이산배당금(더미 포함), 배당적용기준,
									 변동성 타입(0:상수, 1:imp, 2:local), 변동성 기준주가, 변동성 기간구조 개수, 변동성 머니니스 개수, 변동성 기간구조 만기, 변동성 주가 머니니스, 변동성,
									 그릭산출유형(0:total, 1:value, ud, dd, ug,dg, t, c, 2: 1 + vega, va, zo, 3: 1 + rf rho, rd rho)
									 term별 vega 산출 개수, 해당 index, term별 rho 산출 개수, 해당 indx,
									 value & 그릭 출력 배열,
									 term별 베가 출력, term별 바나 출력, term별 좀마 출력,
									 term별 rf 로우 출력, term별 dc 로우 출력,
									 batch file 생성 유무 flag,
									 FDM Sminlevel, Smaxlevel 계수(max(c1 * x/b, c2 * sc/b)의 c1, c2), Smesh 개수(균등, 특이점 전후 개수, 세분), Tmesh 개수(만기, 평일),
									 shiftRAte(방향 염두), spreadRate(방향없이, call은 왼쪽으로 풋은 오른쪽으로 적용))
									 
									 예:
									TwoAssetDigitalFn(1, 1, {276}, {276}, {276}, 1,
									  0.03,
									  1, {252}, 15.2,
									  20, {63, 126, 189, 252, 315, 378, 441, 504, 567, 630, 693, 756, 819, 882, 945, 1008, 1071, 1134, 1197, 1260},
									  {91, 183, 274, 365, 456, 548, 639, 730, 821, 913, 1004, 1095, 1186, 1278, 1369, 1460, 1551, 1643, 1734, 1825},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0.05, 0.05, 0.05, 0.05, 0.055, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
									  {0}, {276}, {2}, {-10, 126}, {0, 3}, 0.75,
									  2, {276}, 7, 5, {63, 84, 105, 126, 252, 504, 756}, {0.85, 0.9, 0.95, 1.0, 1.05},
									  {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
									   0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
									   1,
									   7, {0, 1, 2, 3, 4, 5, 6}, 5, {1, 3, 5, 7, 11}
									   *Ogrkkes,
									   *Otvegas, *Otvannas, *Otzommas,
									   *Otrfrhos, *Otdcrhos,
									   0, {5, 1.5}, {10000, 100, 10}, {24, 6},
									   0.02, 0.01)
*/

{
	if(matu[0] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *insval;
	double *inxval;
	double *psval;

 
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inxval = new double[2];
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];

	if(xvaltype == 100) {
		inxval[0] = xval[0]*bval[0];
		inxval[1] = xval[1]*bval[1];
	}
	else {
		inxval[0] = xval[0];		 
		inxval[1] = xval[1];
	}

	psval = new double[uAssetN*maxevalN];
	for(i=0; i<uAssetN*maxevalN; i++)
		psval[i] = sval[i];

	TwoAssetDigital *twoassetdigital = new TwoAssetDigital(cpFlag, insval, bval, inxval, Cpn[0], 
	      matu[0], ctime,
		  evalN[0], psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, SmaxlevelC, SmeshM, TmeshN,
		  shiftRate, spreadRate);
			  
	/* value & 기본 그릭 산출 */    
	value = twoassetdigital->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = twoassetdigital->m_Ogreeks[i];
	}

	int vegaN;
	vegaN = termgrNs[0];
	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = twoassetdigital->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = twoassetdigital->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = twoassetdigital->m_Otvegas1[i];
			Otvannas1[i] = twoassetdigital->m_Otvannas1[i];
			Otzommas1[i] = twoassetdigital->m_Otzommas1[i];
						
			Otvegas2[i] = twoassetdigital->m_Otvegas2[i];
			Otvannas2[i] = twoassetdigital->m_Otvannas2[i];
			Otzommas2[i] = twoassetdigital->m_Otzommas2[i];
		}
	}


	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = twoassetdigital->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = twoassetdigital->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = twoassetdigital->m_Otrfrhos[i];
			Otdcrhos[i] = twoassetdigital->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)

	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = twoassetdigital->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = twoassetdigital->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = twoassetdigital->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 4)
	  
 
		
	delete insval;
	delete psval;
	delete inxval;
	 

	delete twoassetdigital;

	return value;
		
} 
 

 double ODYSSEY_API  _stdcall StepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
	}
 
 
	StepDownB *stepdownb = new StepDownB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdownb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
	                                     // 2: stepdown ko (원금비보장 ko)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdownb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdownb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdownb->m_Otvegas1[i];
			Otvannas1[i] = stepdownb->m_Otvannas1[i];
			Otzommas1[i] = stepdownb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdownb->m_Otvegas2[i];
			Otvannas2[i] = stepdownb->m_Otvannas2[i];
			Otzommas2[i] = stepdownb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdownb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdownb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdownb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdownb;

	return value;
		
} 
  



 





 double ODYSSEY_API  _stdcall CliquetAFn(int uAssetN, double *sval, double *bval, double prate, double GFloor, double LCap, double LFloor, double NAQ, 	                                 
	                                  int matuN, int *matu, double ctime,									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag,  char *kCode,									 
									  double Xmaxlevel, double Zmesh, int XmeshM, int *TmeshN)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			XmeshM = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	 
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	 

	CliquetA *cliqueta = new CliquetA(sval, bval, prate, GFloor, LCap, LFloor,NAQ,  		  
	      matuN, matu, ctime,		   
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Xmaxlevel, Zmesh, XmeshM, TmeshN);
 
			  
	/* value & 기본 그릭 산출 */    
	value = cliqueta->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = cliqueta->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = cliqueta->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = cliqueta->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = cliqueta->m_Otvegas[i];
			Otvannas[i] = cliqueta->m_Otvannas[i];
			Otzommas[i] = cliqueta->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = cliqueta->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = cliqueta->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = cliqueta->m_Otrfrhos[i];
			Otdcrhos[i] = cliqueta->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
	 
	delete cliqueta;

	return value;
		
} 
  

 



 double ODYSSEY_API  _stdcall ZeroBondFn(int matuN, int *matu, 									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,									 
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos)
{
	if(matu[0] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			iDCrate[0] = iDCrate[0]*0.8;
			 
		}
	}

	double value;	
   
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
 
   

	ZeroBond *zerobond = new ZeroBond(matu[0],
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate);
			  
	/* value & 기본 그릭 산출 */    
	value = zerobond->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = zerobond->m_Ogreeks[i];
	}

	 
	//  Vega 관련 그릭 산출 
	if(greektype == 0 || greektype == 2) {		
		/*
		value = zerobond->getVega(vegaN, vegaidx);
		
		//  Sum vegas 할당 
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = zerobond->m_Ogreeks[i];
		}

		//  term별 vegas 할당 
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = zerobond->m_Otvegas[i];
			Otvannas[i] = zerobond->m_Otvannas[i];
			Otzommas[i] = zerobond->m_Otzommas[i];
		}
		*/
			
		value = 0;
		//  Sum vegas 할당 
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = 0;
		}
		//  term별 vegas 할당 
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = 0;
			Otvannas[i] = 0;
			Otzommas[i] = 0;
		}
	}
	 

	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = zerobond->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = zerobond->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = zerobond->m_Otrfrhos[i];
			Otdcrhos[i] = zerobond->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
	 
	delete zerobond;

	return value;
		
} 


/*
(int *uAssetN, int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
*/
 double ODYSSEY_API  _stdcall FRNFn(double ctime,									 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,									 
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos)
									 
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(PayBDay[CDN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			PayBDay[CDN-1] = 50000;
		 
		}
	}

	double value;	
	 
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	 

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	FRN *frn = new FRN(ctime,		
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate);

	frn->m_modeltype = 1;
			  
	/* value & 기본 그릭 산출 */    
	value = frn->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = frn->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = 0; // frn->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = 0; //frn->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = 0;//frn->m_Otvegas[i];
			Otvannas[i] = 0;//frn->m_Otvannas[i];
			Otzommas[i] = 0;//frn->m_Otzommas[i];
		}


	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = frn->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = frn->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = frn->m_Otrfrhos[i];
			Otdcrhos[i] = frn->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 
 
	delete frn;

	return value;
		
} 
  





 double ODYSSEY_API  _stdcall StepDownCDAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double inkihval;
	double *psval;
	int i;

	double inkohval;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];
		inkohval = hval[1]*bval[0];

	}
	else {
		for(i=0; i<matuN; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
		inkohval = hval[1];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	StepDownCDA *stepdowncda = new StepDownCDA(knockFlag[0], knockFlag[1], sval[0], bval[0], inxval, inkihval, inkohval,
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	stepdowncda->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncda->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = stepdowncda->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncda->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncda->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = stepdowncda->m_Otvegas[i];
			Otvannas[i] = stepdowncda->m_Otvannas[i];
			Otzommas[i] = stepdowncda->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncda->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = stepdowncda->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncda->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncda->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete stepdowncda;

	return value;
		
} 
  



 double ODYSSEY_API  _stdcall StepDownCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; //ko
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // ko
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN[0]];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0]; //ko
		inkohval[1] = hval[3]*bval[1];
	}
	else {
		for(i=0; i<matuN[0]; i++) {		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2]; // ko
		inkohval[1] = hval[3];
	}
 
 
	StepDownCDB *stepdowncdb = new StepDownCDB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,	
		  CDlegInfo,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdowncdb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncdb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdowncdb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncdb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncdb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdowncdb->m_Otvegas1[i];
			Otvannas1[i] = stepdowncdb->m_Otvannas1[i];
			Otzommas1[i] = stepdowncdb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdowncdb->m_Otvegas2[i];
			Otvannas2[i] = stepdowncdb->m_Otvannas2[i];
			Otzommas2[i] = stepdowncdb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncdb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncdb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncdb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncdb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncdb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete psval;

	delete inkohval; //ko

	delete stepdowncdb;

	return value;
		
} 
  



//For MC
 int ODYSSEY_API  _stdcall odysseyRandom_INIT(long SimulNo) //, long iMAXAST)
 {
	 if(RandUp == 0) {
	 opr = odysseyRandom::getInstance();

	 opr->SimulNo = SimulNo;
	 opr->randseq(0,SimulNo,MAXDAYS,NULL);
	 opr->randnum(0,SimulNo,MAXAST,NULL,7);
	 //opr->randnum(0,SimulNo,iMAXAST,NULL,7);
	 RandUp = 1;
	 }
	 
	 return 1;
 }

 int ODYSSEY_API  _stdcall odysseyRandom_END(long SimulNo)//, long iMAXAST)
 {
	 opr->randseq(2,SimulNo,MAXDAYS,NULL);
	 opr->randnum(2,SimulNo,MAXAST,NULL,7);
	 //opr->randnum(2,SimulNo,iMAXAST,NULL,7);

	 opr->removeInstance();
	 RandUp = 0;

	 return 0;
 }


 
  
  


 //knockFlag -> [0]:kiFlag, [1]:koFlag, [2]:ki1Flag, [3]:ki2Flag
 double ODYSSEY_API  _stdcall HStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
	}
 
 
	HStepDownB *hstepdownb = new HStepDownB(knockFlag[0], knockFlag[1], knockFlag[2], knockFlag[3], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	hstepdownb->m_modeltype = uAssetN[1];  
			  
	/* value & 기본 그릭 산출 */    
	value = hstepdownb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = hstepdownb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = hstepdownb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = hstepdownb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = hstepdownb->m_Otvegas1[i];
			Otvannas1[i] = hstepdownb->m_Otvannas1[i];
			Otzommas1[i] = hstepdownb->m_Otzommas1[i];
						
			Otvegas2[i] = hstepdownb->m_Otvegas2[i];
			Otvannas2[i] = hstepdownb->m_Otvannas2[i];
			Otzommas2[i] = hstepdownb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = hstepdownb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = hstepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = hstepdownb->m_Otrfrhos[i];
			Otdcrhos[i] = hstepdownb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = hstepdownb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = hstepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = hstepdownb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete hstepdownb;

	return value;
		
} 
  



 double ODYSSEY_API  _stdcall RStepUpAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag,  char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double inkihval;
	double *psval;
	int i;

	double inkohval;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];
		inkohval = hval[1]*bval[0];

	}
	else {
		for(i=0; i<matuN; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
		inkohval = hval[1];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 
	RStepUpA *rstepupa = new RStepUpA(knockFlag[0],knockFlag[1], sval[0], bval[0], inxval, inkihval, inkohval,
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	rstepupa->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = rstepupa->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = rstepupa->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = rstepupa->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = rstepupa->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = rstepupa->m_Otvegas[i];
			Otvannas[i] = rstepupa->m_Otvannas[i];
			Otzommas[i] = rstepupa->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = rstepupa->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = rstepupa->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = rstepupa->m_Otrfrhos[i];
			Otdcrhos[i] = rstepupa->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete rstepupa;

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall StepDownCPAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag,  char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN[0]+matuN[1]-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
 
	double inkihval;
	double *psval;
	int i;

	double inkohval;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN[0]+matuN[1]];
 
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]+matuN[1]; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];
		inkohval = hval[1]*bval[0];

	}
	else {		 
		for(i=0; i<matuN[0]+matuN[1]; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
		inkohval = hval[1];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 
	StepDownCPA *stepdowncpa = new StepDownCPA(knockFlag[0],knockFlag[1], sval[0], bval[0], inxval, inkihval, inkohval,
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	stepdowncpa->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpa->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = stepdowncpa->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncpa->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpa->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = stepdowncpa->m_Otvegas[i];
			Otvannas[i] = stepdowncpa->m_Otvannas[i];
			Otzommas[i] = stepdowncpa->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncpa->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = stepdowncpa->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncpa->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncpa->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete stepdowncpa;

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall RStepUpCDAFn(int *uAssetN, int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double inkihval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];

	}
	else {
		for(i=0; i<matuN; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	RStepUpCDA *rstepupcda = new RStepUpCDA(knockFlag, sval[0], bval[0], inxval, inkihval, 
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	rstepupcda->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = rstepupcda->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = rstepupcda->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = rstepupcda->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = rstepupcda->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = rstepupcda->m_Otvegas[i];
			Otvannas[i] = rstepupcda->m_Otvannas[i];
			Otzommas[i] = rstepupcda->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = rstepupcda->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = rstepupcda->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = rstepupcda->m_Otrfrhos[i];
			Otdcrhos[i] = rstepupcda->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete rstepupcda;

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall StepDownCPCDAFn(int *uAssetN, int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN[0]+matuN[1]-1]  = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double inkihval;
	double *psval;
	int i;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN[0]+matuN[1]];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]+matuN[1]; i++) 		
			inxval[i] = xval[i]*bval[0];

		inkihval = hval[0]*bval[0];

	}
	else {
		for(i=0; i<matuN[0]+matuN[1]; i++)		
			inxval[i] = xval[i];

		inkihval = hval[0];
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	StepDownCPCDA *stepdowncpcda = new StepDownCPCDA(knockFlag, sval[0], bval[0], inxval, inkihval, 
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	stepdowncpcda->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpcda->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = stepdowncpcda->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncpcda->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpcda->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = stepdowncpcda->m_Otvegas[i];
			Otvannas[i] = stepdowncpcda->m_Otvannas[i];
			Otzommas[i] = stepdowncpcda->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncpcda->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = stepdowncpcda->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncpcda->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncpcda->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete psval;
	delete stepdowncpcda;

	return value;
		
} 
  



 double ODYSSEY_API  _stdcall StepDownCPBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
									   // xval : paycpnbarrier, xval1, xval2 순서임
									   // Cpn : paycpn, cpn 순서
									   // matuN : paycpnN, matuN 순서
									   // matu : paymatu, matu 순서

									   
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]+matuN[1]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*(matuN[0]+matuN[1])];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
		}
		for(i=0; i<matuN[1]; i++) { 		
			inxval[uAssetN[0]*matuN[0]+i] = xval[uAssetN[0]*matuN[0]+i]*bval[0];
			inxval[uAssetN[0]*matuN[0]+ matuN[1]+i] = xval[uAssetN[0]*matuN[0]+matuN[1]+i]*bval[1];
		}

		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
	}
	else {				
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
		}
		for(i=0; i<matuN[1]; i++) { 		
			inxval[uAssetN[0]*matuN[0]+i] = xval[uAssetN[0]*matuN[0]+i];
			inxval[uAssetN[0]*matuN[0]+ matuN[1]+i] = xval[uAssetN[0]*matuN[0]+matuN[1]+i];
		}

		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
	}
 
 
	StepDownCPB *stepdowncpb = new StepDownCPB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdowncpb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
	                                     // 2: stepdown ko (원금비보장 ko)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdowncpb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncpb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdowncpb->m_Otvegas1[i];
			Otvannas1[i] = stepdowncpb->m_Otvannas1[i];
			Otzommas1[i] = stepdowncpb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdowncpb->m_Otvegas2[i];
			Otvannas2[i] = stepdowncpb->m_Otvannas2[i];
			Otzommas2[i] = stepdowncpb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncpb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncpb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncpb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncpb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncpb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdowncpb;

	return value;
		
} 
  


 double ODYSSEY_API  _stdcall StepDownCPCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]+matuN[1]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval;
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2];
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*(matuN[0]+matuN[1])];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
		}
		for(i=0; i<matuN[1]; i++) { 		
			inxval[uAssetN[0]*matuN[0]+i] = xval[uAssetN[0]*matuN[0]+i]*bval[0];
			inxval[uAssetN[0]*matuN[0]+ matuN[1]+i] = xval[uAssetN[0]*matuN[0]+matuN[1]+i]*bval[1];
		}

		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];
		
		inkohval[0] = hval[2]*bval[0];
		inkohval[1] = hval[3]*bval[1];
	}
	else {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
		}
		for(i=0; i<matuN[1]; i++) { 		
			inxval[uAssetN[0]*matuN[0]+i] = xval[uAssetN[0]*matuN[0]+i];
			inxval[uAssetN[0]*matuN[0]+ matuN[1]+i] = xval[uAssetN[0]*matuN[0]+matuN[1]+i];
		}

		inkihval[0] = hval[0];
		inkihval[1] = hval[1];
		
		inkohval[0] = hval[2];
		inkohval[1] = hval[3];
	}
 
 
	StepDownCPCDB *stepdowncpcdb = new StepDownCPCDB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,	
		  CDlegInfo,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdowncpcdb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpcdb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdowncpcdb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdowncpcdb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpcdb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdowncpcdb->m_Otvegas1[i];
			Otvannas1[i] = stepdowncpcdb->m_Otvannas1[i];
			Otzommas1[i] = stepdowncpcdb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdowncpcdb->m_Otvegas2[i];
			Otvannas2[i] = stepdowncpcdb->m_Otvannas2[i];
			Otzommas2[i] = stepdowncpcdb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdowncpcdb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpcdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdowncpcdb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdowncpcdb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncpcdb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpcdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncpcdb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdowncpcdb;

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall StepDownCTBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	double *inCThval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
	inCThval = new double[2];

		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd

		inCThval[0] = hval[4]*bval[0];
		inCThval[1] = hval[5]*bval[1];
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd

		inCThval[0] = hval[4];
		inCThval[1] = hval[5];
	}
 
 
	StepDownCTB *stepdownctb = new StepDownCTB(knockFlag[0], knockFlag[1], knockFlag[2], insval,  bval, inxval, inkihval, inkohval, inCThval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdownctb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
	                                     
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownctb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdownctb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdownctb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownctb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdownctb->m_Otvegas1[i];
			Otvannas1[i] = stepdownctb->m_Otvannas1[i];
			Otzommas1[i] = stepdownctb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdownctb->m_Otvegas2[i];
			Otvannas2[i] = stepdownctb->m_Otvannas2[i];
			Otzommas2[i] = stepdownctb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdownctb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownctb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdownctb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdownctb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownctb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownctb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownctb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete inCThval;
	delete psval;

	delete stepdownctb;

	return value;
		
} 
  





 double ODYSSEY_API  _stdcall StepDownCTCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; //for ko : 현재는 CTCD에서는 사용하지 X
	double *inCThval; //for CT
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2];
	inCThval = new double[2];
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN[0]];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		
		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd

		inCThval[0] = hval[4]*bval[0];
		inCThval[1] = hval[5]*bval[1];

	}
	else {
		for(i=0; i<matuN[0]; i++) {		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		
		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd

		inCThval[0] = hval[4];
		inCThval[1] = hval[5];
	}
 
    //matu = matu + obbday+paybday + CTmatu
	StepDownCTCDB *stepdownctcdb = new StepDownCTCDB(knockFlag[0],knockFlag[1],knockFlag[2], insval,  bval, inxval, inkihval, inkohval, inCThval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,	
		  CDlegInfo,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdownctcdb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownctcdb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdownctcdb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdownctcdb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownctcdb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdownctcdb->m_Otvegas1[i];
			Otvannas1[i] = stepdownctcdb->m_Otvannas1[i];
			Otzommas1[i] = stepdownctcdb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdownctcdb->m_Otvegas2[i];
			Otvannas2[i] = stepdownctcdb->m_Otvannas2[i];
			Otzommas2[i] = stepdownctcdb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdownctcdb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownctcdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdownctcdb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdownctcdb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownctcdb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownctcdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownctcdb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete inCThval;

	delete psval;

	delete stepdownctcdb;

	return value;
		
} 
  

 
 


 double ODYSSEY_API  _stdcall StepDownBCBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[2*uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<2*matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[2*matuN+i] = xval[2*matuN+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
	}
	else {
		for(i=0; i<2*matuN; i++) {		
			inxval[i] = xval[i];
			inxval[2*matuN+i] = xval[2*matuN+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
	}
 
 
	StepDownBCB *stepdownbcb = new StepDownBCB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdownbcb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
	                                     // 2: stepdown ko (원금비보장 ko)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownbcb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdownbcb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdownbcb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownbcb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdownbcb->m_Otvegas1[i];
			Otvannas1[i] = stepdownbcb->m_Otvannas1[i];
			Otzommas1[i] = stepdownbcb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdownbcb->m_Otvegas2[i];
			Otvannas2[i] = stepdownbcb->m_Otvannas2[i];
			Otzommas2[i] = stepdownbcb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdownbcb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownbcb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdownbcb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdownbcb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownbcb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownbcb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownbcb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdownbcb;

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall StepDownBCCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; //ko
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // ko
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[2*uAssetN[0]*matuN[0]];
	if(xvaltype == 100) {
		for(i=0; i<2*matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[2*matuN[0]+i] = xval[2*matuN[0]+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0]; //ko
		inkohval[1] = hval[3]*bval[1];
	}
	else {
		for(i=0; i<2*matuN[0]; i++) {		
			inxval[i] = xval[i];
			inxval[2*matuN[0]+i] = xval[2*matuN[0]+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2]; // ko
		inkohval[1] = hval[3];
	}
 
 
	StepDownBCCDB *stepdownbccdb = new StepDownBCCDB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,	
		  CDlegInfo,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	stepdownbccdb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownbccdb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = stepdownbccdb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = stepdownbccdb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownbccdb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = stepdownbccdb->m_Otvegas1[i];
			Otvannas1[i] = stepdownbccdb->m_Otvannas1[i];
			Otzommas1[i] = stepdownbccdb->m_Otzommas1[i];
						
			Otvegas2[i] = stepdownbccdb->m_Otvegas2[i];
			Otvannas2[i] = stepdownbccdb->m_Otvannas2[i];
			Otzommas2[i] = stepdownbccdb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = stepdownbccdb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownbccdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = stepdownbccdb->m_Otrfrhos[i];
			Otdcrhos[i] = stepdownbccdb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownbccdb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownbccdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownbccdb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete psval;

	delete inkohval; //ko

	delete stepdownbccdb;

	return value;
		
} 
  
 


 double ODYSSEY_API  _stdcall StabilityAFn(int *uAssetN, int knockFlag, double *sval, double erPerf, double *xval, double multiplier, 
	                                  double *CpnFactor, double *Cpn,
	                                  int matuN, int *matu, double ctime,									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, 
									  double *Otrfrhos, double *Otdcrhos,
									  double SLimit, long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
 	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	 

	int i;

		  
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	StabilityA *stabilitya = new StabilityA(knockFlag, sval[0], sval[1], erPerf, xval[0], multiplier, 
		  CpnFactor, Cpn,
	      matuN, matu, ctime,		   			  
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  SLimit, SimulNo);

	 
	stabilitya->m_modeltype = uAssetN[1];  
			  
	for(i=0; i<MAXGREEKS1; i++)
		stabilitya->m_Ogreeks[i] = 0;

	/* value & 기본 그릭 산출 */    
	value = stabilitya->getValue();
	Ogreeks[0] = stabilitya->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 4) {// 0: Total, 1: P,  2:V, 3: R  4:D&G & Theta,
		value = stabilitya->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN1; i++)
			Ogreeks[i] = stabilitya->m_Ogreeks[i];
	}

	 


	/* Vega 관련 그릭 산출 */
 
	if(greektype == 0 || greektype == 2) {		 
		value = stabilitya->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = stabilitya->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		 
		for(i=0; i<vegaN; i++) 
			Otvegas[i] = stabilitya->m_Otvegas[i];
			
			
	}

	 
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
	 
		 
		value = stabilitya->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = stabilitya->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
	 			
		for(i=0; i<rhoN; i++) 
			Otrfrhos[i] = stabilitya->m_Otrfrhos[i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = stabilitya->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = stabilitya->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stabilitya->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stabilitya->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  
	 
	delete stabilitya;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
		
} 
  




 double ODYSSEY_API  _stdcall AirbagStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1,  
									  double *Otvegas2,  
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  long SimulNo)
									 

{
	if(matu[matuN[0]-1] < 0)   //matuN[0] : matuN , matuN[1]: airbagN
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;
		 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;
	double *inairbagxval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
	inairbagxval = new double[2];
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN[0]];
	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd

		inairbagxval[0] = hval[4]*bval[0];
		inairbagxval[1] = hval[5]*bval[1];
	}
	else {
		for(i=0; i<matuN[0]; i++) {		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd

		inairbagxval[0] = hval[4];
		inairbagxval[1] = hval[5];
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}
 
	AirbagStepDownB *airbagstepdownb = new AirbagStepDownB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, inairbagxval, Cpn,
	      matuN[0], matuN[1],  matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	airbagstepdownb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장) 다이나믹 참여율
	                  

	for(i=0; i<MAXGREEKS2; i++) 
		airbagstepdownb->m_Ogreeks[i] = 0;
			  
	/* value & 기본 그릭 산출 */    
	value = airbagstepdownb->getValue();	 
	Ogreeks[twoPriceidx] = airbagstepdownb->m_Ogreeks[twoPriceidx];
	 

		
	if(greektype == 0 || greektype == 5) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = airbagstepdownb->getDeltaGammaTheta();	 
			
		Ogreeks[twoUpDeltaCidx1] = airbagstepdownb->m_Ogreeks[twoUpDeltaCidx1];
		Ogreeks[twoDownDeltaCidx1] = airbagstepdownb->m_Ogreeks[twoDownDeltaCidx1];
		Ogreeks[twoUpDeltaCidx2] = airbagstepdownb->m_Ogreeks[twoUpDeltaCidx2];
		Ogreeks[twoDownDeltaCidx2] = airbagstepdownb->m_Ogreeks[twoDownDeltaCidx2];

		Ogreeks[twoUpGammaCidx1] = airbagstepdownb->m_Ogreeks[twoUpGammaCidx1];
		Ogreeks[twoDownGammaCidx1] = airbagstepdownb->m_Ogreeks[twoDownGammaCidx1];
		Ogreeks[twoUpGammaCidx2] = airbagstepdownb->m_Ogreeks[twoUpGammaCidx2];
		Ogreeks[twoDownGammaCidx2] = airbagstepdownb->m_Ogreeks[twoDownGammaCidx2];

		Ogreeks[twoThetaidx] = airbagstepdownb->m_Ogreeks[twoThetaidx];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = airbagstepdownb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */ 
		Ogreeks[twoVegaidx1] = airbagstepdownb->m_Ogreeks[twoVegaidx1];
		Ogreeks[twoVegaidx2] = airbagstepdownb->m_Ogreeks[twoVegaidx2];

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = airbagstepdownb->m_Otvegas1[i];						
			Otvegas2[i] = airbagstepdownb->m_Otvegas2[i];	 
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = airbagstepdownb->getRho(rhoN, rhoidx);
 
		/* Sum rhos 할당 */ 
		Ogreeks[twoRhorfidx] = airbagstepdownb->m_Ogreeks[twoRhorfidx];
		Ogreeks[twoRhodcidx] = airbagstepdownb->m_Ogreeks[twoRhodcidx];

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = airbagstepdownb->m_Otrfrhos[i];
			Otdcrhos[i] = airbagstepdownb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 

	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	//CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = airbagstepdownb->getCorrDelta(corrdeltaN, corrdeltaidx);


		 //Sum CorrealtionDelta 할당
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = airbagstepdownb->m_Ogreeks[i]; 
		}

		// term별 CorrealtionDelta 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = airbagstepdownb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	*/

	  

  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete inairbagxval;
	delete psval;

	delete airbagstepdownb;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 





 double ODYSSEY_API  _stdcall LizardStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1,  
									  double *Otvegas2,  
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  long SimulNo)
									 

{
	if(matu[matuN-1] < 0)   //matuN[0] : matuN , matuN[1]: airbagN
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
		 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;
	 
	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	inkihval = new double[2];
	inkohval = new double[2]; // koadd
	 
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd

		 
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		}
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
 
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}
 
	LizardStepDownB *lizardstepdownb = new LizardStepDownB(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval,   Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	lizardstepdownb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장) 다이나믹 참여율
	                  

	for(i=0; i<MAXGREEKS2; i++) 
		lizardstepdownb->m_Ogreeks[i] = 0;
			  
	/* value & 기본 그릭 산출 */    
	value = lizardstepdownb->getValue();	 
	Ogreeks[twoPriceidx] = lizardstepdownb->m_Ogreeks[twoPriceidx];
	 

		
	if(greektype == 0 || greektype == 5) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = lizardstepdownb->getDeltaGammaTheta();	 
			
		Ogreeks[twoUpDeltaCidx1] = lizardstepdownb->m_Ogreeks[twoUpDeltaCidx1];
		Ogreeks[twoDownDeltaCidx1] = lizardstepdownb->m_Ogreeks[twoDownDeltaCidx1];
		Ogreeks[twoUpDeltaCidx2] = lizardstepdownb->m_Ogreeks[twoUpDeltaCidx2];
		Ogreeks[twoDownDeltaCidx2] = lizardstepdownb->m_Ogreeks[twoDownDeltaCidx2];

		Ogreeks[twoUpGammaCidx1] = lizardstepdownb->m_Ogreeks[twoUpGammaCidx1];
		Ogreeks[twoDownGammaCidx1] = lizardstepdownb->m_Ogreeks[twoDownGammaCidx1];
		Ogreeks[twoUpGammaCidx2] = lizardstepdownb->m_Ogreeks[twoUpGammaCidx2];
		Ogreeks[twoDownGammaCidx2] = lizardstepdownb->m_Ogreeks[twoDownGammaCidx2];

		Ogreeks[twoThetaidx] = lizardstepdownb->m_Ogreeks[twoThetaidx];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = lizardstepdownb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */ 
		Ogreeks[twoVegaidx1] = lizardstepdownb->m_Ogreeks[twoVegaidx1];
		Ogreeks[twoVegaidx2] = lizardstepdownb->m_Ogreeks[twoVegaidx2];

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = lizardstepdownb->m_Otvegas1[i];						
			Otvegas2[i] = lizardstepdownb->m_Otvegas2[i];	 
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = lizardstepdownb->getRho(rhoN, rhoidx);
 
		/* Sum rhos 할당 */ 
		Ogreeks[twoRhorfidx] = lizardstepdownb->m_Ogreeks[twoRhorfidx];
		Ogreeks[twoRhodcidx] = lizardstepdownb->m_Ogreeks[twoRhodcidx];

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = lizardstepdownb->m_Otrfrhos[i];
			Otdcrhos[i] = lizardstepdownb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 

	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	//CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = lizardstepdownb->getCorrDelta(corrdeltaN, corrdeltaidx);


		 //Sum CorrealtionDelta 할당
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = lizardstepdownb->m_Ogreeks[i]; 
		}

		// term별 CorrealtionDelta 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = lizardstepdownb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	*/

	  

  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	 
	delete psval;

	delete lizardstepdownb;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 




 double ODYSSEY_API  _stdcall BarrierBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	 
	double *inkohval; // koadd
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	 
	inkohval = new double[2]; // koadd
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}
	 

		inkohval[0] = hval[0]*bval[0];  // koadd
		inkohval[1] = hval[1]*bval[1];  // koadd
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		}
		 

		inkohval[0] = hval[0];   // koadd
		inkohval[1] = hval[1];   // koadd
	}
 
 
	BarrierB *barrierb = new BarrierB(knockFlag[0], insval,  bval, inxval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	barrierb->m_modeltype = uAssetN[1]; // 0: UOC with 만기 rebate만 구현되어 있음!
	                                    
			  
	/* value & 기본 그릭 산출 */    
	value = barrierb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = barrierb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = barrierb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = barrierb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = barrierb->m_Otvegas1[i];
			Otvannas1[i] = barrierb->m_Otvannas1[i];
			Otzommas1[i] = barrierb->m_Otzommas1[i];
						
			Otvegas2[i] = barrierb->m_Otvegas2[i];
			Otvannas2[i] = barrierb->m_Otvannas2[i];
			Otzommas2[i] = barrierb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
	 
		 
		value = barrierb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = barrierb->m_Ogreeks[i]; 
		}	 
		 
		


		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = barrierb->m_Otrfrhos[i];
			Otdcrhos[i] = barrierb->m_Otdcrhos[i];
		}
		 
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = barrierb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = barrierb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = barrierb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
 
	delete inkohval;
	delete psval;

	delete barrierb;

	return value;
		
} 
  
 
 double ODYSSEY_API  _stdcall PlainBFn(int *uAssetN,  double *sval, double *bval, double *xval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	  
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	 
	 
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	 
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];
		}		 
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];
		} 
	}
 
 
	PlainB *plainb = new PlainB(insval,  bval, inxval,   Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	plainb->m_modeltype = uAssetN[1]; // 0: UOC with 만기 rebate만 구현되어 있음!
	                                    
			  
	/* value & 기본 그릭 산출 */    
	value = plainb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = plainb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = plainb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = plainb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = plainb->m_Otvegas1[i];
			Otvannas1[i] = plainb->m_Otvannas1[i];
			Otzommas1[i] = plainb->m_Otzommas1[i];
						
			Otvegas2[i] = plainb->m_Otvegas2[i];
			Otvannas2[i] = plainb->m_Otvannas2[i];
			Otzommas2[i] = plainb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
	 
		 
		value = plainb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = plainb->m_Ogreeks[i]; 
		}	 
		 
		


		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = plainb->m_Otrfrhos[i];
			Otdcrhos[i] = plainb->m_Otdcrhos[i];
		}
		 
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = plainb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = plainb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = plainb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
  
	delete psval;

	delete plainb;

	return value;
		
} 
  

//3S

  double ODYSSEY_API  _stdcall StepDownCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int *irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,
									  int *corrtype, int *corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos,long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i,k;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[uAssetN[0]];
	inkihval = new double[uAssetN[0]];
	inkohval = new double[uAssetN[0]]; // koadd
		
	for(k=0; k<uAssetN[0]; k++)
		insval[k] = sval[maxevalN*k];
	/*
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	*/

	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];
	

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN; i++) { 		
				inxval[k*matuN+i] = xval[k*matuN+i]*bval[k];
				//inxval[matuN+i] = xval[matuN+i]*bval[1];
			}
			inkihval[k] = hval[k]*bval[k];

			inkohval[k] = hval[uAssetN[0] + k]*bval[k];
		}
		/*
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
		*/
	}
	else {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN; i++) {		
				inxval[k*matuN+i] = xval[k*matuN+i];
				//inxval[matuN+i] = xval[matuN+i];
			}
			inkihval[k] = hval[k];
			inkohval[k] = hval[uAssetN[0]+k];
		}
		/*
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
		*/
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	StepDownC *stepdownc = new StepDownC(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	stepdownc->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 10  stepdown (원금비보장 noKI)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdownc->getValue();
	Ogreeks[0] = stepdownc->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 2) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = stepdownc->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN3; i++)
			Ogreeks[i] = stepdownc->m_Ogreeks[i];
	}

	 


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 3) {		 
		value = stepdownc->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN3; i<BasicGRsN3+VegaGRsN3; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdownc->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(k=0; k<uAssetN[0]; k++)
			for(i=0; i<vegaN; i++) 
				Otvegas[vegaN*k+i] = stepdownc->m_Otvegas[vegaN*k+i];
			
			
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		for(k=0; k<=uAssetN[0]; k++)		
			stepdownc->m_iRateSeq[k] = irateN[k+1]; // 나라별 금리 코드 그냥 구별되게 0, 1 ..
		 
		 
		value = stepdownc->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN3+VegaGRsN3; i<BasicGRsN3+VegaGRsN3+RhoGRsN3; i++) { // 10, 11
			Ogreeks[i] = stepdownc->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
		for(k=0; k<=uAssetN[0]; k++)		
			for(i=0; i<rhoN; i++) 
			Otrfrhos[rhoN*k+i] = stepdownc->m_Otrfrhos[rhoN*k+i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = stepdownc->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = stepdownc->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdownc->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdownc->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdownc;
	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 
  




  double ODYSSEY_API  _stdcall StepDownCDCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int *irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,
									  int *corrtype, int *corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos,long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i,k;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[uAssetN[0]];
	inkihval = new double[uAssetN[0]];
	inkohval = new double[uAssetN[0]]; // koadd
		
	for(k=0; k<uAssetN[0]; k++)
		insval[k] = sval[maxevalN*k];
	/*
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	*/

	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];
	

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN[0]];
	if(xvaltype == 100) {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) { 		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i]*bval[k];
				//inxval[matuN+i] = xval[matuN+i]*bval[1];
			}
			inkihval[k] = hval[k]*bval[k];

			inkohval[k] = hval[uAssetN[0] + k]*bval[k];
		}
		/*
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
		*/
	}
	else {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) {		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i];
				//inxval[matuN+i] = xval[matuN+i];
			}
			inkihval[k] = hval[k];
			inkohval[k] = hval[uAssetN[0]+k];
		}
		/*
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
		*/
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	StepDownCDC *stepdowncdc = new StepDownCDC(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDlegInfo,
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	stepdowncdc->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 10  stepdown (원금비보장 noKI)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncdc->getValue();
	Ogreeks[0] = stepdowncdc->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 2) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = stepdowncdc->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN3; i++)
			Ogreeks[i] = stepdowncdc->m_Ogreeks[i];
	}

	 


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 3) {		 
		value = stepdowncdc->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN3; i<BasicGRsN3+VegaGRsN3; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncdc->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(k=0; k<uAssetN[0]; k++)
			for(i=0; i<vegaN; i++) 
				Otvegas[vegaN*k+i] = stepdowncdc->m_Otvegas[vegaN*k+i];
			
			
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		for(k=0; k<=uAssetN[0]; k++)		
			stepdowncdc->m_iRateSeq[k] = irateN[k+1]; // 나라별 금리 코드 그냥 구별되게 0, 1 ..
		 
		 
		value = stepdowncdc->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN3+VegaGRsN3; i<BasicGRsN3+VegaGRsN3+RhoGRsN3; i++) { // 10, 11
			Ogreeks[i] = stepdowncdc->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
		for(k=0; k<=uAssetN[0]; k++)		
			for(i=0; i<rhoN; i++) 
			Otrfrhos[rhoN*k+i] = stepdowncdc->m_Otrfrhos[rhoN*k+i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = stepdowncdc->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncdc->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncdc->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncdc->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdowncdc;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 
  



  double ODYSSEY_API  _stdcall StepDownCPCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int *irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,
									  int *corrtype, int *corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos,long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]+matuN[1]-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i,k;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[uAssetN[0]];
	inkihval = new double[uAssetN[0]];
	inkohval = new double[uAssetN[0]]; // koadd
		
	for(k=0; k<uAssetN[0]; k++)
		insval[k] = sval[maxevalN*k];
	/*
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	*/

	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];
	

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*(matuN[0]+matuN[1])];
	if(xvaltype == 100) {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) { 		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i]*bval[k];
				//inxval[matuN+i] = xval[matuN+i]*bval[1];
			}
			for(i=0; i<matuN[1]; i++) {
				inxval[uAssetN[0]*matuN[0] + k*matuN[1] + i] = xval[uAssetN[0]*matuN[0] + k*matuN[1] + i]*bval[k];
			}

			inkihval[k] = hval[k]*bval[k];

			inkohval[k] = hval[uAssetN[0] + k]*bval[k];
		}
		/*
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
		*/
	}
	else {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) {		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i];
				//inxval[matuN+i] = xval[matuN+i];
			}
			for(i=0; i<matuN[1]; i++) {
				inxval[uAssetN[0]*matuN[0] + k*matuN[1] + i] = xval[uAssetN[0]*matuN[0] + k*matuN[1] + i];
			}
			inkihval[k] = hval[k];
			inkohval[k] = hval[uAssetN[0]+k];
		}
		/*
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
		*/
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	StepDownCPC *stepdowncpc = new StepDownCPC(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	stepdowncpc->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 10  stepdown (원금비보장 noKI)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpc->getValue();
	Ogreeks[0] = stepdowncpc->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 2) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = stepdowncpc->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN3; i++)
			Ogreeks[i] = stepdowncpc->m_Ogreeks[i];
	}

	 


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 3) {		 
		value = stepdowncpc->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN3; i<BasicGRsN3+VegaGRsN3; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpc->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(k=0; k<uAssetN[0]; k++)
			for(i=0; i<vegaN; i++) 
				Otvegas[vegaN*k+i] = stepdowncpc->m_Otvegas[vegaN*k+i];
			
			
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		for(k=0; k<=uAssetN[0]; k++)		
			stepdowncpc->m_iRateSeq[k] = irateN[k+1]; // 나라별 금리 코드 그냥 구별되게 0, 1 ..
		 
		 
		value = stepdowncpc->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN3+VegaGRsN3; i<BasicGRsN3+VegaGRsN3+RhoGRsN3; i++) { // 10, 11
			Ogreeks[i] = stepdowncpc->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
		for(k=0; k<=uAssetN[0]; k++)		
			for(i=0; i<rhoN; i++) 
			Otrfrhos[rhoN*k+i] = stepdowncpc->m_Otrfrhos[rhoN*k+i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = stepdowncpc->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncpc->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpc->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncpc->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdowncpc;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 
  




  double ODYSSEY_API  _stdcall StepDownCPCDCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int *irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,
									  int *corrtype, int *corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos,long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]+matuN[1]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]+matuN[1]-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	double *inkohval; // koadd
	double *psval;

	int i,k;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[uAssetN[0]];
	inkihval = new double[uAssetN[0]];
	inkohval = new double[uAssetN[0]]; // koadd
		
	for(k=0; k<uAssetN[0]; k++)
		insval[k] = sval[maxevalN*k];
	/*
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	*/

	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];
	

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*(matuN[0]+matuN[1])];
	if(xvaltype == 100) {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) { 		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i]*bval[k];
				//inxval[matuN+i] = xval[matuN+i]*bval[1];
			}
			for(i=0; i<matuN[1]; i++) {
				inxval[uAssetN[0]*matuN[0] + k*matuN[1] + i] = xval[uAssetN[0]*matuN[0] + k*matuN[1] + i] * bval[k];
			}

			inkihval[k] = hval[k]*bval[k];

			inkohval[k] = hval[uAssetN[0] + k]*bval[k];
		}
		/*
		inkihval[0] = hval[0]*bval[0];
		inkihval[1] = hval[1]*bval[1];

		inkohval[0] = hval[2]*bval[0];  // koadd
		inkohval[1] = hval[3]*bval[1];  // koadd
		*/
	}
	else {
		for(k=0; k<uAssetN[0]; k++) {
			for(i=0; i<matuN[0]; i++) {		
				inxval[k*matuN[0]+i] = xval[k*matuN[0]+i];
				//inxval[matuN+i] = xval[matuN+i];
			}
			for(i=0; i<matuN[1]; i++) {
				inxval[uAssetN[0]*matuN[0] + k*matuN[1] + i] = xval[uAssetN[0]*matuN[0] + k*matuN[1] + i] * bval[k];
			}
			inkihval[k] = hval[k];
			inkohval[k] = hval[uAssetN[0]+k];
		}
		/*
		inkihval[0] = hval[0];
		inkihval[1] = hval[1];

		inkohval[0] = hval[2];   // koadd
		inkohval[1] = hval[3];   // koadd
		*/
	}
 
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	StepDownCPCDC *stepdowncpcdc = new StepDownCPCDC(knockFlag[0], knockFlag[1], insval,  bval, inxval, inkihval, inkohval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDlegInfo,
		  irateN[0], irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr, SimulNo);

	 
	stepdowncpcdc->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 10  stepdown (원금비보장 noKI)
			  
	/* value & 기본 그릭 산출 */    
	value = stepdowncpcdc->getValue();
	Ogreeks[0] = stepdowncpcdc->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 2) {// 0: Total, 1: P, 2:D&G & Theta, 3:V, 4: R
		value = stepdowncpcdc->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN3; i++)
			Ogreeks[i] = stepdowncpcdc->m_Ogreeks[i];
	}

	 


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 3) {		 
		value = stepdowncpcdc->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN3; i<BasicGRsN3+VegaGRsN3; i++) {  // 7, 8, 9
			Ogreeks[i] = stepdowncpcdc->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(k=0; k<uAssetN[0]; k++)
			for(i=0; i<vegaN; i++) 
				Otvegas[vegaN*k+i] = stepdowncpcdc->m_Otvegas[vegaN*k+i];
			
			
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		for(k=0; k<=uAssetN[0]; k++)		
			stepdowncpcdc->m_iRateSeq[k] = irateN[k+1]; // 나라별 금리 코드 그냥 구별되게 0, 1 ..
		 
		 
		value = stepdowncpcdc->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN3+VegaGRsN3; i<BasicGRsN3+VegaGRsN3+RhoGRsN3; i++) { // 10, 11
			Ogreeks[i] = stepdowncpcdc->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
		for(k=0; k<=uAssetN[0]; k++)		
			for(i=0; i<rhoN; i++) 
			Otrfrhos[rhoN*k+i] = stepdowncpcdc->m_Otrfrhos[rhoN*k+i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = stepdowncpcdc->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = stepdowncpcdc->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = stepdowncpcdc->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = stepdowncpcdc->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete inkohval;
	delete psval;

	delete stepdowncpcdc;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
} 
  




 double ODYSSEY_API  _stdcall TomAFn(int *uAssetN, double *sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,	                                   
	                                  int matuN, int *matu, double *prate, double ctime,									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, 
									  double *Otrfrhos, double *Otdcrhos,
									  long SimulNo,int batchFlag, char *kCode)
									   
 

{
 	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	 

	int i;

		  
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}
	
	int m_batchFlag;
	char *m_kCode;
	int outN, batchi;

	outN = 2*SoutRange + 1;
	 
	double *m_SrTemp;
	double *m_SrVTemp;
	 


	m_SrTemp = new double[outN+4];
	m_SrVTemp = new double[outN+4];
	 
	 
	
	m_batchFlag = batchFlag;
	m_kCode = new char[kCodeN+1];
	strcpy_s(m_kCode,kCodeN,kCode);

	TomA *toma = new TomA(sval[0], stom, btom, STDayN, STStartDay, IXT, IndexT, STEndDay, IXTheta, 		  
	      matuN, matu, prate, ctime,		   			  
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  SimulNo);

	 
	toma->m_modeltype = uAssetN[1];  
			  
	for(i=0; i<MAXGREEKS1; i++)
		toma->m_Ogreeks[i] = 0;

	if(m_batchFlag == 1) {
		for(batchi=0; batchi<outN+4; batchi++) {
			m_SrTemp[batchi] = sval[0]*(1-(SoutRange+2)*0.01 + 0.01*batchi);
			m_SrVTemp[batchi] = toma->getBatch(m_SrTemp[batchi]);
		}

		
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
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] << endl;
				}
				batchsf.close();
				
 				 
				delete fname;

	}
	

	/* value & 기본 그릭 산출 */    
	value = toma->getValue();
	Ogreeks[0] = toma->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 4) {// 0: Total, 1: P,  2:V, 3: R  4:D&G & Theta,
		value = toma->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN1; i++)
			Ogreeks[i] = toma->m_Ogreeks[i];
		
		if(m_batchFlag == 1) {	 
		
					using namespace std;
					using std::ofstream;

					ofstream batchsf;
			
					char *fname;
					fname = new char[kCodeN+1];

					strcpy_s(fname,kCodeN,m_kCode);
					strcat_s(fname,kCodeN,"_day+.txt");

					batchsf.open(fname,ios::out);
			
					batchsf << fixed;
					batchsf.precision(15);
					for(batchi=0; batchi<outN + 4; batchi++) {
						batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] + Ogreeks[Thetaidx] << endl;
					}
					batchsf.close();
					
 				 
				delete fname;
					 

		}
	}

	 //m_Ogreeks[Thetaidx]


	/* Vega 관련 그릭 산출 */
 
	if(greektype == 0 || greektype == 2) {		 
		value = toma->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = toma->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		 
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = toma->m_Otvegas[i];


		
			if(m_batchFlag) {

				using namespace std;
				using std::ofstream;

				ofstream batchsf;
			
				char *fname;
				char *tempstring;
				fname = new char[kCodeN+1];
				tempstring = new char[kCodeN+1];


				 
				strcpy_s(fname,kCodeN,m_kCode);			 			
				 
					strcat_s(fname,kCodeN,"_sigma+_");
				 
				_itoa_s(i,tempstring,kCodeN,10);			
				strcat_s(fname,kCodeN,tempstring);			
				strcat_s(fname,kCodeN,".txt");


		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] + Otvegas[i]<< endl;
				}
				batchsf.close();
				
				 
				strcpy_s(fname,kCodeN,m_kCode);			 			
				 
					strcat_s(fname,kCodeN,"_sigma-_");
				 		
				_itoa_s(i,tempstring,kCodeN,10);			
				strcat_s(fname,kCodeN,tempstring);			
				strcat_s(fname,kCodeN,".txt");


		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] - Otvegas[i]<< endl;
				}
				batchsf.close();

 				delete tempstring;
				delete fname;
			
			}// if(m_batchFlag)
		}

			
			
	}

	 
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
	 
		 
		value = toma->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = toma->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
	 			
		for(i=0; i<rhoN; i++) 
			Otrfrhos[i] = toma->m_Otrfrhos[i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = toma->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = toma->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = toma->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = toma->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  
	 
	delete toma;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/
	delete m_kCode;
	delete m_SrTemp;
	delete m_SrVTemp;

	return value;
		
		
} 
  





 double ODYSSEY_API  _stdcall TomCDAFn(int *uAssetN, double *sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,	                                   
	                                  int matuN, int *matu, double *prate, double ctime,					
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  long SimulNo, int batchFlag, char *kCode)
 

{
 	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	 

	int i;

		  
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}
	 
	int m_batchFlag;
	char *m_kCode;
	int outN, batchi;

	outN = 2*SoutRange + 1;
	 
	double *m_SrTemp;
	double *m_SrVTemp;
	 


	m_SrTemp = new double[outN+4];
	m_SrVTemp = new double[outN+4];
	 
	 
	
	m_batchFlag = batchFlag;
	m_kCode = new char[kCodeN+1];
	strcpy_s(m_kCode,kCodeN,kCode);

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	TomCDA *tomcda = new TomCDA(sval[0], stom, btom, STDayN, STStartDay, IXT, IndexT, STEndDay, IXTheta, 		  
	      matuN, matu, prate, ctime,
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		 divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  SimulNo);

	tomcda->m_modeltype = uAssetN[1];
			  
	 
	tomcda->m_modeltype = uAssetN[1];  
			  
	for(i=0; i<MAXGREEKS1; i++)
		tomcda->m_Ogreeks[i] = 0;
	
	if(m_batchFlag == 1) {
		for(batchi=0; batchi<outN+4; batchi++) {
			m_SrTemp[batchi] = sval[0]*(1-(SoutRange+2)*0.01 + 0.01*batchi);
			m_SrVTemp[batchi] = tomcda->getBatch(m_SrTemp[batchi]);
		}

		
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
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] << endl;
				}
				batchsf.close();
				
 				 
				delete fname;

	}
	
	/* value & 기본 그릭 산출 */    
	value = tomcda->getValue();
	Ogreeks[0] = tomcda->m_Ogreeks[0]; // price
	
	if(greektype == 0 || greektype == 4) {// 0: Total, 1: P,  2:V, 3: R  4:D&G & Theta,
		value = tomcda->getDeltaGammaTheta();
		for(i=1; i<BasicGRsN1; i++)
			Ogreeks[i] = tomcda->m_Ogreeks[i];

		
		if(m_batchFlag == 1) {	 
		
					using namespace std;
					using std::ofstream;

					ofstream batchsf;
			
					char *fname;
					fname = new char[kCodeN+1];

					strcpy_s(fname,kCodeN,m_kCode);
					strcat_s(fname,kCodeN,"_day+.txt");

					batchsf.open(fname,ios::out);
			
					batchsf << fixed;
					batchsf.precision(15);
					for(batchi=0; batchi<outN + 4; batchi++) {
						batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] + Ogreeks[Thetaidx] << endl;
					}
					batchsf.close();
					
 				 
				delete fname;
					 

		}

	}

	 


	/* Vega 관련 그릭 산출 */
 
	if(greektype == 0 || greektype == 2) {		 
		value = tomcda->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = tomcda->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		 
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = tomcda->m_Otvegas[i];

			
		
			if(m_batchFlag) {

				using namespace std;
				using std::ofstream;

				ofstream batchsf;
			
				char *fname;
				char *tempstring;
				fname = new char[kCodeN+1];
				tempstring = new char[kCodeN+1];


				 
				strcpy_s(fname,kCodeN,m_kCode);			 			
				 
					strcat_s(fname,kCodeN,"_sigma+_");
				 
				_itoa_s(i,tempstring,kCodeN,10);			
				strcat_s(fname,kCodeN,tempstring);			
				strcat_s(fname,kCodeN,".txt");


		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] + Otvegas[i]<< endl;
				}
				batchsf.close();
				
				 
				strcpy_s(fname,kCodeN,m_kCode);			 			
				 
					strcat_s(fname,kCodeN,"_sigma-_");
				 		
				_itoa_s(i,tempstring,kCodeN,10);			
				strcat_s(fname,kCodeN,tempstring);			
				strcat_s(fname,kCodeN,".txt");


		 
				batchsf.open(fname,ios::out);	  	 
			
				batchsf << fixed;
				batchsf.precision(15);
	 
				for(batchi=0; batchi<outN + 4; batchi++) {
					batchsf << m_SrTemp[batchi] << " " << m_SrVTemp[batchi] - Otvegas[i]<< endl;
				}
				batchsf.close();

 				delete tempstring;
				delete fname;
			
			}// if(m_batchFlag)


		}
			
			
	}

	 
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
	 
		 
		value = tomcda->getRho(rhoN, rhoidx);

		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = tomcda->m_Ogreeks[i]; 
		}	 
		 

		/* term별 rhos 할당 */
	 			
		for(i=0; i<rhoN; i++) 
			Otrfrhos[i] = tomcda->m_Otrfrhos[i];			 
		
		
		for(i=0; i<rhoN; i++) {
			Otdcrhos[i] = tomcda->m_Otdcrhos[i];			
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	/*
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/// CorrealtionDelta 관련 그릭 산출 
	if(greektype == 0 || greektype == 4) {	 
		value = tomcda->getCorrDelta(corrdeltaN, corrdeltaidx);


		///* Sum rhos 할당 
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = tomcda->m_Ogreeks[i]; 
		}

		///* term별 rhos 할당 
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = tomcda->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  */
  
	delete tomcda;
	 

	return value;
		
} 
  





//3S

 
 double ODYSSEY_API  _stdcall MakeHTAFn(char *kCode, int dataN, double baseval, int rangeN, double *Oresult)
 

{
	  
	 
	double value;	

	using namespace std;
	ifstream inFile;

	inFile.open(kCode);

	if(!inFile.is_open()) {
		value = 0;
	}
	else {
		value = 1;

		int i;
		double x[100], y[100];
		double tempx, tempy;

		for(i=0; i<dataN; i++) {
			inFile >> tempx;
			inFile >> tempy;
			x[i] = tempx;
			y[i] = tempy;
		}
		inFile.close();
 


		for(i=-rangeN; i<=rangeN; i++) {
			Oresult[rangeN + i] = Sinterp1(baseval+i,x,y,dataN,1,0);
		}			
		 



	}
	 
	 
	/* 
	double value;

	FILE *fpt;

	fpt = fopen(kCode,"r");

	if(fpt == NULL) {
		value = 0;
	}
	else {
		value = 1;
		int i;
		double *x, *y;
		double tempx, tempy;
		double tempinput;

		x = new double[dataN];
		y = new double[dataN];

		for(i=0; i< dataN; i++) {
			fscanf(fpt,"%f %f", &tempx, &tempy);
			x[i] = tempx;
			y[i] = tempy;
		}

		fclose(fpt);

		for(i=-rangeN; i<=rangeN; i++) {
			tempinput = baseval + i;
			Oresult[i] = Sinterp1(tempinput,x,y,dataN,0,0);
		}

		delete x;
		delete y;

	}
	*/  

	return value;
		
} 
  
  
  
 

//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
double Sinterp1(double x, double *Dx, double *Dy, int N, int meshType, int interpType)
{
	//N은 Dx의 Data #로 Dx[0] ~ Dx[N-1] 까지 사용됨
	double val;
	double ds;
	int idx;
	int maxidx, minidx, mididx;

	if(x<=Dx[0]) { // Xpoint 최소값 보다 작을 때
		if(interpType == 0 || N == 1) { //끝점 flat Data가 하나 일때 자동 flat
			val = Dy[0];
		}
		else { // 끝점 근방으로 선형
			val = (Dy[1]-Dy[0])/(Dx[1]-Dx[0]) *(x-Dx[0]) + Dy[0];
		}
	}
	else {
		if(x>=Dx[N-1]) { // Xpoint 최대값 보다 클 때
			if(interpType == 0 || N == 1) {
				val = Dy[N-1];
			}
			else {
				val = (Dy[N-1]-Dy[N-2])/(Dx[N-1]-Dx[N-2]) *(x-Dx[N-2]) + Dy[N-2];
			}	
		}	
		else { // Xpoint 사이에 존재
			if(meshType == 0) { // 균등일 때
				ds = Dx[1] - Dx[0];
				idx = (int)floor((x-Dx[0])/ds);
			}
			else { // 비균등일 때.. bisection으로 
				idx = -1;
				minidx = 0;
				maxidx = N-2;
				do {	
					mididx = (int)floor((maxidx+minidx)/2.0);
					if(Dx[mididx] <= x && x<Dx[mididx+1])
						idx = mididx;
					else {
						if(x<Dx[mididx])
							maxidx = mididx;
						else
							minidx = mididx+1;
					}
					 
				} while (idx == -1);									 
			}
			val = (Dy[idx+1]-Dy[idx])/(Dx[idx+1]-Dx[idx]) *(x-Dx[idx]) + Dy[idx];
		}
	}
 
	 
	return val;
}


// input x와 Dx가 int 일때..
//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
double Sinterp1(int x, int *Dx, double *Dy, int N, int meshType, int interpType)
{
	//N은 Dx의 Data #로 Dx[0] ~ Dx[N-1] 까지 사용됨
	double val;
	double ds;
	int idx;
	int maxidx, minidx, mididx;

	if(x<=Dx[0]) { // Xpoint 최소값 보다 작을 때
		if(interpType == 0 || N == 1) { //끝점 flat Data가 하마 일때 자동 flat
			val = Dy[0];
		}
		else { // 끝점 근방으로 선형
			val = (Dy[1]-Dy[0])/(Dx[1]-Dx[0]) *(x-Dx[0]) + Dy[0];
		}
	}
	else {
		if(x>=Dx[N-1]) { // Xpoint 최대값 보다 클 때
			if(interpType == 0 || N == 1) {
				val = Dy[N-1];
			}
			else {
				val = (Dy[N-1]-Dy[N-2])/(Dx[N-1]-Dx[N-2]) *(x-Dx[N-2]) + Dy[N-2];
			}	
		}	
		else { // Xpoint 사이에 존재
			if(meshType == 0) { // 균등일 때
				ds = Dx[1] - Dx[0];
				idx = (int)floor((x-Dx[0])/ds);
			}
			else { // 비균등일 때.. bisection으로 
				idx = -1;
				minidx = 0;
				maxidx = N-2;
				do {	
					mididx = (int)floor((maxidx+minidx)/2.0);
					if(Dx[mididx] <= x && x<Dx[mididx+1])
						idx = mididx;
					else {
						if(x<Dx[mididx])
							maxidx = mididx;
						else
							minidx = mididx+1;
					}
					 
				} while (idx == -1);									 
			}
			val = (Dy[idx+1]-Dy[idx])/(Dx[idx+1]-Dx[idx]) *(x-Dx[idx]) + Dy[idx];
		}
	}
 
	 
	return val;
}




 double ODYSSEY_API  _stdcall VBStepDownAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag,  char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double *inkihval;
	double *psval;
	int i;
	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	inkihval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) 	{	
			inxval[i] = xval[i]*bval[0];		
			inkihval[i] = hval[i]*bval[0];
		}
		 

	}
	else {
		for(i=0; i<matuN; i++)		{
			inxval[i] = xval[i];		
			inkihval[i] = hval[i];
		}		 
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];
 
	VBStepDownA *vbstepdowna = new VBStepDownA(knockFlag[0], sval[0], bval[0], inxval, inkihval, 
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	vbstepdowna->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = vbstepdowna->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = vbstepdowna->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = vbstepdowna->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = vbstepdowna->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = vbstepdowna->m_Otvegas[i];
			Otvannas[i] = vbstepdowna->m_Otvannas[i];
			Otzommas[i] = vbstepdowna->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = vbstepdowna->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = vbstepdowna->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = vbstepdowna->m_Otrfrhos[i];
			Otdcrhos[i] = vbstepdowna->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete inkihval;
	delete psval;
	delete vbstepdowna;

	return value;
		
} 
  






 double ODYSSEY_API  _stdcall VBStepDownCDAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 0;
			matu[matuN-1] = 50000;
			SmeshM[0] = 50000;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}

	double value;	
	double *inxval;
	double *inkihval;
	double *psval;
	int i;

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[matuN];
	inkihval = new double[matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i]*bval[0];		
			inkihval[i] = hval[i]*bval[0];
		}

	 

	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inkihval[i] = hval[i];
		}
		 
	}

	psval = new double[maxevalN];
	for(i=0; i<maxevalN; i++)
		psval[i] = sval[i];

 //int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
	VBStepDownCDA *vbstepdowncda = new VBStepDownCDA(knockFlag[0], sval[0], bval[0], inxval, inkihval, 
		  Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,
		  CDN, CDadjval, ObBDay, PayBDay, FlegFactor, ObCDrate,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate[0], divbval[0], divN[0], divBDay, divCash, divApy,
		  voltype, volbval[0], volTN, volSN, volBDay, volSmness, vol,		  		
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	vbstepdowncda->m_modeltype = uAssetN[1];
			  
	/* value & 기본 그릭 산출 */    
	value = vbstepdowncda->getValue();	 
	for(i=0; i<BasicGRsN1; i++) {	
		Ogreeks[i] = vbstepdowncda->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	if(greektype == 0 || greektype == 2) {		 
		value = vbstepdowncda->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN1; i<BasicGRsN1+VegaGRsN1; i++) {  // 7, 8, 9
			Ogreeks[i] = vbstepdowncda->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas[i] = vbstepdowncda->m_Otvegas[i];
			Otvannas[i] = vbstepdowncda->m_Otvannas[i];
			Otzommas[i] = vbstepdowncda->m_Otzommas[i];
		}
	}


	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = vbstepdowncda->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
			Ogreeks[i] = vbstepdowncda->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = vbstepdowncda->m_Otrfrhos[i];
			Otdcrhos[i] = vbstepdowncda->m_Otdcrhos[i];
		}
		
///////////////////////////////////////////--------------Debug	Mode s
		int debugFlag;
		debugFlag = 0;
		if(debugFlag) {

			using namespace std;
			using std::ofstream;

			ofstream logsf; 
			logsf.open("c:/logland/debug_gtype3.txt",ios::out);
 
			logsf << "value = "<< value << endl;

			for(i=BasicGRsN1+VegaGRsN1; i<BasicGRsN1+VegaGRsN1+RhoGRsN1; i++) { // 10, 11
				logsf << "Ogreeks[" << i <<"]" <<"= "<< Ogreeks[i] << endl;			 
			}		
			
			for(i=0; i<rhoN; i++) {
				logsf << "Otrfrhos[" << i <<"]" <<"= "<< Otrfrhos[i] << endl;			 
			}		
					
			for(i=0; i<rhoN; i++) {
				logsf << "Otdcrhos[" << i <<"]" <<"= "<< Otdcrhos[i] << endl;			 
			}		
	 
	 
			logsf.close(); 
		
		}
///////////////////////////////////////////--------------Debug	Mode e	

	}//if(greektype == 0 || greektype == 3)

 

	delete inxval;
	delete inkihval;
	delete psval;
	delete vbstepdowncda;

	return value;
		
} 
  






 double ODYSSEY_API  _stdcall VBStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	 
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	 
	 
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN];
	inkihval = new double[uAssetN[0]*matuN];
	if(xvaltype == 100) {
		for(i=0; i<matuN; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN+i] = xval[matuN+i]*bval[1];

			inkihval[i] = hval[i]*bval[0];		
			inkihval[matuN+i] = hval[matuN+i]*bval[1];
		}
 
	}
	else {
		for(i=0; i<matuN; i++) {		
			inxval[i] = xval[i];
			inxval[matuN+i] = xval[matuN+i];

			inkihval[i] = hval[i];		
			inkihval[matuN+i] = hval[matuN+i];
		}
		 
	}
 
 
	VBStepDownB *vbstepdownb = new VBStepDownB(knockFlag[0],  insval,  bval, inxval, inkihval,   Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,			  
		  irateN, irateBDay, irateCDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	vbstepdownb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
	                                     // 2: stepdown ko (원금비보장 ko)
			  
	/* value & 기본 그릭 산출 */    
	value = vbstepdownb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = vbstepdownb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = vbstepdownb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = vbstepdownb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = vbstepdownb->m_Otvegas1[i];
			Otvannas1[i] = vbstepdownb->m_Otvannas1[i];
			Otzommas1[i] = vbstepdownb->m_Otzommas1[i];
						
			Otvegas2[i] = vbstepdownb->m_Otvegas2[i];
			Otvannas2[i] = vbstepdownb->m_Otvannas2[i];
			Otzommas2[i] = vbstepdownb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = vbstepdownb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = vbstepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = vbstepdownb->m_Otrfrhos[i];
			Otdcrhos[i] = vbstepdownb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = vbstepdownb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = vbstepdownb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = vbstepdownb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	 
	delete psval;

	delete vbstepdownb;

	return value;
		
} 
  



 


 double ODYSSEY_API  _stdcall VBStepDownCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
	if(matu[matuN[0]-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
		if(veriipq == 0) { // 부적격
			greektype = 4;
			matu[matuN[0]-1] = 50000;
			SmeshM[0] = 200;
			TmeshN[1] = 1000;
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	double *insval;
	double *inxval;
	double *inkihval;
	 
	double *psval;

	int i;


	 
	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	insval = new double[2];
	 
	 
		
	insval[0] = sval[0];
	insval[1] = sval[maxevalN];
	
	psval = new double[uAssetN[0]*maxevalN];
	for(i=0; i<uAssetN[0]*maxevalN; i++)
		psval[i] = sval[i];

	 

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	inxval = new double[uAssetN[0]*matuN[0]];
	inkihval = new double[uAssetN[0]*matuN[0]];

	if(xvaltype == 100) {
		for(i=0; i<matuN[0]; i++) { 		
			inxval[i] = xval[i]*bval[0];
			inxval[matuN[0]+i] = xval[matuN[0]+i]*bval[1];
					
			inkihval[i] = hval[i]*bval[0];
			inkihval[matuN[0]+i] = hval[matuN[0]+i]*bval[1];		 
		}
		 
	}
	else {
		for(i=0; i<matuN[0]; i++) {		
			inxval[i] = xval[i];
			inxval[matuN[0]+i] = xval[matuN[0]+i];
			
			inkihval[i] = hval[i];
			inkihval[matuN[0]+i] = hval[matuN[0]+i];
		}
	 
	}
 
 
	VBStepDownCDB *vbstepdowncdb = new VBStepDownCDB(knockFlag[0], insval,  bval, inxval, inkihval, Cpn,
	      matuN, matu, ctime,
		  maxevalN, evalN, psval,	
		  CDlegInfo,
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  corrtype, corrN, corrBDay, corr,
		  batchFlag, kCode,
		  Sminlevel, Smaxlevel, SmeshM, TmeshN,
		  shiftHRate, spreadHRate, spreadXRate);

	 
	vbstepdowncdb->m_modeltype = uAssetN[1]; // 0: stepdown (원금비보장)
	                                    // 1  stepdown (원금보장 KI 있는)
			  
	/* value & 기본 그릭 산출 */    
	value = vbstepdowncdb->getValue();	 
	for(i=0; i<BasicGRsN2; i++) {	
		Ogreeks[i] = vbstepdowncdb->m_Ogreeks[i];
	}


	/* Vega 관련 그릭 산출 */
	int vegaN;
	vegaN = termgrNs[0];
	if(greektype == 0 || greektype == 2) {		 
		value = vbstepdowncdb->getVega(vegaN, vegaidx);
		
		/* Sum vegas 할당 */
		for(i=BasicGRsN2; i<BasicGRsN2+VegaGRsN2; i++) {  // 7, 8, 9
			Ogreeks[i] = vbstepdowncdb->m_Ogreeks[i];
		}

		/* term별 vegas 할당 */
		for(i=0; i<vegaN; i++) {
			Otvegas1[i] = vbstepdowncdb->m_Otvegas1[i];
			Otvannas1[i] = vbstepdowncdb->m_Otvannas1[i];
			Otzommas1[i] = vbstepdowncdb->m_Otzommas1[i];
						
			Otvegas2[i] = vbstepdowncdb->m_Otvegas2[i];
			Otvannas2[i] = vbstepdowncdb->m_Otvannas2[i];
			Otzommas2[i] = vbstepdowncdb->m_Otzommas2[i];
		}
	}

	int rhoN;
	rhoN = termgrNs[1];
	/* Rho 관련 그릭 산출 */
	if(greektype == 0 || greektype == 3) {	 
		value = vbstepdowncdb->getRho(rhoN, rhoidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2; i++) { // 10, 11
			Ogreeks[i] = vbstepdowncdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<rhoN; i++) {
			Otrfrhos[i] = vbstepdowncdb->m_Otrfrhos[i];
			Otdcrhos[i] = vbstepdowncdb->m_Otdcrhos[i];
		}
 
	}//if(greektype == 0 || greektype == 3)
	 
	int corrdeltaN;
	corrdeltaN = termgrNs[2];
	/* CorrealtionDelta 관련 그릭 산출 */
	if(greektype == 0 || greektype == 4) {	 
		value = vbstepdowncdb->getCorrDelta(corrdeltaN, corrdeltaidx);


		/* Sum rhos 할당 */
		for(i=BasicGRsN2+VegaGRsN2+RhoGRsN2; i<BasicGRsN2+VegaGRsN2+RhoGRsN2+CorrDeltaGRsN2; i++) { // 10, 11
			Ogreeks[i] = vbstepdowncdb->m_Ogreeks[i]; 
		}

		/* term별 rhos 할당 */
		for(i=0; i<corrdeltaN; i++) {
			Otcorrdelta[i] = vbstepdowncdb->m_Otcorrdelta[i];			 
		}
 
	}//if(greektype == 0 || greektype == 3)
	  
  

	delete insval;
	delete inxval;
	delete inkihval;
	delete psval;
 

	delete vbstepdowncdb;

	return value;
		
} 
  




 
 double ODYSSEY_API  _stdcall mcPlainFn(int cpFlag, double *sval, double *xval, 	                                  
	                                  int matuN, int *matu, 							   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  long SimulNo)
									  /* 주의 
									   SmaxlevelC 배열이 아니라 Smaxlevel 임
									   S축 Adaptive mesh는 만기기준으로만 적용
									   shiftXRate 는 없고, spreadXRate가 배열임
									   */
 

{
 	if(matu[matuN-1] < 0) 
		return -999999;
		
	 
		
	//ip check
	veriipq = veriipqFlag;
	if( veriipq == 0) {
		/*
		if(::WSAStartup(MAKEWORD(1, 0), &WSAData))
		{
		  // Error handling
		}
		if(::gethostname(szHostName, sizeof(szHostName)))
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		pHost = ::gethostbyname(szHostName);
		if(!pHost)
		{
		  // Error handling -> call 'WSAGetLastError()'
		}
		*/
		while(::WSAStartup(MAKEWORD(1, 0), &WSAData)) {}

		while(::gethostname(szHostName, sizeof(szHostName))) {}

		while(!(pHost = ::gethostbyname(szHostName))) {}

		veriipq = 0;

		for(ipCnt = 0; ((pHost->h_addr_list[ipCnt]) && (ipCnt < 10)); ++ipCnt)
		{
		  strcpy_s(aszIPAddresses[ipCnt],16, inet_ntoa(*(struct in_addr*)pHost->h_addr_list[ipCnt]));
		  ipfindi = 0;
		  for(ipi=0; ipi<16; ipi++) {
			  if(aszIPAddresses[ipCnt][ipi] == '.') {
				  ipfindi++;
				  if(ipfindi == 2) {
					  iptargeti = ipi;
					  break;
				  }
			  }
		  }	 
		  strcpy_s(temps,16,aszIPAddresses[ipCnt]); 
		  temps[iptargeti+1]='\0';	   
	   
		  for(ipi=0; ipi<10; ipi++) {
			  if(strcmp(temps,veriip[ipi]) == 0)
				  veriipq = 1;
		  } 
		}
		WSACleanup();		 
	
	
		if(veriipq == 0) { // 부적격			
			matu[matuN-1] = 50000;			 
			sval[0] = sval[0]*0.5;
			vol[0] = 0.00001;
		}
	}
	
	double value;	
	 
 
		  
	int ownrand;


	ownrand = 0;

	if(RandUp == 0) {	
		//odysseyRandom_INIT(SimulNo);
		ownrand = 1;
	}

	mcPlain *mcplain = new mcPlain(cpFlag, sval[0], xval[0], 		  
	      matuN, matu, 	   			  
		  irateN, irateBDay, iRFrate, iDCrate, iIRSrate,
		  divrate, divbval, divN, divBDay, divCash, divApy,
		  voltype, volbval, volTN, volSN, volBDay, volSmness, vol,		
		  SimulNo);

	 
	/* value & 기본 그릭 산출 */    
	value = mcplain->getValue();
		 



	 
	delete mcplain;

	/*
	if(ownrand == 1)	
		odysseyRandom_END(SimulNo);
	*/

	return value;
		
		
} 
  




 int ODYSSEY_API  _stdcall SVolTableCal(double sval,  		   
									  int irateN, int *irateCDay, double *iRFrate, 
									  int divrateN, int *divrateCDay, double *divrate, 
									  int *divN, int *divBDay, double *divCash, double divApy,
									  int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int *mkcpFlag, double *mkdata, double *mkvega, double *outvol, double *outdata,									  
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN, int maxloop, double *inmytol)
 
{
	 

// 변수 간소화
	  
	 
 


//	double value;	
	 
	int i,j,k,ai,aj,ak;
	int tj, si;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	 
 


	//Plainq **plainq;	 
	int loopcnt;
	double pa = 0.0001;
	double vegazeromytol;
	double mytol;
	double lambda;

	double *cy;
	double *newcy;
	double oldss, newss;
	double **pcya;
	double *beta;
	double **alpha, **alphad;
	double *da;
	double *tempvol;
	double *tempa;
	double maxerror;
	int rankA;
	int maxraw;
	double tempsum;
	double tempvalue;
	double *value2;

	vegazeromytol = inmytol[0];
	mytol = inmytol[1];

	cy = new double[volSN];
	newcy = new double[volSN];
	tempvol = new double[volSN*volTN];
	tempa = new double[volSN];

	beta = new double[volSN];
	da = new double[volSN];

	pcya = new double*[volSN];
	alpha = new double*[volSN];
	alphad = new double*[volSN];
	for(i=0; i<volSN; i++) {
		pcya[i] = new double[volSN];
		alpha[i] = new double[volSN];
		alphad[i] = new double[volSN];
	}

	
		
	using namespace std;	
	using std::ofstream;
	ofstream logs;
 
  
	lambda = 0.001;

	 
 

	 
		//plainq = new Plainq*[calnum];
	 
	value2 = new double[volSN*volTN];
 
		// pricing s
		int matu;
		int cpFlag;
		double xval;
		int voltype;		 
		int impvolTN=1;
		int impvolSN=1;
		double *impvol;
		int *impvolBDay;
		double *impvolSmness;

	 
		impvol = new double[1];
		impvolBDay = new int[1];
		impvolSmness = new double[1];
 
		impvolBDay[0] = 365;
		impvolSmness[0] = 1.0;

		voltype = 0; // impvol로 각 grid mk price 계산!
		for(i=0; i<volSN; i++) {
			xval = volSmness[i]*sval;
			for(j=0; j<volTN; j++) {
				matu = volBDay[j];
				cpFlag = mkcpFlag[i*volTN+j];
				impvol[0] = vol[i*volTN+j];				 
				 
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
		                            irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);
				mkdata[i*volTN+j] = plainq->getValue();
				delete plainq;
  
			}
		}

	  
		 
		//memcpy(mkdata,value2,(size_t)((volSN*volTN)*sizeof(double)));
/////////////////////////////////////////////////////////////////////////////////// mk price 계산!
 
  

////////////// for vega s
 
		 
		for(i=0; i<volSN; i++) {
			xval = volSmness[i]*sval;
			for(j=0; j<volTN; j++) {
				matu = volBDay[j];
				cpFlag = mkcpFlag[i*volTN+j];
				impvol[0] = vol[i*volTN+j] + 0.01;
				 
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
		                            irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);

				mkvega[i*volTN+j] = plainq->getValue();
				delete plainq;
  
			}
		}

	 

		 
		//memcpy(mkvega,value2,(size_t)((volSN*volTN)*sizeof(double)));
///// -
//		 
//		 
//		for(i=0; i<volSN; i++) {
//			xval = volSmness[i]*sval;
//			for(j=0; j<volTN; j++) {
//				matu = volBDay[j];
//				cpFlag = mkcpFlag[i*volTN+j];
//				impvol[0] = vol[i*volTN+j] - 0.01;
//				impvolBDay[0] = matu;
//				impvolSmness[0] = volSmness[i]; // =1.0 이나 결과는 같은?
//				  
//				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
//		                            irateN, irateCDay, iRFrate,  
//									divrateN, divrateCDay, divrate, 
//									sval, divN[0], divBDay, divCash, divApy,
//		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
//		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);
//				value2[i*volTN+j] = plainq->getValue();
//				delete plainq;  
//			}
//		}
//
//	  	
// 
//
//		for(i=0; i<volSN*volTN; i++)
//			value2[i] = (mkvega[i] - value2[i])/2.0;
//		

		for(i=0; i<volSN*volTN; i++)
			mkvega[i] = (mkvega[i] - mkdata[i]);

		//memcpy(mkvega,value2,(size_t)((volSN*volTN)*sizeof(double)));
		delete value2;
		delete impvol;
		delete impvolBDay;
		delete impvolSmness;
////////////////////// for vega end

			
 



		memcpy(outvol,vol,(size_t)((volSN*volTN)*sizeof(double)));

		for(i=0; i<volSN; i++) {
			for(j=0; j<volTN; j++) {
				outvol[i*volTN+j] = vol[12*volTN+j];
			}
		}
 
///////// start here!! #1

		double *tvalue; // 각 tj 마다 찾을 mk option pr
		int *pidx;
		int targetN;

		tvalue = new double[volSN]; // nonzero vega mk price 저장
		value2 = new double[volSN];
		pidx = new int[volSN]; // nonzero vega index 저장
		voltype = 2; // outvol로 각 grid cy price 계산!
		for(tj=0; tj<volTN; tj++) {
			
 
			//vega 0 = mytol/10 수준으로 설정 찾기
			targetN = 0;
			for(si=0; si<volSN; si++) {
				if(mkvega[si*volTN + tj] >= vegazeromytol) {
					//tvalue[targetN] = mkdata[si*volTN + tj]; 아래의 표현으로 test
					pidx[targetN] = si;
					targetN++; // vega 0 가 아닌 갯수
				} // if
			} // for si

			if(targetN == 0) return maxloop;

			// outvol 준비
			for(si=0; si<targetN; si++) {
				tvalue[si] = mkdata[pidx[si]*volTN + tj];
				outvol[pidx[si]*volTN + tj] = vol[pidx[si]*volTN + tj]; // imp vol 수준의 initial guess vol start!!
			}
			for(i=0; i<volSN; i++) {
				if(i<pidx[0]) {
					outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
				}
				if(i>pidx[targetN-1]) {
					outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기
				}
			}
			 			
						
if(1) {

			using namespace std;
			using std::ofstream;

			ofstream batchsf;
			
			char *fname;
			char *tempstring;
			fname = new char[kCodeN+1];
			tempstring = new char[kCodeN+1];
			//strcpy_s(m_kCode,kCodeN,kCode); 
			strcpy_s(fname,kCodeN,"c:/logland/tj_");
			_itoa_s(tj,tempstring,kCodeN,10);			
			strcat_s(fname,kCodeN,tempstring);			
			strcat_s(fname,kCodeN,"s.txt");


		 
			batchsf.open(fname,ios::out);	  	 
			
			batchsf << fixed;
			batchsf.precision(15);
	  
	
			for(i=0; i<volSN; i++)
				batchsf << "vol[" << i << "]=  "  <<  vol[i*volTN+tj]  << endl;

			for(i=0; i<volSN; i++)
				batchsf << "mkdata[" << i << "]=  "  << mkdata[i*volTN+tj] << endl;
			
			for(i=0; i<volSN; i++)
				batchsf << "mkvega[" << i << "]=  " << mkvega[i*volTN+tj] << endl;
			
			for(i=0; i<volSN; i++)
				batchsf << "outvol[" << i << "]=  "  << outvol[i*volTN+tj] << endl;
			
			for(i=0; i<targetN; i++)
				batchsf << "pidx[" << i << "]=  " << pidx[i] << endl;
			
			for(i=0; i<targetN; i++)
				batchsf << "tvalue[" << i << "]=  " << tvalue[i] << endl;
			 
			batchsf.close();

 			delete tempstring;
			delete fname;
}
			// guess vol 로 일차 cy 계산하기!!
 
			 
			matu = volBDay[tj];
			for(si=0; si<targetN; si++) {
				xval = volSmness[pidx[si]]*sval;
				cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
									irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
									voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
									Sminlevel, SmaxlevelC, SmeshM, TmeshN);
				cy[si] = plainq->getValue();
				delete plainq;  
   
			}

	  	   
			// LM 시작
			oldss = 0.0;
			for(si=0; si<targetN; si++) 
				oldss = oldss + (tvalue[si]-cy[si])*(tvalue[si]-cy[si]);
					
		
			//pcya 계산
			memcpy(tempvol,outvol,(size_t)((volSN*volTN)*sizeof(double)));					
			for(ai=0; ai<targetN; ai++) {

				memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));			
				outvol[pidx[ai]*volTN + tj] = tempvol[pidx[ai]*volTN + tj] + pa;						
				for(i=0; i<volSN; i++) {
					if(i<pidx[0]) {
						outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
					}
					if(i>pidx[targetN-1]) {
						outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기
					}
				}				
			
				//pcya 계산 하기위함
				  
				 
				for(si=0; si<targetN; si++) {
					xval = volSmness[pidx[si]]*sval;
					cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
					Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
										irateN, irateCDay, iRFrate,  
										divrateN, divrateCDay, divrate, 
										sval, divN[0], divBDay, divCash, divApy,
										voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
										Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   
					value2[si] = plainq->getValue();				
					delete plainq;  
				}

  
				for(si=0; si<targetN; si++)						
					pcya[ai][si] = (value2[si]-cy[si])/pa;				 
						 
			}
			memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
			//// pcya 계산 end!!

			

			for(ai=0; ai<targetN; ai++) {
				beta[ai] = 0.0;
				for(k=0; k<targetN; k++)
					beta[ai] = beta[ai] + (tvalue[k]-cy[k])*pcya[ai][k];
			}

	 

			for(ai=0; ai<targetN; ai++) {
				for(aj=0; aj<targetN; aj++) {
					alpha[ai][aj] = 0.0;
					for(k=0; k<targetN; k++) 
						alpha[ai][aj] = alpha[ai][aj] + pcya[ai][k]*pcya[aj][k];
				}
			}
		
	 

			loopcnt = 0;
			maxerror = 0;
			for(ai=0; ai<targetN; ai++) 
				if(maxerror < abs(tvalue[ai]-cy[ai]))
					maxerror = abs(tvalue[ai]-cy[ai]);
	 
		
			while (loopcnt < maxloop && maxerror > mytol ) {

				for(ai=0; ai<targetN; ai++) {
					for(aj=0; aj<targetN; aj++) {
						if(ai == aj)
							alphad[ai][aj] = alpha[ai][aj]*(1+lambda);
						else
							alphad[ai][aj] = alpha[ai][aj];
					}
				}

		 
				//da = inv(alphad)*beta;
				// alphad * da = beta 풀기 s
				memcpy(tempa, beta,(size_t)((volSN)*sizeof(double)));
				rankA = targetN;
				for(ai=0; ai<targetN; ai++)  {
					maxraw = ai;
					for(aj=ai; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++)
							tempsum = tempsum + alphad[aj][ak]*alphad[ak][ai];
						if((alphad[aj][ai]-tempsum) != 0) {
							maxraw = aj;
							break;
						}
					}

					if(maxraw != ai) {
						for(ak=0; ak<targetN; ak++) {
							tempvalue = alphad[ai][ak];
							alphad[ai][ak] = alphad[maxraw][ak];
							alphad[maxraw][ak] = tempvalue;
						}
						tempvalue = beta[ai];
						beta[ai] = beta[maxraw];
						beta[maxraw] = tempvalue;
					}

					for(aj=ai; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[ai][ak]*alphad[ak][aj];
						alphad[ai][aj] = alphad[ai][aj]-tempsum;
					}

					if(alphad[ai][ai] == 0) {
						rankA = ai;				
						break;
					}

					for(aj=ai+1; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[aj][ak]*alphad[ak][ai];
						alphad[aj][ai] = (alphad[aj][ai]-tempsum)/alphad[ai][ai];
					}
				} // for(ai.. ) // end lu

				if(rankA == targetN) {
					for(ai=0; ai<targetN; ai++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[ai][ak]*da[ak];
						da[ai] = beta[ai] - tempsum;
					}

					for(ai=targetN-1; ai>=0; ai--) {
						tempsum = 0;
						for(ak=targetN-1; ak>ai; ak--)
							tempsum = tempsum + alphad[ai][ak]*da[ak];
						da[ai] = (da[ai]-tempsum)/alphad[ai][ai];
					}			
				}
				else {
					loopcnt = 2*maxloop; //
				}
				memcpy(beta, tempa,(size_t)((volSN)*sizeof(double)));
				// alphad * da = beta 풀기 e
		 

				memcpy(tempvol, outvol,(size_t)((volSN*volTN)*sizeof(double)));
				for(ai=0; ai<targetN; ai++) 
					outvol[pidx[ai]*volTN + tj] = outvol[pidx[ai]*volTN + tj] + da[ai];
							
				for(i=0; i<volSN; i++) {
					if(i<pidx[0]) {
						outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
					}
					if(i>pidx[targetN-1]) {
						outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기
					}
				}				
		 
			 //vector<PlainqRunnable*> runnable;  // multi
		 
			   
				for(si=0; si<targetN; si++) {
					xval = volSmness[pidx[si]]*sval;
					cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
					Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
										irateN, irateCDay, iRFrate,  
										divrateN, divrateCDay, divrate, 
										sval, divN[0], divBDay, divCash, divApy,
										voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
										Sminlevel, SmaxlevelC, SmeshM, TmeshN);
					
					newcy[si] = plainq->getValue();				
					delete plainq;  
   
				}   
	 	
				memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
 			
				newss = 0;
				for(si=0; si<volSN; si++)					
					newss = newss + (tvalue[si] - newcy[si])*(tvalue[si] - newcy[si]);


				if(newss >= oldss)
					lambda = 10*lambda;
				else {
					lambda = 0.1*lambda;

					for(ai=0; ai<targetN; ai++) 				
						outvol[pidx[ai]*volTN + tj] = outvol[pidx[ai]*volTN + tj] + da[ai];
									
					for(i=0; i<volSN; i++) {
						if(i<pidx[0]) {
							outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
						}
						if(i>pidx[targetN-1]) {
							outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기
						}
					}		
					/// here !!
					memcpy(cy,newcy,(size_t)((volSN)*sizeof(double)));									
					maxerror = 0;		
					for(ai=0; ai<targetN; ai++) 			
						if(maxerror < abs(tvalue[ai]-cy[ai]))				
							maxerror = abs(tvalue[ai]-cy[ai]);
 		
					oldss = newss;

					memcpy(tempvol,outvol,(size_t)((volSN*volTN)*sizeof(double)));			
					for(ai=0; ai<targetN; ai++) {

						memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));			
						outvol[pidx[ai]*volTN + tj] = tempvol[pidx[ai]*volTN + tj] + pa;						
						for(i=0; i<volSN; i++) {
							if(i<pidx[0]) {
								outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
							}
							if(i>pidx[targetN-1]) {
								outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기
							}
						}				
			
					 
						for(si=0; si<targetN; si++) {
							xval = volSmness[pidx[si]]*sval;
							cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
							Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
												irateN, irateCDay, iRFrate,  
												divrateN, divrateCDay, divrate, 
												sval, divN[0], divBDay, divCash, divApy,
												voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
												Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   					
							value2[si] = plainq->getValue();									
							delete plainq;  
						}
 
  
  
						for(si=0; si<targetN; si++)						
							pcya[ai][si] = (value2[si]-cy[si])/pa;						
						 
					}
					memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
				//// pcya 계산 end!!
 

					for(ai=0; ai<targetN; ai++) {
						beta[ai] = 0.0;
						for(k=0; k<targetN; k++)
							beta[ai] = beta[ai] + (tvalue[k]-cy[k])*pcya[ai][k];
					}

					for(ai=0; ai<targetN; ai++) {
						for(aj=0; aj<targetN; aj++) {
							alpha[ai][aj] = 0.0;
							for(k=0; k<targetN; k++) 
								alpha[ai][aj] = alpha[ai][aj] + pcya[ai][k]*pcya[aj][k];
						}
					}


				} // if(newss.. ) .. else
				 
				loopcnt = loopcnt + 1;

			} // while
		
		//memcpy(outvol, vol,(size_t)((volSN*volTN)*sizeof(double)));
		//memcpy(outdata, cy,(size_t)((volSN)*sizeof(double)));
	 
			for(si=0; si<volSN; si++) {
				xval = volSmness[si]*sval;
				cpFlag = mkcpFlag[si*volTN + tj];
				  
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
									irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
									voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
									Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   
							
				outdata[si*volTN + tj] = plainq->getValue();													
				delete plainq;  
			}
 if(1) {

			using namespace std;
			using std::ofstream;

			ofstream batchsf;
			
			char *fname;
			char *tempstring;
			fname = new char[kCodeN+1];
			tempstring = new char[kCodeN+1];
			 
			strcpy_s(fname,kCodeN,"c:/logland/tj_");
			_itoa_s(tj,tempstring,kCodeN,10);			
			strcat_s(fname,kCodeN,tempstring);			
			strcat_s(fname,kCodeN,".txt");


		 
			batchsf.open(fname,ios::out);	  	 
			
			batchsf << fixed;
			batchsf.precision(15);
	  
		
			for(i=0; i<volSN; i++)
				batchsf << "vol[" << i << "]=  "  <<  vol[i*volTN+tj]  << endl;

			for(i=0; i<volSN; i++)
				batchsf << "mkdata[" << i << "]=  "  << mkdata[i*volTN+tj] << endl;
			
			for(i=0; i<volSN; i++)
				batchsf << "mkvega[" << i << "]=  " << mkvega[i*volTN+tj] << endl;
			
			for(i=0; i<volSN; i++)
				batchsf << "outvol[" << i << "]=  "  << outvol[i*volTN+tj] << endl;
			
			for(i=0; i<volSN; i++)
				batchsf << "outdata[" << i << "]=  " << outdata[i*volTN+tj] << endl;


			batchsf.close();

 			delete tempstring;
			delete fname;

}

		} // for tj..



		delete tvalue;
		delete value2;
		delete pidx;		

///////////////////////////////////////////////////////////////////// end here!! #1

  
	

	delete cy;
	delete newcy;

	delete beta;
	delete da;
 
	for(i=0; i<volSN; i++) {
		delete pcya[i];
		delete alpha[i];
		delete alphad[i];
	}	
	delete pcya;
	delete alpha;
	delete alphad;
 
 
	return loopcnt;
		
} 
  


 int ODYSSEY_API  _stdcall SVolTableCalCF(double sval,  		   
									  int irateN, int *irateCDay, double *iRFrate, 
									  int divrateN, int *divrateCDay, double *divrate, 
									  int *divN, int *divBDay, double *divCash, double divApy,
									  int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int *mkcpFlag, double *mkdata, double *mkvega, double *outvol, double *outdata,									  
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN, int maxloop, double *inmytol)
 
{
	 

// 변수 간소화
	  
	 
 


//	double value;	
	 
	int i,j,k,ai,aj,ak;
	int tj, si;

	//쿼트 할때는.. Sv난 Sd의 값과 비슷한 sval, bval등을 입력 또는
	//sval = 100, bval = 100, Sv = 100, Sd = 100, D는 D/Sd *100에 해당하는 값 입력으로..!!
	 
 


	//Plainq **plainq;	 
	int loopcnt;
	double pa = 0.0001;
	double vegazeromytol;
	double mytol;
	double lambda;

	double *cy;
	double *newcy;
	double oldss, newss;
	double **pcya;
	double *beta;
	double **alpha, **alphad;
	double *da;
	double *tempvol;
	double *tempa;
	double maxerror;
	int rankA;
	int maxraw;
	double tempsum;
	double tempvalue;
	double *value2;

	vegazeromytol = inmytol[0];
	mytol = inmytol[1];

	cy = new double[volSN];
	newcy = new double[volSN];
	tempvol = new double[volSN*volTN];
	tempa = new double[volSN];

	beta = new double[volSN];
	da = new double[volSN];

	pcya = new double*[volSN];
	alpha = new double*[volSN];
	alphad = new double*[volSN];
	for(i=0; i<volSN; i++) {
		pcya[i] = new double[volSN];
		alpha[i] = new double[volSN];
		alphad[i] = new double[volSN];
	}

	
		
	using namespace std;	
	using std::ofstream;
	ofstream logs;
 
  
	lambda = 0.001;

	 
 

	 
		//plainq = new Plainq*[calnum];
	 
	value2 = new double[volSN*volTN];
 
		// pricing s
		int matu;
		int cpFlag;
		double xval;
		int voltype;		 
		int impvolTN=1;
		int impvolSN=1;
		double *impvol;
		int *impvolBDay;
		double *impvolSmness;

	 
		impvol = new double[1];
		impvolBDay = new int[1];
		impvolSmness = new double[1];
 
		impvolBDay[0] = 365;
		impvolSmness[0] = 1.0;

		voltype = 0; // impvol로 각 grid mk price 계산!
		for(i=0; i<volSN; i++) {
			xval = volSmness[i]*sval;
			for(j=0; j<volTN; j++) {
				matu = volBDay[j];
				cpFlag = mkcpFlag[i*volTN+j];
				impvol[0] = vol[i*volTN+j];				 
				 
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
		                            irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);
				mkdata[i*volTN+j] = plainq->getcfValue();
				delete plainq;
  
			}
		}

				 			
 
		 
		//memcpy(mkdata,value2,(size_t)((volSN*volTN)*sizeof(double)));
/////////////////////////////////////////////////////////////////////////////////// mk price 계산!
 
  

////////////// for vega s
 
		 
		for(i=0; i<volSN; i++) {
			xval = volSmness[i]*sval;
			for(j=0; j<volTN; j++) {
				matu = volBDay[j];
				cpFlag = mkcpFlag[i*volTN+j];
				impvol[0] = vol[i*volTN+j] + 0.01;
				 
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
		                            irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);

				mkvega[i*volTN+j] = plainq->getcfValue();
				delete plainq;
  
			}
		}

	 

		 
		//memcpy(mkvega,value2,(size_t)((volSN*volTN)*sizeof(double)));
///// -
//		 
//		 
//		for(i=0; i<volSN; i++) {
//			xval = volSmness[i]*sval;
//			for(j=0; j<volTN; j++) {
//				matu = volBDay[j];
//				cpFlag = mkcpFlag[i*volTN+j];
//				impvol[0] = vol[i*volTN+j] - 0.01;
//				impvolBDay[0] = matu;
//				impvolSmness[0] = volSmness[i]; // =1.0 이나 결과는 같은?
//				  
//				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
//		                            irateN, irateCDay, iRFrate,  
//									divrateN, divrateCDay, divrate, 
//									sval, divN[0], divBDay, divCash, divApy,
//		                            voltype, sval, impvolTN, impvolSN, impvolBDay, impvolSmness, impvol,		  		
//		                            Sminlevel, SmaxlevelC, SmeshM, TmeshN);
//				value2[i*volTN+j] = plainq->getValue();
//				delete plainq;  
//			}
//		}
//
//	  	
// 
//
//		for(i=0; i<volSN*volTN; i++)
//			value2[i] = (mkvega[i] - value2[i])/2.0;
//		

		for(i=0; i<volSN*volTN; i++)
			mkvega[i] = (mkvega[i] - mkdata[i]);

		//memcpy(mkvega,value2,(size_t)((volSN*volTN)*sizeof(double)));
		delete value2;
		delete impvol;
		delete impvolBDay;
		delete impvolSmness;
////////////////////// for vega end

			
 



		memcpy(outvol,vol,(size_t)((volSN*volTN)*sizeof(double)));

		//for(i=0; i<volSN; i++) {
		//	for(j=0; j<volTN; j++) {
		//		outvol[i*volTN+j] = vol[12*volTN+j];
		//	}
		//}
 
///////// start here!! #1

		double *tvalue; // 각 tj 마다 찾을 mk option pr
		int *pidx;
		int targetN;

		tvalue = new double[volSN]; // nonzero vega mk price 저장
		value2 = new double[volSN];
		pidx = new int[volSN]; // nonzero vega index 저장
		voltype = 2; // outvol로 각 grid cy price 계산!
		for(tj=0; tj<volTN; tj++) {
			
 
			//vega 0 = mytol/10 수준으로 설정 찾기
			targetN = 0;
			for(si=0; si<volSN; si++) {
				if(mkvega[si*volTN + tj] >= vegazeromytol) {
					//tvalue[targetN] = mkdata[si*volTN + tj]; 아래의 표현으로 test
					pidx[targetN] = si;
					targetN++; // vega 0 가 아닌 갯수
				} // if
			} // for si

			if(targetN == 0) return maxloop;

			// outvol 준비
			for(si=0; si<targetN; si++) {
				tvalue[si] = mkdata[pidx[si]*volTN + tj];
				outvol[pidx[si]*volTN + tj] = vol[pidx[si]*volTN + tj]; // imp vol 수준의 initial guess vol start!!
			}
			for(i=0; i<volSN; i++) {
				if(i<pidx[0]) {
					outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]*(volSmness[i]-volSmness[pidx[1]]) /(volSmness[pidx[0]]-volSmness[pidx[1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					outvol[i*volTN + tj] += outvol[pidx[1]*volTN+tj]*(volSmness[i]-volSmness[pidx[0]])/(volSmness[pidx[1]]-volSmness[pidx[0]]); // 아래쪽 vega 0 인 부분 ex 채우기!
				}
				if(i>pidx[targetN-1]) {
					//outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 ex 채우기
					outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-2]]) /(volSmness[pidx[targetN-1]]-volSmness[pidx[targetN-2]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					outvol[i*volTN + tj] += outvol[pidx[targetN-2]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-1]]) /(volSmness[pidx[targetN-2]]-volSmness[pidx[targetN-1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
				}
			}
			 			
 
			// guess vol 로 일차 cy 계산하기!!
 
			 
			matu = volBDay[tj];
			for(si=0; si<targetN; si++) {
				xval = volSmness[pidx[si]]*sval;
				cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
									irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
									voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
									Sminlevel, SmaxlevelC, SmeshM, TmeshN);
				cy[si] = plainq->getValue();
				delete plainq;  
   
			}

	  	   
			// LM 시작
			oldss = 0.0;
			for(si=0; si<targetN; si++) 
				oldss = oldss + (tvalue[si]-cy[si])*(tvalue[si]-cy[si]);
					
		
			//pcya 계산
			memcpy(tempvol,outvol,(size_t)((volSN*volTN)*sizeof(double)));					
			for(ai=0; ai<targetN; ai++) {

				memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));			
				outvol[pidx[ai]*volTN + tj] = tempvol[pidx[ai]*volTN + tj] + pa;						
				for(i=0; i<volSN; i++) {
					if(i<pidx[0]) {
						//outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
						outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]*(volSmness[i]-volSmness[pidx[1]]) /(volSmness[pidx[0]]-volSmness[pidx[1]]); // 아래쪽 vega 0 인 부분 ex 채우기!					
						outvol[i*volTN + tj] += outvol[pidx[1]*volTN+tj]*(volSmness[i]-volSmness[pidx[0]])/(volSmness[pidx[1]]-volSmness[pidx[0]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					}
					if(i>pidx[targetN-1]) {
						//outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기					
						outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-2]]) /(volSmness[pidx[targetN-1]]-volSmness[pidx[targetN-2]]); // 아래쪽 vega 0 인 부분 ex 채우기!
						outvol[i*volTN + tj] += outvol[pidx[targetN-2]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-1]]) /(volSmness[pidx[targetN-2]]-volSmness[pidx[targetN-1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					}
				}				
			
				//pcya 계산 하기위함
				  
				 
				for(si=0; si<targetN; si++) {
					xval = volSmness[pidx[si]]*sval;
					cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
					Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
										irateN, irateCDay, iRFrate,  
										divrateN, divrateCDay, divrate, 
										sval, divN[0], divBDay, divCash, divApy,
										voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
										Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   
					value2[si] = plainq->getValue();				
					delete plainq;  
				}

  
				for(si=0; si<targetN; si++)						
					pcya[ai][si] = (value2[si]-cy[si])/pa;				 
						 
			}
			memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
			//// pcya 계산 end!!

			

			for(ai=0; ai<targetN; ai++) {
				beta[ai] = 0.0;
				for(k=0; k<targetN; k++)
					beta[ai] = beta[ai] + (tvalue[k]-cy[k])*pcya[ai][k];
			}

	 

			for(ai=0; ai<targetN; ai++) {
				for(aj=0; aj<targetN; aj++) {
					alpha[ai][aj] = 0.0;
					for(k=0; k<targetN; k++) 
						alpha[ai][aj] = alpha[ai][aj] + pcya[ai][k]*pcya[aj][k];
				}
			}
		
	 

			loopcnt = 0;
			maxerror = 0;
			for(ai=0; ai<targetN; ai++) 
				if(maxerror < abs(tvalue[ai]-cy[ai]))
					maxerror = abs(tvalue[ai]-cy[ai]);
	 
		
			while (loopcnt < maxloop && maxerror > mytol ) {

				for(ai=0; ai<targetN; ai++) {
					for(aj=0; aj<targetN; aj++) {
						if(ai == aj)
							alphad[ai][aj] = alpha[ai][aj]*(1+lambda);
						else
							alphad[ai][aj] = alpha[ai][aj];
					}
				}

		 
				//da = inv(alphad)*beta;
				// alphad * da = beta 풀기 s
				memcpy(tempa, beta,(size_t)((volSN)*sizeof(double)));
				rankA = targetN;
				for(ai=0; ai<targetN; ai++)  {
					maxraw = ai;
					for(aj=ai; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++)
							tempsum = tempsum + alphad[aj][ak]*alphad[ak][ai];
						if((alphad[aj][ai]-tempsum) != 0) {
							maxraw = aj;
							break;
						}
					}

					if(maxraw != ai) {
						for(ak=0; ak<targetN; ak++) {
							tempvalue = alphad[ai][ak];
							alphad[ai][ak] = alphad[maxraw][ak];
							alphad[maxraw][ak] = tempvalue;
						}
						tempvalue = beta[ai];
						beta[ai] = beta[maxraw];
						beta[maxraw] = tempvalue;
					}

					for(aj=ai; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[ai][ak]*alphad[ak][aj];
						alphad[ai][aj] = alphad[ai][aj]-tempsum;
					}

					if(alphad[ai][ai] == 0) {
						rankA = ai;				
						break;
					}

					for(aj=ai+1; aj<targetN; aj++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[aj][ak]*alphad[ak][ai];
						alphad[aj][ai] = (alphad[aj][ai]-tempsum)/alphad[ai][ai];
					}
				} // for(ai.. ) // end lu

				if(rankA == targetN) {
					for(ai=0; ai<targetN; ai++) {
						tempsum = 0;
						for(ak=0; ak<ai; ak++) 
							tempsum = tempsum + alphad[ai][ak]*da[ak];
						da[ai] = beta[ai] - tempsum;
					}

					for(ai=targetN-1; ai>=0; ai--) {
						tempsum = 0;
						for(ak=targetN-1; ak>ai; ak--)
							tempsum = tempsum + alphad[ai][ak]*da[ak];
						da[ai] = (da[ai]-tempsum)/alphad[ai][ai];
					}			
				}
				else {
					loopcnt = 2*maxloop; //
				}
				memcpy(beta, tempa,(size_t)((volSN)*sizeof(double)));
				// alphad * da = beta 풀기 e
		 

				memcpy(tempvol, outvol,(size_t)((volSN*volTN)*sizeof(double)));
				for(ai=0; ai<targetN; ai++) 
					outvol[pidx[ai]*volTN + tj] = outvol[pidx[ai]*volTN + tj] + da[ai];
							
				for(i=0; i<volSN; i++) {
					if(i<pidx[0]) {
						//outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
						outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]*(volSmness[i]-volSmness[pidx[1]]) /(volSmness[pidx[0]]-volSmness[pidx[1]]); // 아래쪽 vega 0 인 부분 ex 채우기!					
						outvol[i*volTN + tj] += outvol[pidx[1]*volTN+tj]*(volSmness[i]-volSmness[pidx[0]])/(volSmness[pidx[1]]-volSmness[pidx[0]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					}
					if(i>pidx[targetN-1]) {
						//outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기					
						outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-2]]) /(volSmness[pidx[targetN-1]]-volSmness[pidx[targetN-2]]); // 아래쪽 vega 0 인 부분 ex 채우기!
						outvol[i*volTN + tj] += outvol[pidx[targetN-2]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-1]]) /(volSmness[pidx[targetN-2]]-volSmness[pidx[targetN-1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
					}
				}				
		 
			 //vector<PlainqRunnable*> runnable;  // multi
		 
			   
				for(si=0; si<targetN; si++) {
					xval = volSmness[pidx[si]]*sval;
					cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
					Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
										irateN, irateCDay, iRFrate,  
										divrateN, divrateCDay, divrate, 
										sval, divN[0], divBDay, divCash, divApy,
										voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
										Sminlevel, SmaxlevelC, SmeshM, TmeshN);
					
					newcy[si] = plainq->getValue();				
					delete plainq;  
   
				}   
	 	
				memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
 			
				newss = 0;
				for(si=0; si<volSN; si++)					
					newss = newss + (tvalue[si] - newcy[si])*(tvalue[si] - newcy[si]);


				if(newss >= oldss)
					lambda = 10*lambda;
				else {
					lambda = 0.1*lambda;

					for(ai=0; ai<targetN; ai++) 				
						outvol[pidx[ai]*volTN + tj] = outvol[pidx[ai]*volTN + tj] + da[ai];
									
					for(i=0; i<volSN; i++) {
						if(i<pidx[0]) {
							//outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
							outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]*(volSmness[i]-volSmness[pidx[1]]) /(volSmness[pidx[0]]-volSmness[pidx[1]]); // 아래쪽 vega 0 인 부분 ex 채우기!					
							outvol[i*volTN + tj] += outvol[pidx[1]*volTN+tj]*(volSmness[i]-volSmness[pidx[0]])/(volSmness[pidx[1]]-volSmness[pidx[0]]); // 아래쪽 vega 0 인 부분 ex 채우기!
						}
						if(i>pidx[targetN-1]) {
							//outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기												
							outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-2]]) /(volSmness[pidx[targetN-1]]-volSmness[pidx[targetN-2]]); // 아래쪽 vega 0 인 부분 ex 채우기!												
							outvol[i*volTN + tj] += outvol[pidx[targetN-2]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-1]]) /(volSmness[pidx[targetN-2]]-volSmness[pidx[targetN-1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
						}
					}		
					/// here !!
					memcpy(cy,newcy,(size_t)((volSN)*sizeof(double)));									
					maxerror = 0;		
					for(ai=0; ai<targetN; ai++) 			
						if(maxerror < abs(tvalue[ai]-cy[ai]))				
							maxerror = abs(tvalue[ai]-cy[ai]);
 		
					oldss = newss;

					memcpy(tempvol,outvol,(size_t)((volSN*volTN)*sizeof(double)));			
					for(ai=0; ai<targetN; ai++) {

						memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));			
						outvol[pidx[ai]*volTN + tj] = tempvol[pidx[ai]*volTN + tj] + pa;						
						for(i=0; i<volSN; i++) {
							if(i<pidx[0]) {
								//outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]; // 아래쪽 vega 0 인 부분 flat 채우기!
								outvol[i*volTN + tj] = outvol[pidx[0]*volTN + tj]*(volSmness[i]-volSmness[pidx[1]]) /(volSmness[pidx[0]]-volSmness[pidx[1]]); // 아래쪽 vega 0 인 부분 ex 채우기!					
								outvol[i*volTN + tj] += outvol[pidx[1]*volTN+tj]*(volSmness[i]-volSmness[pidx[0]])/(volSmness[pidx[1]]-volSmness[pidx[0]]); // 아래쪽 vega 0 인 부분 ex 채우기!
							}
							if(i>pidx[targetN-1]) {
								//outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]; // 위쪽 vega 0 인부분 flat 채우기					
								outvol[i*volTN + tj] = outvol[pidx[targetN-1]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-2]]) /(volSmness[pidx[targetN-1]]-volSmness[pidx[targetN-2]]); // 아래쪽 vega 0 인 부분 ex 채우기!
								outvol[i*volTN + tj] += outvol[pidx[targetN-2]*volTN + tj]*(volSmness[i]-volSmness[pidx[targetN-1]]) /(volSmness[pidx[targetN-2]]-volSmness[pidx[targetN-1]]); // 아래쪽 vega 0 인 부분 ex 채우기!
							}
						}				
			
					 
						for(si=0; si<targetN; si++) {
							xval = volSmness[pidx[si]]*sval;
							cpFlag = mkcpFlag[pidx[si]*volTN + tj];
				  
							Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
												irateN, irateCDay, iRFrate,  
												divrateN, divrateCDay, divrate, 
												sval, divN[0], divBDay, divCash, divApy,
												voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
												Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   					
							value2[si] = plainq->getValue();									
							delete plainq;  
						}
 
  
  
						for(si=0; si<targetN; si++)						
							pcya[ai][si] = (value2[si]-cy[si])/pa;						
						 
					}
					memcpy(outvol, tempvol,(size_t)((volSN*volTN)*sizeof(double)));
				//// pcya 계산 end!!
 

					for(ai=0; ai<targetN; ai++) {
						beta[ai] = 0.0;
						for(k=0; k<targetN; k++)
							beta[ai] = beta[ai] + (tvalue[k]-cy[k])*pcya[ai][k];
					}

					for(ai=0; ai<targetN; ai++) {
						for(aj=0; aj<targetN; aj++) {
							alpha[ai][aj] = 0.0;
							for(k=0; k<targetN; k++) 
								alpha[ai][aj] = alpha[ai][aj] + pcya[ai][k]*pcya[aj][k];
						}
					}


				} // if(newss.. ) .. else
				 
				loopcnt = loopcnt + 1;

			} // while
		
		//memcpy(outvol, vol,(size_t)((volSN*volTN)*sizeof(double)));
		//memcpy(outdata, cy,(size_t)((volSN)*sizeof(double)));
	 
			for(si=0; si<volSN; si++) {
				xval = volSmness[si]*sval;
				cpFlag = mkcpFlag[si*volTN + tj];
				  
				Plainq *plainq = new Plainq(cpFlag,  sval ,  xval, matu,  
									irateN, irateCDay, iRFrate,  
									divrateN, divrateCDay, divrate, 
									sval, divN[0], divBDay, divCash, divApy,
									voltype, sval, volTN, volSN, volBDay, volSmness, outvol,		  		
									Sminlevel, SmaxlevelC, SmeshM, TmeshN);
   
							
				outdata[si*volTN + tj] = plainq->getValue();													
				delete plainq;  
			}
  

		} // for tj..



		delete tvalue;
		delete value2;
		delete pidx;		

///////////////////////////////////////////////////////////////////// end here!! #1

  
	

	delete cy;
	delete newcy;

	delete beta;
	delete da;
 
	for(i=0; i<volSN; i++) {
		delete pcya[i];
		delete alpha[i];
		delete alphad[i];
	}	
	delete pcya;
	delete alpha;
	delete alphad;
 
 
	return loopcnt;
		
} 
  

  