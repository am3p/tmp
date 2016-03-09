#ifdef ODYSSEY_EXPORTS
#define ODYSSEY_API   __declspec(dllexport)
#else
#define ODYSSEY_API   __declspec(dllimport)
#endif
 
double ODYSSEY_API  __stdcall PlainFn(int uAssetN, int cpFlag, double prate,
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
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN);
 
double ODYSSEY_API  __stdcall PlainqFn( int cpFlag,  double *sval,  double *xval, int *matu, 
									  int irateN,  int *irateCDay, double *iRFrate,  
									  int divrateN, int *divrateCDay, double *divrate, 
									  double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,									  
									  double *Ogreeks,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN);

 
double ODYSSEY_API  __stdcall DigitalFn(int uAssetN, int cpFlag, double *sval, double *bval, double *xval, int xvaltype,
	                                  double *Cpn, // Digital option 액면 1원(계약 수 1) 기준의 실제 지급하는 금액 예: 3% 지급일 때 0.03 입력!
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
									  double shiftRate, double spreadRate);

 
double ODYSSEY_API  __stdcall BarrierFn(int uAssetN, int duFlag, int ioFlag, int cpFlag, int knockFlag, double prate,
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
									  double shiftRate, double spreadRate);


double ODYSSEY_API  __stdcall BinaryBarrierFn(int uAssetN, int duFlag, int ioFlag, int cpFlag, int knockFlag,
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
									  double shiftHRate, double shiftXRate, double spreadHRate, double spreadXRate);


double ODYSSEY_API  __stdcall StepDownAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 double ODYSSEY_API  __stdcall TwoAssetDigitalFn(int uAssetN, int cpFlag, double *sval, double *bval, double *xval, int xvaltype, 
									  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN,
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int corrtype, int corrN, int *corrBDay, double *corr,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx, int corrdeltaN, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas1, double *Otvannas1, double *Otzommas1,
									  double *Otvegas2, double *Otvannas2, double *Otzommas2,
									  double *Otrfrhos, double *Otdcrhos,
									  double *Otcorrdelta,
									  int batchFlag, char *kCode,
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN,
									  double shiftRate, double spreadRate);
 
 
 double ODYSSEY_API  __stdcall StepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);



 

double ODYSSEY_API  __stdcall CliquetAFnCliquetAFn(int uAssetN, double *sval, double *bval, double prate, double GFloor, double LCap, double LFloor, double NAQ,	                                 
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
									  double Xmaxlevel, double Zmesh, int XmeshM, int *TmeshN);




double ODYSSEY_API  __stdcall ZeroBondFnZeroBondFn(int matuN, int *matu, 									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,									 
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos);



double ODYSSEY_API  __stdcall FRNFn(double ctime,									 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate,  // cd
									  int irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,									 
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, double *Otvannas, double *Otzommas,
									  double *Otrfrhos, double *Otdcrhos);



double ODYSSEY_API  __stdcall StepDownCDAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);

 
 double ODYSSEY_API  __stdcall StepDownCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


  
 
 //For MC
 int ODYSSEY_API  _stdcall odysseyRandom_INIT(long SimulNo);//, long iMAXAST);
 
 int ODYSSEY_API  _stdcall odysseyRandom_END(long SimulNo);//, long iMAXAST);

  

 
 double ODYSSEY_API  __stdcall HStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);

 
double ODYSSEY_API  __stdcall RStepUpAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate);




double ODYSSEY_API  __stdcall StepDownCPAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  int batchFlag, char *kCode,
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate);




double ODYSSEY_API  __stdcall RStepUpCDAFn(int *uAssetN, int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);





double ODYSSEY_API  __stdcall StepDownCPCDAFn(int *uAssetN, int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);

 


 double ODYSSEY_API  __stdcall StepDownCPBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
 double ODYSSEY_API  __stdcall StepDownCPCDBFn(int *uAssetN,  int knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
 double ODYSSEY_API  __stdcall StepDownCTBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
 
 double ODYSSEY_API  __stdcall StepDownCTCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
 double ODYSSEY_API  __stdcall StepDownBCBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);

 
 
 double ODYSSEY_API  __stdcall StepDownBCCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
double ODYSSEY_API  __stdcall  StabilityAFn(int *uAssetN, int knockFlag, double *sval, double erPerf, double *xval, double multiplier, 
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
									  double SLimit, long SimulNo);




 double ODYSSEY_API  __stdcall AirbagStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  long SimulNo);

 
 double ODYSSEY_API  __stdcall LizardStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  long SimulNo);


 

  double ODYSSEY_API  __stdcall BarrierBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


  
  double ODYSSEY_API  __stdcall PlainBFn(int *uAssetN,  double *sval, double *bval, double *xval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


  //3S
  
 double ODYSSEY_API  __stdcall StepDownCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos, long SimulNo);
 


 double ODYSSEY_API  __stdcall StepDownCDCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos, long SimulNo);

 
 double ODYSSEY_API  __stdcall StepDownCPCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int*matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int *irateN, int *irateBDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int *voltype, double *volbval, int *volTN, int *volSN, int *volBDay, double *volSmness, double *vol,
									  int *corrtype, int *corrN, int *corrBDay, double *corr,
									  int greektype,
									  int *termgrNs, int *vegaidx, int *rhoidx, int *corrdeltaidx,  
									  double *Ogreeks,
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos, long SimulNo);
 


 double ODYSSEY_API  __stdcall StepDownCDCPCFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double *Otvegas, double *Otrfrhos, double *Otdcrhos, long SimulNo);




 
 
double ODYSSEY_API  __stdcall  TomAFn(int *uAssetN, double *sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,	                                   
	                                  int matuN, int *matu, double *prate, double ctime,									   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int greektype,
									  int vegaN, int *vegaidx, int rhoN, int *rhoidx,
									  double *Ogreeks,
									  double *Otvegas, 
									  double *Otrfrhos, double *Otdcrhos,
									  long SimulNo, int batchFlag, char *kCode);




double ODYSSEY_API  __stdcall TomCDAFn(int *uAssetN, double *sval, double stom, double btom, int STDayN, int *STStartDay, double *IXT, double *IndexT, int *STEndDay, double *IXTheta,	                                   
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
									  long SimulNo, int batchFlag, char *kCode);





  
  double ODYSSEY_API  __stdcall MakeHTAFn(char *kCode, int dataN, double baseval, int rangeN, double *Oresult);

  	
  double Sinterp1(double x, double *Dx, double *Dy, int N, int meshType, int interpType);
  double Sinterp1(int x, int *Dx, double *Dy, int N, int meshType, int interpType);
	  


  
double ODYSSEY_API  __stdcall VBStepDownAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double Sminlevel, double Smaxlevel, int *SmeshM, int *TmeshN,
									  double shiftHRate, double spreadHRate, double *spreadXRate);


double ODYSSEY_API  __stdcall VBStepDownCDAFn(int *uAssetN, int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  int CDN, double CDadjval, int *ObBDay, int *PayBDay, double *FlegFactor, double *ObCDrate, // cd
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 double ODYSSEY_API  __stdcall VBStepDownBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);


 
 double ODYSSEY_API  __stdcall VBStepDownCDBFn(int *uAssetN,  int *knockFlag, double *sval, double *bval, double *xval, double *hval, int xvaltype, 
	                                  double *Cpn,
	                                  int *matuN, int *matu, double ctime,
									  int maxevalN, int *evalN, 
									  double *CDlegInfo,
									  int irateN, int *irateBDay,  double *iRFrate, double *iDCrate, double *iIRSrate,
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
									  double shiftHRate, double spreadHRate, double *spreadXRate);

 
 

 
double ODYSSEY_API  __stdcall  mcPlainFn(int cpFlag, double *sval, double *xval, 	                                  
	                                  int matuN, int *matu,							   
									  int irateN, int *irateBDay, int *irateCDay, double *iRFrate, double *iDCrate, double *iIRSrate,
									  double *divrate, double *divbval, int *divN, int *divBDay, double *divCash, double divApy,
									  int voltype, double *volbval, int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									 long SimulNo);



 
int ODYSSEY_API  __stdcall SVolTableCal(double sval,  
									  int irateN,  int *irateCDay, double *iRFrate,  
									  int divrateN, int *divrateCDay, double *divrate, 
									  int *divN, int *divBDay, double *divCash, double divApy,
									  int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int *mkcpFlag, double *mkdata, double *mkvega, double *outvol, double *outdata,									  
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN, int maxloop, double *inmytol);


 
int ODYSSEY_API  __stdcall SVolTableCalCF(double sval,  
									  int irateN,  int *irateCDay, double *iRFrate,  
									  int divrateN, int *divrateCDay, double *divrate, 
									  int *divN, int *divBDay, double *divCash, double divApy,
									  int volTN, int volSN, int *volBDay, double *volSmness, double *vol,
									  int *mkcpFlag, double *mkdata, double *mkvega, double *outvol, double *outdata,									  
									  double Sminlevel, double *SmaxlevelC, int *SmeshM, int *TmeshN, int maxloop, double *inmytol);

      