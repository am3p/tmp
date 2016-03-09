#include "StepDownBMC.h"
#include <iostream>
#include <fstream>

 
extern odysseyRandom *opr;

//#define NASSETS 2L


StepDownBMC::StepDownBMC(long NAST, long Nst, long simi, long sidx, long rndidx)
		   
{		
	 
	NASSETS = NAST;
	nostep = Nst;
	m_sidx = sidx;
	m_rndidx = rndidx;
	m_simi = simi;
 
	rnd = opr;

	m_noiters = rnd->SimulNo;

	  
}

StepDownBMC::~StepDownBMC()
{	
 
}

double StepDownBMC::getValue() {
 
	 
	double value;
	 
 
	value = mc_price();
	 
	 

	return value;

}
 

double StepDownBMC::cf_price()
{	
	return 7777;
}
    


double StepDownBMC::mc_price()
{	
	double value;
	long i;

	
	long   *rsequns;
	double **randoms;   
 
	//nostep = 252;

				
	rsequns =(long *)malloc((size_t)((nostep)*sizeof(long ))); 
	randoms =(double **)malloc((size_t)((NASSETS)*sizeof(double *))); 
	 

	
	for ( i = 0; i < NASSETS; i++ ) {
		randoms[i] = (double *)malloc((size_t)((m_noiters)*sizeof(double)));
		 
	}

	if ( rnd->randnum(1,m_noiters,NASSETS,randoms,0) != 0 ) return(-6);

	//for sim
	if(rnd->randseq(1,m_simi,nostep,rsequns)!=0)  return(-6);

	value = randoms[m_sidx][rsequns[m_rndidx]];

	
	/* free of memory */
	for (int k = 0; k < NASSETS; k++ )  { 
		free(randoms[k]); 	
		 
	} 
	free(randoms); 
	free(rsequns);	
	 




	return value;
}
    