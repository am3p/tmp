#include "Basic.h"
#include "odysseyRandom.h"

 

#pragma once
using namespace std;
 
class StepDownBMC : public Basic
{
public:
	StepDownBMC(long NAST, long Nst, long simi, long sidx, long rndidx);

	virtual ~StepDownBMC();

	 
	double getValue();
	 
	 
private:	

	long m_sidx;
	long m_rndidx;
	
	long m_simi;
	long m_noiters;
	long nostep;
	long NASSETS;

	odysseyRandom *rnd;

	
	 
	double cf_price();
	double mc_price();
	 
};


