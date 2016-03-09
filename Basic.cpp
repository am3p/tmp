#include "Basic.h"

Basic::Basic()
{

}

Basic::~Basic()
{

}

double Basic::cnorm(double x)
{
		
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;

	double a = fabs(x);
	double t = 1.0/(1.0 + a*p);
	double b = c2*exp((-x)*(x/2.0));
	double out = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;

	
	if(x > 6.0) {return 1.0;}
	if(x < -6.0) {return 0.0;}

	out = 1.0 - b*out;
	if(x < 0.0) out = 1.0 - out;
	
	return out;


}


//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
double Basic::interp1(double x, double *Dx, double *Dy, int N, int meshType, int interpType)
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
double Basic::interp1(int x, int *Dx, double *Dy, int N, int meshType, int interpType)
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



int  Basic::IsinArr(int x, int N, int *Arr) 
{
	int i;
	for(i=0; i<N; i++) 
		if(x == Arr[i]) 
			return i+1;

	return 0;
}


//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
double Basic::interp2(double x, double y, double *Dx, double *Dy, double **Dz, int Nx, int Ny) // Nx는 x의 분할 개수로 Dx는 0~Nx까지 Nx+1개의 Data임
{
	double val;
	int idx, idy, i, j;

	for(i=0; i<=Nx; i++) {
		if(x<Dx[i]) {
			idx = i-1;
			break;
		}
	}
	if(x<Dx[0]) idx = 0;
	if(x>=Dx[Nx]) idx = Nx-1;

	for(j=0; j<=Ny; j++) {
		if(y<Dy[j]) {
			idy = j-1;
			break;
		}
	}
	if(y<Dy[0]) idy = 0;
	if(y>=Dy[Ny]) idy = Ny-1;
			
	val = (x-Dx[idx+1])*(y-Dy[idy+1])/((Dx[idx]-Dx[idx+1])*(Dy[idy]-Dy[idy+1])) * Dz[idx][idy];	
	val = val + (x-Dx[idx])*(y-Dy[idy+1])/((Dx[idx+1]-Dx[idx])*(Dy[idy]-Dy[idy+1])) * Dz[idx+1][idy];	
	val = val + (x-Dx[idx+1])*(y-Dy[idy])/((Dx[idx]-Dx[idx+1])*(Dy[idy+1]-Dy[idy])) * Dz[idx][idy+1];	
	val = val + (x-Dx[idx])*(y-Dy[idy])/((Dx[idx+1]-Dx[idx])*(Dy[idy+1]-Dy[idy])) * Dz[idx+1][idy+1];		 
	 
	return val;
}


//(target x, Xpoints, Yvalue, Data#, (0:균등, 1:비균등), (0:flat, 1:두점 선형))
double Basic::interp2(int x, int y, int *Dx, int *Dy, double **Dz, int Nx, int Ny) // Nx는 x의 분할 개수로 Dx는 0~Nx까지 Nx+1개의 Data임
{
	double val;
	int idx, idy, i, j;

	for(i=0; i<=Nx; i++) {
		if(x<Dx[i]) {
			idx = i-1;
			break;
		}
	}
	if(x<Dx[0]) idx = 0;
	if(x>=Dx[Nx]) idx = Nx-1;

	for(j=0; j<=Ny; j++) {
		if(y<Dy[j]) {
			idy = j-1;
			break;
		}
	}
	if(y<Dy[0]) idy = 0;
	if(y>=Dy[Ny]) idy = Ny-1;
			
	val = (x-Dx[idx+1])*(y-Dy[idy+1])/((Dx[idx]-Dx[idx+1])*(Dy[idy]-Dy[idy+1])) * Dz[idx][idy];	
	val = val + (x-Dx[idx])*(y-Dy[idy+1])/((Dx[idx+1]-Dx[idx])*(Dy[idy]-Dy[idy+1])) * Dz[idx+1][idy];	
	val = val + (x-Dx[idx+1])*(y-Dy[idy])/((Dx[idx]-Dx[idx+1])*(Dy[idy+1]-Dy[idy])) * Dz[idx][idy+1];	
	val = val + (x-Dx[idx])*(y-Dy[idy])/((Dx[idx+1]-Dx[idx])*(Dy[idy+1]-Dy[idy])) * Dz[idx+1][idy+1];		 
	 
	return val;
}

 