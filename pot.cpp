#include <cmath>
#include "class.h"


void pot :: init(double E, double L)
{
	double aa,bb;
	ene = E;
	min = L;
	aa = ene / pow(min,4);
	bb = 2*ene / pow(min,2);
	x2 = y2 = -aa;
	xy = 0;
	x4 = y4 = bb;
	x2y2 = 2*bb;
	
}

void pot :: init(double a0,double a1,double a2,double a3, double a4, double a5)
{
	x2=a0,y2=a1,xy=a2,x4=a3,y4=a4,x2y2=a5;
}

double pot :: get_E(vec &p)
{
	return x2*pow(p.x[0],2)+y2*pow(p.x[1],2)+xy*p.x[0]*p.x[1]+x4*pow(p.x[0],4)+y4*pow(p.x[1],4)+x2y2*pow(p.x[0]*p.x[1],2);
}
