#include <cmath>
#include "class.h"


void pot :: init(double E, double L)
{
	ene = E;
	min = L;
	aa = ene / pow(min,4);
	bb = 2*ene / pow(min,2);
}

void pot :: init_2d(double a, double b, double c)
{
	ax2 = a, axy = b, bx4 = c;
}

double pot :: get_E(double r)
{
	return aa*pow(r,4) - bb*pow(r,2);
}

double pot :: get_E_2d(vec &p)
{
	return -ax2*(pow(p.x[0],2)+pow(p.x[1],2)) - 2*axy*p.x[0]*p.x[1] + bx4*(pow(p.x[0],4)+pow(p.x[1],4));
}
