#include <cmath>
#include "class.h"


void pot :: init(double E, double L)
{
	ene = E;
	min = L;
	aa = ene / pow(min,4);
	bb = 2*ene / pow(min,2);
}

double pot :: get_E(double r)
{
	return aa*pow(r,4) - bb*pow(r,2);
}
