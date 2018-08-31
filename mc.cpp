#include <cstdlib>
#include <cmath>
#include "class.h"

void mc :: init(double Scale, int Check_scale, double Temperature)
{
	scale = Scale;
	check_scale = Check_scale;
	accept_scale = 0;
	T = Temperature;
}

void mc :: update_scale()
{
	if (accept_scale/check_scale < 0.2)
		scale = scale / 2;
	else if (accept_scale/check_scale > 0.4)
		scale = scale * 1.616;
	
	accept_scale = 0;
}

void mc :: mv_atm(int atom_num, int &target, vec & d_dipole)
{
	double len[3];
	target=rand()%atom_num;
	len[0] = ((rand()/(double)RAND_MAX)*2-1)*scale;
	len[1] = ((rand()/(double)RAND_MAX)*2-1)*scale;
	len[2] = ((rand()/(double)RAND_MAX)*2-1)*scale;
	d_dipole = len;
}

int mc :: if_accept(double dE)
{
	if (dE < 0)
		return 1;
	else if (rand()/(double)RAND_MAX < exp(dE/T))
		return 1;
	else
		return 0;
}
