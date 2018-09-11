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

void mc :: update_scale(pot& dwp)
{
	double l_max = sqrt(dwp.ax2/2/dwp.bx4);
	if (accept_scale/check_scale < 0.2)
		scale = scale / 2;
	else if (accept_scale/check_scale > 0.4)
		scale = scale * 1.616;
	accept_scale = 0;
	if (scale > l_max*4)
		scale = l_max*4;
}

void mc :: mv_atm(cell &sys1, int &target, vec & d_dipole)
{
	double len[3];
	target=rand()%sys1.num;
	len[0] = 0;//((rand()/(double)RAND_MAX)*2-1)*scale;
	len[1] = 0;//((rand()/(double)RAND_MAX)*2-1)*scale;
	if (rand()/(double)RAND_MAX < 0.5)
		len[2] = ((rand()/(double)RAND_MAX)*2-1)*scale;
	else 
		len[2] = -2*sys1.a_l[target].dipole.x[2] + ((rand()/(double)RAND_MAX)*2-1)*scale;
	d_dipole = len;
}

void mc :: mv_atm_2d(cell &sys1, int &target, vec & d_dipole)
{
	double len[3];
	target=rand()%sys1.num;
	len[0] = 0;
	len[1] = 0;
	len[2] = 0;
	if (rand()/(double)RAND_MAX < 0.5)
		len[0] = -2*sys1.a_l[target].dipole.x[0];
	if (rand()/(double)RAND_MAX < 0.5)
		len[1] = -2*sys1.a_l[target].dipole.x[1];
	len[0] += ((rand()/(double)RAND_MAX)*2-1)*scale;
	len[1] += ((rand()/(double)RAND_MAX)*2-1)*scale;
	d_dipole = len;
}

int mc :: if_accept(double dE)
{
	if (dE < 0)
		return 1;
	else if (rand()/(double)RAND_MAX < exp(-dE/T))
		return 1;
	else
		return 0;
}
