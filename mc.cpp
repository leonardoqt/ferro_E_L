#include "class.h"

void mc :: mc_init(double Scale, int Check_scale)
{
	scale = Scale;
	check_scale = Check_scale;
	accept_scale = 0;
}

void mc :: update_scale()
{
	if (accept_scale/check_scale < 0.2)
		scale = scale / 2;
	else if (accept_scale/check_scale > 0.4)
		scale = scale * 1.616;
	
	accept_scale = 0;
}

void mc :: mv_atm(int atom_num, int &target, double *&vec)
{
	int t1;
	double res;
	target=rand()%atom_num;
	for (t1=0; t1<DIM; t1++)
	{
		res = rand()%200001 - 100000;
		vec[t1] = scale*res/100000;
	}
}
