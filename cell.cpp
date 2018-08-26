#include "class.h"

void cell :: init(double length, int NN)
{	
	vec pos;
	lat = length;
	numx = NN;
	num = NN*NN*NN;
	a_l = new atom[num];

	for(int t1=0; t1<NN; t1++)
	for(int t2=0; t2<NN; t2++)
	for(int t3=0; t3<NN; t3++)
	{
		pos.x[0] = t1 * lat;
		pos.x[1] = t2 * lat;
		pos.x[2] = t3 * lat;
		a_l[t1*NN*NN+t2*NN+t3].get_pos0(pos);
		a_l[t1*NN*NN+t2*NN+t3].get_pos(pos);
		a_l[t1*NN*NN+t2*NN+t3].get_dipole();
	}
	ene_onsite0 = ene_dipole0 = ene_tot0 = 0;
}

void cell :: get_ene(pot & dwp)
{
	// onsite energy
	ene_onsite1 = 0;
	for(int t1=0; t1<num; t1++)
		ene_onsite1 += dwp.get_E(a_l[t1].dipole.norm());
	// dipole-dipole interaction
}
