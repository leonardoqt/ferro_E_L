#include "class.h"
#include <cmath>
#include <iostream>

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

void cell :: update_pos(int n0, vec &  d_new_pos)
{
	a_l[n0].pos = a_l[n0].pos + d_new_pos;
	a_l[n0].dipole = a_l[n0].pos - a_l[n0].pos0;
}

double cell :: get_d_ene(pot & dwp, int n0, vec & d_new_pos)
{
	double de_onsite, de_short, de_long;
	double e_temp;
	double norm_temp;
	const int n_max=1, k_max=2;
	const double sigma = 1;
	vec pos_old;
	vec LL;
	// change in onsite energy
	e_temp = dwp.get_E(a_l[n0].dipole.norm());
	de_onsite = dwp.get_E((a_l[n0].dipole+d_new_pos).norm()) - e_temp;
	// dipole-dipole interaction
	// change in short range energy
	e_temp = 0;
	for(int t1=-n_max; t1<=n_max; t1++)
	for(int t2=-n_max; t2<=n_max; t2++)
	for(int t3=-n_max; t3<=n_max; t3++)
	{
		LL.x[0]=t1*lat*numx, LL.x[1]=t2*lat*numx, LL.x[2]=t3*lat*numx;
		for(int nn=0; nn<num; nn++)
		{
			if(t1!=0 || t2!=0 || t3!=0 || nn!=n0)	
			{
				norm_temp = (a_l[n0].pos - a_l[nn].pos0 + LL).norm();
				e_temp -= erfc(norm_temp/sqrt(2)/sigma)/norm_temp;
				norm_temp = (a_l[n0].pos - a_l[nn].pos + LL).norm();
				e_temp += erfc(norm_temp/sqrt(2)/sigma)/norm_temp;
			}
		}
	}
//	std::cout<<"old E= "<<e_temp<<std::endl;
	a_l[n0].pos = a_l[n0].pos + d_new_pos;
	de_short = 0;
	for(int t1=-n_max; t1<=n_max; t1++)
	for(int t2=-n_max; t2<=n_max; t2++)
	for(int t3=-n_max; t3<=n_max; t3++)
	{
		LL.x[0]=t1*lat*numx, LL.x[1]=t2*lat*numx, LL.x[2]=t3*lat*numx;
		for(int nn=0; nn<num; nn++)
		{
			if(t1!=0 || t2!=0 || t3!=0 || nn!=n0)
			{
				norm_temp = (a_l[n0].pos - a_l[nn].pos0 + LL).norm();
				de_short -= erfc(norm_temp/sqrt(2)/sigma)/norm_temp;
				norm_temp = (a_l[n0].pos - a_l[nn].pos + LL).norm();
				de_short += erfc(norm_temp/sqrt(2)/sigma)/norm_temp;
			}
		}
	}
//	std::cout<<"new E= "<<de_short<<std::endl;
	a_l[n0].pos = a_l[n0].pos - d_new_pos;
	de_short -= e_temp;
	// change in long range energy
	e_temp = 0;

	return de_short;
}
