#include <cmath>
#include <cstdlib>
#include "class.h"
#include <iostream>
#include <fstream>

void cell :: init(double length, int NN, pot & dwp)
{	
	const int nei_max = 50;
	int current;
	vec pos;
	vec dp;
	vec p2;
	double l_max = dwp.x2*dwp.x4*sqrt(fabs(dwp.x2/2/dwp.x4))/fabs(dwp.x2*dwp.x4);
	dp.x[2] = 0;

	lat = length;
	numx = NN;
	num = NN*NN;
	a_l = new atom[num];
	nei = new int*[num];
	num_nei = new int[num];

	for(int t1=0; t1<num; t1++)
	{
		nei[t1] = new int[nei_max];
	}

	current = 0;
	for(int t1=0; t1<NN; t1++)
	for(int t2=0; t2<NN; t2++)
	{
		// initiate atom position and dipole
		pos.x[0] = t1 * lat;
		pos.x[1] = t2 * lat;
		pos.x[2] = 0;
		dp.x[0] = ((rand()/(double)RAND_MAX)*2-1)*l_max;
		dp.x[1] = ((rand()/(double)RAND_MAX)*2-1)*l_max;
		p2 = pos + dp;
		a_l[current].get_pos0(pos);
		a_l[current].get_pos(p2);
		a_l[current].get_dipole();
		current++;
	}
}

void cell :: get_neighbor(double max_d)
{
	vec xx,yy;

	xx.clean(), xx.x[0] = lat*numx;
	yy.clean(), yy.x[1] = lat*numx;

	for(int t1=0; t1<num; t1++)
	{
		num_nei[t1]=0;
		for(int t2=0; t2<num; t2++)
			if(t2!=t1)
			{
				for(int nx=-1;nx<=1;nx++)
				for(int ny=-1;ny<=1;ny++)
				{
					if ((a_l[t2].pos0+xx*nx+yy*ny-a_l[t1].pos0).norm() <= max_d)
					{
						nei[t1][num_nei[t1]] = t2;
						num_nei[t1]++;
					}
				}
			}
	}
}

void cell :: update_pos(int n0, vec & d_pos)
{
	a_l[n0].pos = a_l[n0].pos + d_pos;
	a_l[n0].dipole = a_l[n0].pos - a_l[n0].pos0;
}

double cell :: get_d_ene_cross(pot & dwp, int n0, vec & d_pos, double lambda)
{
	double de_onsite, de_short;
	double e_temp;
	double pi_x,pi_y;
	// change in onsite energy
	e_temp = dwp.get_E(a_l[n0].dipole);
	a_l[n0].dipole = a_l[n0].dipole+d_pos;
	de_onsite = dwp.get_E(a_l[n0].dipole) - e_temp;
	a_l[n0].dipole = a_l[n0].dipole-d_pos;
	// dipole-dipole interaction
	// change in short range energy
	pi_x = pi_y = 0;
	for(int t1=0; t1<num_nei[n0]; t1++)
	{
		pi_x += a_l[nei[n0][t1]].dipole.x[0];
		pi_y += a_l[nei[n0][t1]].dipole.x[1];
	}
	de_short = lambda*(d_pos.x[0]*pi_y+d_pos.x[1]*pi_x);
	return de_onsite+de_short;
}

vec cell :: find_dipole()
{
	vec res;
	res.clean();
	for (int t1=0; t1<num; t1++)
		res = res + a_l[t1].dipole;
	return res/num;
}

std::ofstream& operator<<(std::ofstream& out,cell& sys1)
{
	for (int t1=0; t1<sys1.num; t1++)
		out<<sys1.a_l[t1].pos0<<'\t'<<sys1.a_l[t1].dipole<<std::endl;
	return out;
}
