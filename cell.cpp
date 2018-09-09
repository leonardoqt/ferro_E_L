#include "class.h"
#include <cmath>
#include <iostream>

void cell :: init(double length, int NN, pot & dwp)
{	
	const int nei_max = 50;
	int current;
	vec pos;
	vec dp;
	vec p2;
	dp.x[0] = 0;
	dp.x[1] = 0;
	dp.x[2] = dwp.min;
//	dp.x[2] = 0;

	lat = length;
	numx = NN;
	num = NN*NN*NN;
	a_l = new atom[num];
	ind = new int*[num];
	nei = new int*[num];
	num_nei = new int[num];

	for(int t1=0; t1<num; t1++)
	{
		ind[t1] = new int[3];
		nei[t1] = new int[nei_max];
	}

	current = 0;
	for(int t1=0; t1<NN; t1++)
	for(int t2=0; t2<NN; t2++)
	for(int t3=0; t3<NN; t3++)
	{
		// initiate atom position and dipole
		pos.x[0] = t1 * lat;
		pos.x[1] = t2 * lat;
		pos.x[2] = t3 * lat;
		p2 = pos + dp;
		a_l[current].get_pos0(pos);
		a_l[current].get_pos(p2);
		a_l[current].get_dipole();
		// register index of each atom
		ind[current][0] = t1;
		ind[current][1] = t2;
		ind[current][2] = t3;
		current++;
	}
//	ene_onsite0 = ene_dipole0 = ene_tot0 = 0;
}

void cell :: get_neighbor(double max_d)
{
	vec xx,yy,zz;

	xx.clean(), xx.x[0] = lat*numx;
	yy.clean(), yy.x[1] = lat*numx;
	zz.clean(), zz.x[2] = lat*numx;

	for(int t1=0; t1<num; t1++)
	{
		num_nei[t1]=0;
		for(int t2=0; t2<num; t2++)
			if(t2!=t1)
			{
				for(int nx=-1;nx<=1;nx++)
				for(int ny=-1;ny<=1;ny++)
				for(int nz=-1;nz<=1;nz++)
				{
					if ((a_l[t2].pos0+xx*nx+yy*ny+zz*nz-a_l[t1].pos0).norm() <= max_d)
					{
						nei[t1][num_nei[t1]] = t2;
						num_nei[t1]++;
					}
				}
			}
	}
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
	double re_s, im_s, k2, dk;
	const double pi = 3.141592653589793238462643383279502884;
	const int n_max=2, k_max=2;
	const double sigma = 3.6;
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
	//
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
	dk = 2*pi / lat / numx;
	e_temp = 0;
	for(int t1=-k_max; t1<=k_max; t1++)
	for(int t2=-k_max; t2<=k_max; t2++)
	for(int t3=-k_max; t3<=k_max; t3++)
	{
		LL.x[0]=t1*dk, LL.x[1]=t2*dk, LL.x[2]=t3*dk;
		k2 = LL.norm() * LL.norm();
		re_s = im_s = 0;
		if(t1!=0 || t2!=0 || t3!=0)
		{
			for(int nn=0; nn<num; nn++)
			{
				re_s += (cos(LL*a_l[nn].pos) - cos(LL*a_l[nn].pos0));
				im_s += (sin(LL*a_l[nn].pos) - sin(LL*a_l[nn].pos0));
			}
			e_temp += exp(-sigma*sigma*k2/2)/k2*(re_s*re_s+im_s*im_s);
		}
	}
	e_temp /= (2*pow(lat*numx,3));
//	std::cout<<"old E= "<<e_temp<<std::endl;
	//
	a_l[n0].pos = a_l[n0].pos + d_new_pos;
	de_long = 0;
	for(int t1=-k_max; t1<=k_max; t1++)
	for(int t2=-k_max; t2<=k_max; t2++)
	for(int t3=-k_max; t3<=k_max; t3++)
	{
		LL.x[0]=t1*dk, LL.x[1]=t2*dk, LL.x[2]=t3*dk;
		k2 = LL.norm() * LL.norm();
		re_s = im_s = 0;
		if(t1!=0 || t2!=0 || t3!=0)
		{
			for(int nn=0; nn<num; nn++)
			{
				re_s += (cos(LL*a_l[nn].pos) - cos(LL*a_l[nn].pos0));
				im_s += (sin(LL*a_l[nn].pos) - sin(LL*a_l[nn].pos0));
			}
			de_long += exp(-sigma*sigma*k2/2)/k2*(re_s*re_s+im_s*im_s);
		}
	}
	de_long /= (2*pow(lat*numx,3));
//	std::cout<<"new E= "<<de_long<<std::endl;
	a_l[n0].pos = a_l[n0].pos - d_new_pos;
	de_long -= e_temp;

	return de_onsite+de_short+de_long;
}

double cell :: get_d_ene_short(pot & dwp, int n0, vec & d_new_pos, double lambda)
{
	double de_onsite, de_short;
	double e_temp;
	double norm_temp;
	// change in onsite energy
	e_temp = dwp.get_E(a_l[n0].dipole.norm());
	de_onsite = dwp.get_E((a_l[n0].dipole+d_new_pos).norm()) - e_temp;
	// dipole-dipole interaction
	// change in short range energy
	e_temp = 0;
	for(int t1=0; t1<num_nei[n0]; t1++)
		e_temp += lambda*(a_l[n0].dipole*a_l[nei[n0][t1]].dipole);
	a_l[n0].dipole = a_l[n0].dipole + d_new_pos;
	de_short = 0;
	for(int t1=0; t1<num_nei[n0]; t1++)
		de_short += lambda*(a_l[n0].dipole*a_l[nei[n0][t1]].dipole);
	de_short -= e_temp;
	a_l[n0].dipole = a_l[n0].dipole - d_new_pos;

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
