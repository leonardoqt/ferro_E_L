#include <iostream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "class.h"

using namespace std;

int main()
{
	srand(time(0));
	pot dvp;
	cell sys1;
	mc run1;
	double barrier = 0.01, x0 = 0.1, length = 1.0, temperature = 0.1;
	int num_x = 10,num_tot;
	int tot_run = 100000;
	int check_scale = 100;
	vec d_dipole;
	int atom_to_mv;
	double dE;
/*
	double temp;
	double dd[3]={0.0,0.0,0.1};
	vec ddd;

	dvp.init(barrier,x0);
	sys1.init(length,num_x);
	
	for(int t1=0; t1<num_x*num_x*num_x; t1++)
	{
		dd[0] = rand()/(double)RAND_MAX * 0.1;
		dd[1] = rand()/(double)RAND_MAX * 0.1;
		dd[2] = rand()/(double)RAND_MAX * 0.1;
		ddd = dd;
		sys1.update_pos(t1, ddd);
	}
	dd[0] = dd[1] = dd[2] = 0.1;
*/

	// get potential, cell parameters
	cin>>barrier>>x0>>length>>num_x>>temperature;
	cin>>tot_run>>check_scale;
	num_tot = num_x*num_x*num_x;

	// initiate potential, cell, and mc
	dvp.init(barrier,x0);
	sys1.init(length,num_x,dvp);
	run1.init(x0,check_scale,temperature);

	// start mc
	for(int t1=1; t1<tot_run; t1++)
	{
		run1.mv_atm(num_tot,atom_to_mv,d_dipole);
		dE = sys1.get_d_ene(dvp,atom_to_mv,d_dipole);
		if(run1.if_accept(dE))
		{
			sys1.update_pos(atom_to_mv,d_dipole);
			run1.accept_scale += 1;
		}
		// update scale
		if(t1 % run1.check_scale == 0)
		{
			run1.update_scale();
//			cout<<"new_scale: "<<run1.scale<<endl;
		}
//		cout<<dE<<'\t'<<run1.accept_scale<<endl;
		if(t1 % 10 == 0)
		{
			cout<<t1<<'\t'<<sys1.find_dipole()<<'\t'<<run1.scale<<endl;
//			sys1.find_dipole().print();
		}

	}
	return 0;
}
