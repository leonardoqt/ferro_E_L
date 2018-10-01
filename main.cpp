#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "class.h"

using namespace std;

int main()
{
	srand(time(0));
	ofstream output;
	pot dwp;
	cell sys1;
	mc run1;
	double barrier,p_min;
	double length = 1.0, max_length_nei = 1.1, lambda = 0.1;
	int num_x = 10;
	double temperature = 0.01;
	int tot_run = 100000, dump_factor = 100;
	int check_scale = 100;

	vec d_dipole;
	int atom_to_mv;
	double dE;

	output.open("last_dipole.dat");

	// get potential, cell parameters
	cin>>barrier>>p_min;
	cin>>length>>num_x>>max_length_nei>>lambda;
	cin>>temperature;
	cin>>tot_run>>dump_factor>>check_scale;

	// initiate potential, cell, and mc
	dwp.init(barrier,p_min);
	sys1.init(length,num_x,dwp);
	output<<sys1<<endl;
	sys1.get_neighbor(max_length_nei);
	run1.init(length/5,check_scale,temperature);

	// start mc
	for(int t1=1; t1<tot_run; t1++)
	{
		run1.mv_atm(sys1,atom_to_mv,d_dipole);
		dE = sys1.get_d_ene_cross(dwp,atom_to_mv,d_dipole,lambda);
		if(run1.if_accept(dE))
		{
			sys1.update_pos(atom_to_mv,d_dipole);
			run1.accept_scale += 1;
		}
		// update scale

		if(t1 % run1.check_scale == 0)
		{
			run1.update_scale(dwp);
		}

//		cout<<dE<<'\t'<<run1.accept_scale<<endl;
		if(t1 % dump_factor == 1)
		{
			cout<<t1<<'\t'<<sys1.find_dipole()<<'\t'<<run1.scale<<endl;
//			sys1.find_dipole().print();
		}

	}
	output<<sys1<<endl;
	output.close();
	return 0;
}
