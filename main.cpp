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
	pot dwp;
	cell sys1;
	mc run1;
	double dw_a, dw_b, dw_c, length = 1.0, max_length_nei = 1.1, temperature = 0.1, l_xx = 0.1, l_xy = 0.1;
	int num_x = 10,num_tot;
	int tot_run = 100000, dump_factor = 100;
	int check_scale = 100;
	vec d_dipole;
	int atom_to_mv;
	double dE;

//	double dd[3]={0.0,0.0,-0.1};

	// get potential, cell parameters
	cin>>dw_a>>dw_b>>dw_c>>l_xx>>l_xy>>length>>max_length_nei>>num_x>>temperature;
	cin>>tot_run>>dump_factor>>check_scale;
	num_tot = num_x*num_x*num_x;

	// initiate potential, cell, and mc
	dwp.init_2d(dw_a, dw_b, dw_c);
	sys1.init_2d(length,num_x,dwp);
	run1.init(length/5,check_scale,temperature);
	sys1.get_neighbor(max_length_nei);

	// start mc
	for(int t1=1; t1<tot_run; t1++)
	{
		run1.mv_atm_2d(sys1,atom_to_mv,d_dipole);
		dE = sys1.get_d_ene_short_2d(dwp,atom_to_mv,d_dipole,l_xx, l_xy);
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
	return 0;
}
