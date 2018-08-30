#include <iostream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "class.h"

using namespace std;

int main()
{
	srand(111);
	pot dvp;
	cell sys1;
	double barrier = 1.0, x0 = 0.1, length = 1.0;
	int num_x = 10;

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
	ddd = dd;
	temp = sys1.get_d_ene(dvp, 1, ddd);
	cout<<setw(13)<<setprecision(8)<<fixed<<temp<<endl;
	return 0;
}
