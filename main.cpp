#include <iostream>
#include "class.h"

using namespace std;

int main()
{
	pot dvp;
	cell sys1;

	dvp.init(1,1);
	sys1.init(10,10);
	sys1.get_ene(dvp);
	cout<<sys1.ene_onsite1<<endl;
	return 0;
}
