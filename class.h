#include <iostream>
#include <cstring>
#include <fstream>
#ifndef MY_CLASS
#define MY_CLASS

class vec;
class pot;
class atom;
class cell;
class mc;

class vec
{
private:
	double x[3];
friend atom;
friend cell;
public:
	void import(double *);
	void clean();			// reset value to zero
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec operator/(const double&);
	vec & operator=(const vec&);
	vec & operator=(double*);	// can replace import
	double operator*(const vec&);
	double norm();
	friend std::ostream& operator<<(std::ostream&,vec);
	//debug
	void print();
};

class pot
{
private:
	double ene, min;
	double aa, bb;
friend cell;
friend mc;
public:
	void init(double E, double L);	//E, L
	double get_E(double r);
};

class atom
{
private:
	vec pos, pos0, dipole;
friend cell;
public:
	void get_pos0(vec &);
	void get_pos(vec &);
	void get_dipole();
	atom & operator=(const atom&);
};

class cell
{
private:
	double lat;
	int numx, num;	//along one direction, tot number
	atom * a_l;
public:
//	double ene_onsite0, ene_onsite1;
//	double ene_dipole0, ene_dipole1;
//	double ene_tot0, ene_tot1;

	void init(double length, int NN, pot& dwp);		//unit lattice, number of atoms on each axis
	void update_pos(int n0, vec& d_new_pos);
	double get_d_ene(pot& dwp, int n0, vec& d_new_pos);
	vec find_dipole();
};

class mc
{
public:
    double scale;       // scale random move vector
    int check_scale;    // number of runs to adjust scale
    double accept_scale;    // update scale
	double T;				// temperature of system
public:
    void init(double Scale, int Check_scale, double Temperature);
    void update_scale(pot& dwp);
    void mv_atm(int atom_num, int &target, vec & d_dipole);
	int if_accept(double dE);
};

#endif
