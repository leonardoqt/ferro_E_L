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
friend mc;
friend pot;
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
	friend std::ofstream& operator<<(std::ofstream&,vec);
	//debug
	void print();
};

class pot
{
private:
	double ene, min;
	double x2,y2,xy,x4,y4,x2y2;
friend cell;
friend mc;
public:
	void init(double E, double L);
	void init(double a0,double a1,double a2,double a3, double a4, double a5);
	double get_E(vec & p);
};

class atom
{
public:
	vec pos, pos0, dipole;
friend cell;
friend mc;
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
	int ** nei;		// 1d index of neighbor of each site
	int  * num_nei;	// number of neighbors
friend mc;
public:
	void init(double length, int NN, pot& dwp);		//unit lattice, number of atoms on each axis
	void get_neighbor(double max_d);
	void update_pos(int n0, vec& d_pos);
	double get_d_ene_cross(pot& dwp, int n0, vec& d_pos, double lambda);
	vec find_dipole();
	friend std::ofstream& operator<<(std::ofstream&,cell&);
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
    void mv_atm(cell &sys1, int &target, vec &d_dipole);
	int if_accept(double dE);
};

#endif
