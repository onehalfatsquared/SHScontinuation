#pragma once
#include "dlib/optimization.h"

using namespace dlib;
typedef matrix<double,0,1> column_vector;

#define DIMENSION 3 //set to 3d optimization

//general stuff
bool isDoF(int i);
int delta(int i, int j);

//morse stuff
double morseP(double r, double rho, double E);
double morsePD(double r, double rho, double E);
double morseEval(double* particles, double rho, double E, int N);
void morseGrad(double* particles, double rho, double E, int N, column_vector& g);
void hessMorse(column_vector cluster, double rho, double E, int N, matrix<double>& H);

//lj stuff
double ljP(double r, double rho, double E);
double ljPD(double r, double rho, double E);
double ljEval(double* particles, double rho, double E, int N);
void ljGrad(double* particles, double rho, double E, int N, column_vector& g);
void hessLJ(double* particles, double rho, double E, int N, double* H);

//sticky stuff
void stickyF(double E, double rho, double beta, double k0, double& f, double& fprime);
double stickyNewton(double E, double rho, double k0, double beta);

//particle stuff
void c2p(column_vector cluster, double* particles, int N);
double euDist(double* particles, int i, int j, int N, double* Z);