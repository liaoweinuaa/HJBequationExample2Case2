#pragma once

const double PI = 3.141592654;

//Computational domain
const double xmin = -4;
const double xmax = 4;
const double ymin = -4;
const double ymax = 4;
const double thetamin = 0;
const double thetamax = 2 * PI;
const int threadnum = 10;

//Number of grid points
const int Nx = 257;
const int Ny = 257;
const int Ntheta = 257;

const int Ntotal = Nx * Ny * Ntheta;
const double gridx = (xmax - xmin) / (Nx - 1);
const double gridy = (ymax - ymin) / (Ny - 1);
const double gridtheta = (thetamax - thetamin) / (Ntheta - 1);

//Admissible control inputs
const double umin = -0.5;
const double umax = 0.5;

const int Nu = 21;
const double gridu = (umax - umin) / (Nu - 1);
const double v = 1;

//Time step length
const double dt = 0.02;

//The grid resolution used for saving data
const int savestep = 2;

struct State
{
	double x;
	double y;
	double theta;
};
