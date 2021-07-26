#include "Constant.h"
#include <iostream>
#include <string>
#include <mutex>
#include <thread>
#include "Functions.h"
#include <fstream>
using namespace std;

//Arraies to hold the values of the solution at the grid points
extern double(*Arr_StateValue)[Ny][Ntheta];
extern double(*Arr_StateValueNew)[Ny][Ntheta];

//The states represented by the grid points
extern State(*Arr_State)[Ny][Ntheta];


//Endpoint cost function
double Func_EndPointCost(State s0)
{
	double x = s0.x, y = s0.y, theta = s0.theta;
	if (theta > PI)
	{
		theta = 2 * PI - theta;
	}

	return -exp(-x * x - y * y - theta);
}

void Func_Init()
{
	

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Ntheta; k++)
			{
				Arr_State[i][j][k].x = xmin + gridx * i;
				Arr_State[i][j][k].y = ymin + gridy * j;
				Arr_State[i][j][k].theta = thetamin + gridtheta * k;
				//Terminal condition of HJB equation
				Arr_StateValue[i][j][k] = Func_EndPointCost(Arr_State[i][j][k]);
			}
		}
	}
}

//Target set
bool Func_IsTargetSet(State s)
{
	if (s.x >= -0.5 && s.x <= 0.5 && s.y >= 1.5 && s.y <= 2.5)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//Dynamic system, s(t+dt)=s(t)+f(s,u)dt+ 0.5 * df/ds * f(s,u)*dt*dt
State Func_Transition(State s, double u)
{
	if (Func_IsTargetSet(s))
	{
		return s;
	}

	double x = s.x, y = s.y, theta = s.theta;

	double dotx = v * cos(theta) + 1 * y + 0.1 * y * y * y;
	double doty = v * sin(theta) - 1 * x - 0.1 * x * x * x;
	double dottheta = u;

	double ax = (1 + 0.3 * y * y) * doty - v * sin(theta) * dottheta;
	double ay = (-1 - 0.3 * x * x) * dotx + v * cos(theta) * dottheta;
	double atheta = 0;

	double newx = x + dotx * dt + 0.5 * ax * dt * dt;
	double newy = y + doty * dt + 0.5 * ay * dt * dt;
	double newtheta = theta + dottheta * dt + 0.5 * atheta * dt * dt;


	if (newx >= xmax)
	{
		newx = xmax - 0.0001;
	}

	if (newx <= xmin)
	{
		newx = xmin + 0.0001;
	}

	if (newy >= ymax)
	{
		newy = ymax - 0.0001;
	}

	if (newy <= ymin)
	{
		newy = ymin + 0.0001;
	}

	if (newtheta > 2 * PI)
	{
		newtheta = newtheta - 2 * PI;
	}

	if (newtheta < 0)
	{
		newtheta = newtheta + 2 * PI;
	}

	return State{
		newx,newy,newtheta
	};
}



//Trilinear interpolation see https://en.wikipedia.org/wiki/Trilinear_interpolation
double Func_InterPolation(State s0)
{
	if (Func_IsTargetSet(s0))
	{
		return Func_EndPointCost(s0);
	}
	double ialpha_double = ((s0.x - xmin) / gridx);
	double itheta_double = ((s0.y - ymin) / gridy);
	double iq_double = ((s0.theta - thetamin) / gridtheta);

	int ialpha = int(ialpha_double);
	int itheta = int(itheta_double);
	int iq = int(iq_double);

	double v000 = Arr_StateValue[ialpha][itheta][iq];
	double v001 = Arr_StateValue[ialpha][itheta][iq + 1];
	double v010 = Arr_StateValue[ialpha][itheta + 1][iq];
	double v011 = Arr_StateValue[ialpha][itheta + 1][iq + 1];
	double v100 = Arr_StateValue[ialpha + 1][itheta][iq];
	double v101 = Arr_StateValue[ialpha + 1][itheta][iq + 1];
	double v110 = Arr_StateValue[ialpha + 1][itheta + 1][iq];
	double v111 = Arr_StateValue[ialpha + 1][itheta + 1][iq + 1];

	//if (v000==0 && v001 == 0 && v010==0 && v011 == 0 && v100 == 0 && v101 == 0 && v110 == 0 && v111 == 0)
	//{
	//	return 0;
	//}

	double dalpha0 = ialpha_double - ialpha;
	double dalpha1 = 1 - dalpha0;
	double dtheta0 = itheta_double - itheta;
	double dtheta1 = 1 - dtheta0;
	double dq0 = iq_double - iq;
	double dq1 = 1 - dq0;

	double V000 = dalpha0 * dtheta0 * dq0;
	double V001 = dalpha0 * dtheta0 * dq1;
	double V010 = dalpha0 * dtheta1 * dq0;
	double V011 = dalpha0 * dtheta1 * dq1;
	double V100 = dalpha1 * dtheta0 * dq0;
	double V101 = dalpha1 * dtheta0 * dq1;
	double V110 = dalpha1 * dtheta1 * dq0;
	double V111 = dalpha1 * dtheta1 * dq1;


	return (v000 * V111 + v001 * V110 + v010 * V101 + v011 * V100 + v100 * V011 + v101 * V010 + v110 * V001 + v111 * V000);
}

//Running cost function, c((s0+s1)/2,u), where s1 is the state transferred from s0 with one time step 
double Func_RunningCost(State s0, State s1, double u)
{
	if (Func_IsTargetSet(s0))
	{
		return 0;
	}

	double x = 0.5 * (s0.x + s1.x);
	double y = 0.5 * (s0.y + s1.y);
	double theta = 0.5 * (s0.theta + s1.theta);

	double dotx = v * cos(theta) + 1 * y + 0.1 * y * y * y;
	double doty = v * sin(theta) - 1 * x - 0.1 * x * x * x;

	return 1 + 0.1 * sqrt(dotx * dotx + doty * doty);
}

//Bellman's principle of optimality
double Func_ValueUnderOptInput(State s0)
{
	double minvalue = 1000000;
	for (int i = 0; i < Nu; i++)
	{

		double u = umin + i * gridu;

		State snew = Func_Transition(s0, u);
		double snewvalue = Func_InterPolation(snew);
		double value = Func_RunningCost(s0, snew,u) * dt + snewvalue;
		if (value < minvalue)
		{
			minvalue = value;
		}

	}
	return minvalue;
}


//Single-threaded recursion
void Func_RecursionST(int threadid)
{
	int index = threadid;
	while (true)
	{
		if (index >= Ntotal)
		{
			break;
		}

		int ix = index / (Ny * Ntheta);
		int iy = (index - ix * (Ny * Ntheta)) / Ntheta;
		int iz = index - ix * (Ny * Ntheta) - iy * Ntheta;

		Arr_StateValueNew[ix][iy][iz] = Func_ValueUnderOptInput(Arr_State[ix][iy][iz]);
		index += threadnum;
	}
}

//Ten-thread recursion
void Func_RecursionMT()
{

	thread t0(Func_RecursionST, 0);
	thread t1(Func_RecursionST, 1);
	thread t2(Func_RecursionST, 2);
	thread t3(Func_RecursionST, 3);
	thread t4(Func_RecursionST, 4);
	thread t5(Func_RecursionST, 5);
	thread t6(Func_RecursionST, 6);
	thread t7(Func_RecursionST, 7);
	thread t8(Func_RecursionST, 8);
	thread t9(Func_RecursionST, 9);


	t0.join();
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	t5.join();
	t6.join();
	t7.join();
	t8.join();
	t9.join();

	memcpy(Arr_StateValue, Arr_StateValueNew, sizeof(double) * Nx * Ny * Ntheta);
}

void Func_SavaData(string filename0, string filename1)
{
	//Save data in binary form
	ofstream ofs(filename0, ios::binary | ios::out);
	ofs.write((const char*)Arr_StateValue, sizeof(double) * Nx * Ny * Ntheta);
	ofs.close();

	//Save data in text form
	ofstream ofs1;
	ofs1.open(filename1);
	for (int i = 0; i < Nx; i += savestep)
	{
		for (int j = 0; j < Ny; j += savestep)
		{
			for (int k = 0; k < Ntheta; k += savestep)
			{
				ofs1 << Arr_StateValue[i][j][k] << endl;
			}
		}
	}
	ofs1.close();
}