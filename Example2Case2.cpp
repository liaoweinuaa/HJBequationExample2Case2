
#include <iostream>
#include "Constant.h"
#include "Functions.h"
#include <string>

using namespace std;

double(*Arr_StateValue)[Ny][Ntheta] = new double[Nx][Ny][Ntheta];
double(*Arr_StateValueNew)[Ny][Ntheta] = new double[Nx][Ny][Ntheta];
State(*Arr_State)[Ny][Ntheta] = new State[Nx][Ny][Ntheta];



int main()
{
	Func_Init();
	int filenameindex = 0;
	for (int i = 0; i < 205; i++)
	{
		Func_RecursionMT();
		cout << "The " << i + 1 << "-th recursion is completed" << endl;
	}
	Func_SavaData("solution.dat", "solution.csv");
}

