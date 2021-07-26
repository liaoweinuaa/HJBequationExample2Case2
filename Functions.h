#pragma once

#include "Constant.h"
#include <string>
using namespace std;
void Func_Init();

void Func_RecursionMT();

void Func_SavaData(string filename0, string filename1);

bool Func_IsTargetSet(State s);

State Func_Transition(State s, double u);
