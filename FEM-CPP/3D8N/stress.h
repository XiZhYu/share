#pragma once
#include <string>
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/*********************************************************
成果整理。
得到单元中心和结点处的应力值；
[P]=[D][B][f]
**********************************************************/
extern int elementNumber;
extern int nodeNumber;
extern string solve_choice;

extern vector<double> displacement;

//各单元中心处应力值：x、y、z方向正应力、xy、yz、zx方向切应力；
extern vector<vector<double>> centerStress;
//各结点应力值：x、y、z方向正应力、xy、yz、zx方向切应力；
extern vector<vector<double>> nodeStress;

class stress
{
public:
	stress();
	~stress();

	//计算应力值，函数doStress（）；
	void doStress();

	//得到各单元中心处应力值，getCenterStress();
	vector<vector<double>> getCenterStress();

	//得到各结点应力值，getNodeStress();
	vector<vector<double>> getNodeStress();

	//将结点对应的单元存入数组，函数nodeLinkElement()
	//nodeLinkElement[i][j]为i号结点的第j个相关单元的ID
	void nodeLinkElement();

private:
	stiffness objS;
	matrix objM;

	vector<vector<int>> v_nodeLinkElement;
};


// xiezhuoyu
// mechanics_xzy@163.com