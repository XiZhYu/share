#pragma once
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/*******************************************************
预处理共轭梯度法precondition conjugate gradient。
根据Jacobian Precondition，预处理矩阵取主元素对角阵。
********************************************************/
extern int nodeNumber;
extern int elementNumber;
extern vector<vector<int>> element;
extern int boundNumber;
extern vector<vector<double>> bound;
extern string PCG_choice;

extern vector <double> globalStiffnessV;
extern vector <int> globalStiffnessR;
extern vector <int> globalStiffnessC;
extern vector<double> loadR;

//整体结点位移列阵
extern vector<double> displacement;

class solve_PCG
{
public:
	solve_PCG();
	~solve_PCG();

	//得到各节点位移值，函数getDisplacement()；
	//返回	位移列阵displacement[];
	vector<double> getDisplacement();

	//PCG实现；
	void doPCG_Normal();
	void doPCG_EBE();

	//得到预处理矩阵，函数（整体刚度矩阵值，整体刚度矩阵值行号，整体刚度矩阵值列号）；
	//返回预处理矩阵，即主对角元素；
	vector<double> getMN_N();
	vector<double> getMN_EBE();

	//用于计算Z列阵；
	vector<double> getZ(vector<double> M, vector<double> R);

	//用于计算Q列阵；
	vector<double> getQ_N(vector<double> P);
	vector<double> getQ_EBE(vector<double> P);

private:
	matrix objM;
	stiffness objS;
};


// xiezhuoyu
// mechanics_xzy@163.com