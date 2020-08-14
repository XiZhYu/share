#pragma once
#include <vector>
using namespace std;
/*********************************************************
用于求解支配方程[K]*{f}={R}；
利用LU分解求解线性代数方程组；
LU分解步骤主要为：分解得到L、U矩阵，前代，回代。
最终得到的结果为结点位移值。
**********************************************************/
extern int nodeNumber;

extern vector <double> globalStiffnessV;
extern vector <int> globalStiffnessR;
extern vector <int> globalStiffnessC;
extern vector<double> loadR;

//整体结点位移列阵
extern vector<double> displacement;

class solve_LU
{
public:
	solve_LU();
	~solve_LU(void);

	//得到各节点位移值，函数getDisplacement()；
	//返回	位移列阵displacement[];
	vector<double> getDisplacement();

	//分解得到L、U矩阵中的L矩阵，函数doLU()；
	//将计算结果存入LU[]矩阵中；
	void doLU();

	//前代求解[L][g]=[R]中的[g]，函数forward（）；
	void forward();

	//回代求解[U][displacement]=[g]，函数backwards（）；
	void backwards();

	//用于操作LU，函数getLUij(行号，列号)；
	int getLUij(int row, int col);

private:
	vector<double> LU;//用于存储LU分解矩阵，只存下三角的L矩阵；
	
};


// xiezhuoyu
// mechanics_xzy@163.com