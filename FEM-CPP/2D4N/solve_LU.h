#pragma once
#include <string>
#include <vector>
using namespace std;
/*********************************************************
	用于求解支配方程[K]*{f}={R}；
	利用LU分解求解线性代数方程组；
	LU分解步骤主要为：分解得到L、U矩阵，前代，回代。
	最终得到的结果为结点位移值。
**********************************************************/
extern int nodeNumber;

class solve_LU
{
public:
	solve_LU(string file);
	~solve_LU(void);

	//得到各节点位移值，函数getDisplacement()；
	//返回	位移列阵displacement[];
	vector<double> getDisplacement();

private:
	string filename;

	//分解得到L、U矩阵中的L矩阵，函数LUglobalK()；
	//将计算结果存入LU[]矩阵中；
	void doLU();

	//前代求解[L][g]=[R]中的[g]，函数forward（）；
	void forward();

	//回代求解[U][displacement]=[g]，函数backwards（）；
	void backwards();

	vector<double> LU;//用于存储LU分解矩阵，只存下三角的L矩阵；

	//用于操作LU，函数getLUij(行号，列号)；
	int getLUij(int row, int col);

	vector<double> displacement;

	//输出位移值到FEMout.txt中
	void out();
};


// xiezhuoyu
// mechanics_xzy@163.com