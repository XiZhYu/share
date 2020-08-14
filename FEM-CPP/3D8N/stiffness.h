#pragma once
#include <string>
#include <vector>
#include "matrix.h"
using namespace std;
/**********************************************************
计算单元刚度矩阵，并组集成整体刚度矩阵。
***********************************************************/
extern int elementNumber;
extern int boundNumber;
extern vector<vector<double>> bound;
extern vector<vector<int>> element;
extern vector<vector<double>> node;
extern vector<vector<double>> material;

//整体劲度矩阵非0元素值（上三角）；
extern vector <double> globalStiffnessV;
//整体劲度矩阵非0元素的行号；
extern vector <int> globalStiffnessR;
//整体劲度矩阵非0元素的列号；
extern vector <int> globalStiffnessC;

class stiffness
{
public:
	stiffness();
	~stiffness();

	//形函数对局部坐标的偏导数，函数（局部坐标r，局部坐标s，局部坐标t，局部坐标t）
	vector<vector<vector<double>>> Nri(double r, double s, double t);

	//雅各比矩阵J，函数（局部坐标r，局部坐标s，局部坐标t，单元编号ID）；
	//返回3*3的雅各比矩阵J；
	vector<vector<double>> J(double r, double s, double t,  int elementID);

	//应变矩阵B，函数（局部坐标r，局部坐标s，局部坐标t，单元编号ID）；
	//返回6*24的应变矩阵B；
	vector<vector<double>> B(double r, double s, double t, int elementID);

	//弹性矩阵D，函数（单元编号ID）；
	//返回6*6的弹性矩阵D；
	vector<vector<double>> D(int elementID);

	//计算单元刚度矩阵；
	//返回24*24的单元劲度矩阵；
	vector<vector<double>> elementStiffness(int elementID);

	//组集整体劲度矩阵，函数（）；
	//以COO格式存储整体劲度矩阵；
	void globalStiffness();

	//局部结点序号与全局结点序号匹配,函数（单元号，局部结点序号）；
	//返回全局的结点序号；
	int match(int elementID, int i);

	//COO存储，函数（元素值，行号，列号）；
	void COO(double value, int i, int j);

private:
	matrix objM;
};


// xiezhuoyu
// mechanics_xzy@163.com