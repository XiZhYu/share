#pragma once
#include <string>
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/***********************************************************************************
计算集中力、分布力、体力对应的结点荷载列阵{nodeLoadR}、{faceLoadR}、{bodyLoadR}
并组装成整体结点荷载列阵{loadR}
************************************************************************************/
extern vector<vector<double>> node;
extern int nodeNumber;
extern vector<vector<int>> element;
extern int elementNumber;
extern vector<vector<double>> material;
extern vector<vector<double>> bound;
extern int boundNumber;
extern vector<vector<double>> faceload;
extern int faceloadNumber;
extern vector<vector<double>> nodeload;
extern int nodeloadNumber;

//整体结点荷载列阵
extern vector<double> loadR;

class load
{
public:
	load();
	~load();

	//组集整体结点荷载列阵loadR，函数doLoadR()；
	void getLoadR();
	
	//形函数矩阵N，函数(r,s,t)
	vector<vector<double>> N(double r, double s, double t);

	//处理集中力形成的结点荷载列阵{NLR}，函数nodeloadR()；
	void nodeloadR();

	//处理面力形成的结点荷载列阵{FLR}，函数doFaceload()；
	void faceloadR();

	//处理体力形成的结点荷载列阵{BLR}，函数doBodyforce()；
	void bodyloadR();

	//确定分布力与单元以及面的关系，函数（faceloadID）
	//返回分布力所在的单元ID和面ID
	vector<int> faceload_ElementID_FaceID(int faceloadID);

	//计算3*3情况下的高斯点处的荷载值
	vector<vector<double>> GaussPiontP(int faceloadID);

	//计算公式中的列阵项，函数（）
	//返回{A}
	vector<vector<double>> A(int faceID, int i, int j, int elementID);

private:
	stiffness objS;
	matrix objM;
};

// xiezhuoyu
// mechanics_xzy@163.com