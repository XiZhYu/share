#pragma once
#include <string>
#include <vector>
using namespace std;
/**********************************************************
	计算整体刚度矩阵；
***********************************************************/
extern int elementNumber;
extern int boundNumber;
extern vector<vector<double>> bound;
extern vector<vector<int>> element;
extern vector<vector<double>> node;
extern vector<vector<double>> material;
class elementK
{
public:
	elementK(string file);//构造函数；
	~elementK();

	//计算单元刚度矩阵；
	//返回8*8的单元劲度矩阵；
	vector<vector<double>> getElementK(int elementID);

	//组集整体劲度矩阵，函数（）；
	//以COO格式存储整体劲度矩阵；
	void doGlobalK();//createGlobalK；

	//应变矩阵B，函数（局部坐标r，局部坐标s，单元编号ID）；
	//返回3*8的应变矩阵B；
	vector<vector<double>> getStrB(double r, double s, int elementID);//getStrainB

	//（平面应力问题）弹性矩阵D，函数（单元编号ID）；
	//返回3*3的弹性矩阵D；
	vector<vector<double>> getElaD(int elementID);//getElasticityD

	//雅各比矩阵J，函数（局部坐标r，局部坐标s，单元编号ID）；
	//返回2*2的雅各比矩阵J；
	vector<vector<double>> getJacJ(double r, double s, int elementID);

	//整体劲度矩阵非0元素值（上三角）；
	vector <double> globalKV;
	//整体劲度矩阵非0元素的行号；
	vector <int> globalKR;
	//整体劲度矩阵非0元素的列号；
	vector <int> globalKC;

private:
	string filename;//数据文件名；

	//COO存储，函数（元素值，行号，列号）；
	void useCOO(double value, int i, int j);
};


// xiezhuoyu
// mechanics_xzy@163.com