#pragma once
#include <string>
#include <vector>
using namespace std;
/*******************************************
声明：所有矩阵使用vector进行储存。
*******************************************/
extern vector<vector<int>> element;
class often
{
public:
	often();//构造函数；
	~often();//析构函数；

	//矩阵乘法，函数（矩阵A，矩阵B）；
	vector<vector<double>> matrixM(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

	//矩阵转置，函数（矩阵A）；
	vector<vector<double>> matrixT(vector<vector<double>> matrixA);

	//（二阶）矩阵的逆，函数（矩阵A）；
	vector<vector<double>> matrixN2(vector<vector<double>> matrixA);

	//矩阵的逆，函数（矩阵A）；
	vector<vector<double>> matrixN(vector<vector<double>> matrixA);

	//矩阵数乘，函数（矩阵A，数B）；
	vector<vector<double>> matrixM2(vector<vector<double>> matrixA, double B);

	//矩阵相加，函数（矩阵A，矩阵B）；
	vector<vector<double>> matrixP(vector<vector<double>> matrixA, vector<vector<double>>matrixB);

	//矩阵相减，函数（矩阵A，矩阵B）；
	vector<vector<double>> matrixP2(vector<vector<double>>matrixA, vector<vector<double>>matrixB);

	//（二阶）矩阵行列式，函数（矩阵A）；
	double matrixDet2(vector<vector<double>> matrixA);

	//（三阶）矩阵行列式，函数（矩阵A）；
	double matrixDet3(vector<vector<double>> matrixA);

	//矩阵的行列式，函数（矩阵A）；
	double matrixDet(vector<vector<double>> matrixA);

	//局部节点号与全局节点号匹配，函数（单元编号，局部结点对应的行列号，数据文件名）
	//返回										  全局结点对应的行列号
	int match(int elementID, int i, string filename);

	//形函数矩阵N，函数（局部坐标r，局部坐标s）；
	//返回2*8的应变矩阵B；
	vector<vector<double>> getN(double r, double s);

private:
	//矩阵的提取，函数（矩阵A，row，col）；
	vector<vector<double>> matrixSon(vector<vector<double>> matrixA, int row, int col);

};


// xiezhuoyu
// mechanics_xzy@163.com