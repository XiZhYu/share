#pragma once
#include <string>
#include <vector>
using namespace std;
/*******************************************
所有矩阵使用vector进行储存。
*******************************************/
class matrix
{
public:
	matrix();
	~matrix();

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

	//矩阵的提取，函数（矩阵A，row，col）；
	vector<vector<double>> matrixSon(vector<vector<double>> matrixA, int row, int col);

private:

};


// xiezhuoyu
// mechanics_xzy@163.com