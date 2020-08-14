#include "matrix.h"
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

matrix::matrix()
{
}


matrix::~matrix()
{
}


//矩阵乘法，函数（矩阵A，矩阵B）；
vector<vector<double>> matrix::matrixM(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	unsigned int ARow = matrixA.size(), ACol = matrixA[0].size();//vector.size()类型为unsigned int；
	unsigned int BRow = matrixB.size(), BCol = matrixB[0].size();
	if (ACol != BRow)
	{
		cout << "matrix::matrixM input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> matrixC(ARow, vector<double>(BCol, 0.0));//将所有元素初始化为0；
	for (unsigned int i = 0; i < ARow; i++)//相乘算法；
	{
		for (unsigned int j = 0; j < BCol; j++)
		{
			double sum = 0;
			for (unsigned int k = 0; k < ACol; k++)
			{
				//对应二维数组应为：C[i,j]+=A[i,k]*B[k,j]；
				sum += matrixA[i][k] * matrixB[k][j];
			}
			matrixC[i][j] = sum;
		}
	}
	return matrixC;
}


//矩阵转置，函数（矩阵A）
vector<vector<double>> matrix::matrixT(vector<vector<double>> matrixA)
{
	vector<vector<double>> matrixAT(matrixA[0].size(), vector<double>(matrixA.size(), 0));
	for (unsigned int i = 0; i < matrixA.size(); i++)
	{
		for (unsigned int j = 0; j < matrixA[0].size(); j++)
		{
			matrixAT[j][i] = matrixA[i][j];
		}
	}
	return matrixAT;
}


//（二阶）矩阵的逆，函数（矩阵A）；
vector<vector<double>> matrix::matrixN2(vector<vector<double>> matrixA)
{
	if ((matrixA.size() != 2) || (matrixA[0].size() != 2))
	{
		cout << "matrix::matrixN2 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double detA = matrixA[0][0] * matrixA[1][1] - matrixA[0][1] * matrixA[1][0];
	vector<vector<double>> matrixAN{ { (matrixA[1][1] / detA),(-matrixA[0][1] / detA) },
	{ (-matrixA[1][0] / detA),(matrixA[0][0] / detA) } };
	return matrixAN;
}


//矩阵数乘，函数（矩阵A，数B）；
vector<vector<double>> matrix::matrixM2(vector<vector<double>> matrixA, double B = 0)
{
	vector<vector<double>> matrixC(matrixA.size(), vector<double>(matrixA[0].size(), 0));
	for (unsigned int i = 0; i < matrixA.size(); i++)
	{
		for (unsigned int j = 0; j < matrixA[0].size(); j++)
		{
			matrixC[i][j] = matrixA[i][j] * B;
		}
	}
	return matrixC;
}


//矩阵相加，函数（矩阵A，矩阵B）；
vector<vector<double>> matrix::matrixP(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size()))
	{
		cout << "matrix::matrixP input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> matrixC(matrixA.size(), vector<double>(matrixA[0].size(), 0));
	for (unsigned int i = 0; i < matrixA.size(); i++)
	{
		for (unsigned int j = 0; j < matrixA[0].size(); j++)
		{
			matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}
	return matrixC;
}


//矩阵相减，函数（矩阵A，矩阵B）；
vector<vector<double>> matrix::matrixP2(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size()))
	{
		cout << "matrix::matrixP2 input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> matrixC(matrixA.size(), vector<double>(matrixA[0].size(), 0));
	for (unsigned int i = 0; i < matrixA.size(); i++)
	{
		for (unsigned int j = 0; j < matrixA[0].size(); j++)
		{
			matrixC[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}
	return matrixC;
}


//（二阶）矩阵行列式，函数（矩阵A）；
double matrix::matrixDet2(vector<vector<double>> matrixA)
{
	if ((matrixA.size() != 2) || (matrixA[0].size() != 2))
	{
		cout << "matrix::matrixDet2 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = matrixA[0][0] * matrixA[1][1] - matrixA[0][1] * matrixA[1][0];
	return Det;
}


//矩阵的行列式，函数（矩阵A）；
//尽量不用递归，速度慢！！！！！！！！！！！！！；
double matrix::matrixDet(vector<vector<double>> matrixA)
{
	if (matrixA.size() != matrixA[0].size())
	{
		cout << "matrix::matrixDet input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = 0;
	if (1 == matrixA.size())//递归；
	{
		Det = matrixA[0][0];
	}
	else
	{
		for (unsigned int i = 0; i < matrixA.size(); i++)//按第一行展开求行列式；
		{
			vector<vector<double>> ASon = matrixSon(matrixA, 0, i);//余子式对应的子矩阵；
			Det += pow(double(-1), int(i))*matrixA[0][i] * matrixDet(ASon);
		}
	}
	return Det;
}


/*****************************************************************
n*n阶矩阵的提取，函数（矩阵A，row，col）；
将矩阵的第（row+1）行、第（col+1）列剔除,余下子矩阵存入son中。
行列数和行列号仔细区分，如row行对应第（row+1）行。
主要思路：将即不在row行、也不在col列的元素按原顺序存入son中。
主要用于求矩阵行列式和矩阵的逆。
******************************************************************/
vector<vector<double>> matrix::matrixSon(vector<vector<double>> matrixA, int row = 0, int col = 0)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != ACol) || (1 == ARow) || (row > (ARow - 1)))
	{
		cout << "matrix::matrixSon input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> son((ARow - 1), vector<double>((ACol - 1), 0));
	for (unsigned int i = 0; i < ARow; i++)
	{
		for (unsigned int j = 0; j < ACol; j++)
		{
			if ((i < row) && (j < col))
			{
				son[i][j] = matrixA[i][j];
			}
			else if ((i < row) && (j > col))
			{
				son[i][j - 1] = matrixA[i][j];
			}
			else if ((i > row) && (j > col))
			{
				son[i - 1][j - 1] = matrixA[i][j];
			}
			else if ((i > row) && (j < col))
			{
				son[i - 1][j] = matrixA[i][j];
			}
			else
			{
			}
		}
	}
	return son;
}


/**************************************************************
矩阵的逆，函数（矩阵A）；
求出矩阵的伴随矩阵和矩阵的行列式，即可得到矩阵的逆。
***************************************************************/
vector<vector<double>> matrix::matrixN(vector<vector<double>> matrixA)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != ACol) || (0 == matrixDet(matrixA)))
	{
		cout << "matrix::matrixN input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> AN01(ARow, vector<double>(ACol, 0));
	for (unsigned int i = 0; i < ARow; i++)
	{
		for (unsigned int j = 0; j < ACol; j++)
		{
			vector<vector<double>> son = matrixSon(matrixA, i, j);
			AN01[i][j] = (pow(double(-1), (i + j))
				*matrixDet(son))
				/ matrixDet(matrixA);
		}
	}
	vector<vector<double>>AN = matrixT(AN01);
	return AN;
}


//（三阶）矩阵行列式，函数（矩阵A）；
double matrix::matrixDet3(vector<vector<double>> matrixA)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != 3) || (ACol != 3))
	{
		cout << "matrix::matrixDet3 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = matrixA[0][0] * (matrixA[1][1] * matrixA[2][2] - matrixA[2][1] * matrixA[1][2])
		- matrixA[0][1] * (matrixA[1][0] * matrixA[2][2] - matrixA[2][0] * matrixA[1][2])
		+ matrixA[0][2] * (matrixA[1][0] * matrixA[2][1] - matrixA[2][0] * matrixA[1][1]);//按第一行展开；
	return Det;
}


// xiezhuoyu
// mechanics_xzy@163.com