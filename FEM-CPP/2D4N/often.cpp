#include <iostream>
#include <stdlib.h>
#include <math.h>//pow(),�˷����㣻
#include <vector>
#include <string>
#include "often.h"
using namespace std;

often::often()
{
}


//����˷�������������A������B����
vector<vector<double>> often::matrixM(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	unsigned int ARow = matrixA.size(), ACol = matrixA[0].size();//vector.size()����Ϊunsigned int��
	unsigned int BRow = matrixB.size(), BCol = matrixB[0].size();
	if (ACol != BRow)
	{
		cout << "often::matrixM input error!!!" << endl;
		system("pause"); exit(0);
	}
	vector<vector<double>> matrixC(ARow, vector<double>(BCol, 0.0));//������Ԫ�س�ʼ��Ϊ0��
	for (unsigned int i = 0; i < ARow; i++)//����㷨��
	{
		for (unsigned int j = 0; j < BCol; j++)
		{
			double sum = 0;
			for (unsigned int k = 0; k < ACol; k++)
			{
				//��Ӧ��ά����ӦΪ��C[i,j]+=A[i,k]*B[k,j]��
				sum += matrixA[i][k] * matrixB[k][j];
			}
			matrixC[i][j] = sum;
		}
	}
	return matrixC;
}


//����ת�ã�����������A��
vector<vector<double>> often::matrixT(vector<vector<double>> matrixA)
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


//�����ף�������棬����������A����
vector<vector<double>> often::matrixN2(vector<vector<double>> matrixA)
{
	if ((matrixA.size() != 2) || (matrixA[0].size() != 2))
	{
		cout << "often::matrixN2 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double detA = matrixA[0][0] * matrixA[1][1] - matrixA[0][1] * matrixA[1][0];
	vector<vector<double>> matrixAN{ { (matrixA[1][1] / detA),(-matrixA[0][1] / detA) },
										{(-matrixA[1][0] / detA),(matrixA[0][0] / detA)} };
	return matrixAN;
}


//�������ˣ�����������A����B����
vector<vector<double>> often::matrixM2(vector<vector<double>> matrixA, double B = 0)
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


//������ӣ�����������A������B����
vector<vector<double>> often::matrixP(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size()))
	{
		cout << "often::matrixP input error!!!" << endl;
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


//�������������������A������B����
vector<vector<double>> often::matrixP2(vector<vector<double>> matrixA, vector<vector<double>> matrixB)
{
	if ((matrixA.size() != matrixB.size()) || (matrixA[0].size() != matrixB[0].size()))
	{
		cout << "often::matrixP2 input error!!!" << endl;
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


//�����ף���������ʽ������������A����
double often::matrixDet2(vector<vector<double>> matrixA)
{
	if ((matrixA.size() != 2) || (matrixA[0].size() != 2))
	{
		cout << "often::matrixDet2 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = matrixA[0][0] * matrixA[1][1] - matrixA[0][1] * matrixA[1][0];
	return Det;
}


//���������ʽ������������A����
//�������õݹ飬�ٶ�������������������������������
double often::matrixDet(vector<vector<double>> matrixA)
{
	if (matrixA.size() != matrixA[0].size())
	{
		cout << "often::matrixDet input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = 0;
	if (1 == matrixA.size())//�ݹ飻
	{
		Det = matrixA[0][0];
	}
	else
	{
		for (unsigned int i = 0; i < matrixA.size(); i++)//����һ��չ��������ʽ��
		{
			vector<vector<double>> ASon = matrixSon(matrixA, 0, i);//����ʽ��Ӧ���Ӿ���
			Det += pow(double(-1), int(i))*matrixA[0][i] * matrixDet(ASon);
		}
	}
	return Det;
}


/*****************************************************************
n*n�׾������ȡ������������A��row��col����
	������ĵڣ�row+1���С��ڣ�col+1�����޳�,�����Ӿ������son�С�
	�����������к���ϸ���֣���row�ж�Ӧ�ڣ�row+1���С�
	��Ҫ˼·����������row�С�Ҳ����col�е�Ԫ�ذ�ԭ˳�����son�С�
	��Ҫ�������������ʽ�;�����档
******************************************************************/
vector<vector<double>> often::matrixSon(vector<vector<double>> matrixA, int row = 0, int col = 0)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != ACol) || (1 == ARow) || (row > (ARow - 1)))
	{
		cout << "often::matrixSon input error!!!" << endl;
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
������棬����������A����
	�������İ������;��������ʽ�����ɵõ�������档
***************************************************************/
vector<vector<double>> often::matrixN(vector<vector<double>> matrixA)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != ACol) || (0 == matrixDet(matrixA)))
	{
		cout << "often::matrixN input error!!!" << endl;
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


//�����ף���������ʽ������������A����
double often::matrixDet3(vector<vector<double>> matrixA)
{
	unsigned int ARow = matrixA.size();
	unsigned int ACol = matrixA[0].size();
	if ((ARow != 3) || (ACol != 3))
	{
		cout << "often::matrixDet3 input error!!!" << endl;
		system("pause"); exit(0);
	}
	double Det = matrixA[0][0] * (matrixA[1][1] * matrixA[2][2] - matrixA[2][1] * matrixA[1][2])
		- matrixA[0][1] * (matrixA[1][0] * matrixA[2][2] - matrixA[2][0] * matrixA[1][2])
		+ matrixA[0][2] * (matrixA[1][0] * matrixA[2][1] - matrixA[2][0] * matrixA[1][1]);//����һ��չ����
	return Det;
}


//�κ�������N���������ֲ�����r���ֲ�����s����
vector<vector<double>> often::getN(double r = 0, double s = 0)
{
	double N1 = 0.25*(1 - r)*(1 - s), N2 = 0.25*(1 + r)*(1 - s), N3 = 0.25*(1 + r)*(1 + s), N4 = 0.25*(1 - r)*(1 + s);
	vector<vector<double>> N{ {N1,0,N2,0,N3,0,N4,0} ,
								{0,N1,0,N2,0,N3,0,N4} };
	return N;
}


//�ֲ��ڵ����ȫ�ֽڵ��ƥ�䣻
/*******************************************************************************
�������徢�Ⱦ����������㣬���Ӧ���ĸ�Ԫ��:						j���
																   2j  2j+1  �к�
													i���	2i     xx   xy
															2i+1   yx   yy
															�к�
*******************************************************************************/
int often::match(int elementID = 0, int i = 0, string filename = "02.txt")
{
	int ii;
	if (0 == (i % 2))
	{
		ii = 2 * element[elementID][i / 2 + 2];
	}
	else
	{
		ii = 2 * element[elementID][i / 2 + 2] + 1;
	}
	return ii;
}


often::~often()
{
}


// xiezhuoyu
// mechanics_xzy@163.com