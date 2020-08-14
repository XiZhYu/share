#include "stiffness.h"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

stiffness::stiffness()
{
}


stiffness::~stiffness()
{
}


//形函数对局部坐标的偏导数，函数（局部坐标r，局部坐标s，局部坐标t）
vector<vector<vector<double>>> stiffness::Nri(double r = 0.0, double s = 0.0, double t = 0.0)
{
	vector<vector<vector<double>>> Nri(8, vector<vector<double>>(3, vector<double>(1, 0)));
	Nri[0][0][0] = -(1 - s)*(1 - t) / 8, Nri[1][0][0] = (1 - s)*(1 - t) / 8, Nri[2][0][0] = (1 + s)*(1 - t) / 8, Nri[3][0][0] = -(1 + s)*(1 - t) / 8, Nri[4][0][0] = -(1 - s)*(1 + t) / 8, Nri[5][0][0] = (1 - s)*(1 + t) / 8, Nri[6][0][0] = (1 + s)*(1 + t) / 8, Nri[7][0][0] = -(1 + s)*(1 + t) / 8;
	Nri[0][1][0] = -(1 - r)*(1 - t) / 8, Nri[1][1][0] = -(1 + r)*(1 - t) / 8, Nri[2][1][0] = (1 + r)*(1 - t) / 8, Nri[3][1][0] = (1 - r)*(1 - t) / 8, Nri[4][1][0] = -(1 - r)*(1 + t) / 8, Nri[5][1][0] = -(1 + r)*(1 + t) / 8, Nri[6][1][0] = (1 + r)*(1 + t) / 8, Nri[7][1][0] = (1 - r)*(1 + t) / 8;
	Nri[0][2][0] = -(1 - r)*(1 - s) / 8, Nri[1][2][0] = -(1 + r)*(1 - s) / 8, Nri[2][2][0] = -(1 + r)*(1 + s) / 8, Nri[3][2][0] = -(1 - r)*(1 + s) / 8, Nri[4][2][0] = (1 - r)*(1 - s) / 8, Nri[5][2][0] = (1 + r)*(1 - s) / 8, Nri[6][2][0] = (1 + r)*(1 + s) / 8, Nri[7][2][0] = (1 - r)*(1 + s) / 8;
	return Nri;
}


//雅各比矩阵J，函数（局部坐标r，局部坐标s，局部坐标t，单元编号ID）；
vector<vector<double>> stiffness::J(double r = 0.0, double s = 0.0, double t = 0.0, int elementID = 0)
{
	vector<vector<double>> P{
		{-(1 - s)*(1 - t) / 8,(1 - s)*(1 - t) / 8,(1 + s)*(1 - t) / 8,-(1 + s)*(1 - t) / 8,-(1 - s)*(1 + t) / 8,(1 - s)*(1 + t) / 8,(1 + s)*(1 + t) / 8,-(1 + s)*(1 + t) / 8 },
		{-(1 - r)*(1 - t) / 8,-(1 + r)*(1 - t) / 8,(1 + r)*(1 - t) / 8,(1 - r)*(1 - t) / 8,-(1 - r)*(1 + t) / 8,-(1 + r)*(1 + t) / 8,(1 + r)*(1 + t) / 8,(1 - r)*(1 + t) / 8 },
		{-(1 - r)*(1 - s) / 8,-(1 + r)*(1 - s) / 8,-(1 + r)*(1 + s) / 8,-(1 - r)*(1 + s) / 8,(1 - r)*(1 - s) / 8,(1 + r)*(1 - s) / 8,(1 + r)*(1 + s) / 8,(1 - r)*(1 + s) / 8 }
	};
	vector<vector<double>> C{
		{ node[element[elementID][2]][0], node[element[elementID][2]][1], node[element[elementID][2]][2] },
		{ node[element[elementID][3]][0], node[element[elementID][3]][1], node[element[elementID][3]][2] },
		{ node[element[elementID][4]][0], node[element[elementID][4]][1], node[element[elementID][4]][2] },
		{ node[element[elementID][5]][0], node[element[elementID][5]][1], node[element[elementID][5]][2] },
		{ node[element[elementID][6]][0], node[element[elementID][6]][1], node[element[elementID][6]][2] },
		{ node[element[elementID][7]][0], node[element[elementID][7]][1], node[element[elementID][7]][2] },
		{ node[element[elementID][8]][0], node[element[elementID][8]][1], node[element[elementID][8]][2] },
		{ node[element[elementID][9]][0], node[element[elementID][9]][1], node[element[elementID][9]][2] },
	};
	vector<vector<double>> J = objM.matrixM(P, C);
	return J;
}


//应变矩阵B，函数（局部坐标r，局部坐标s，局部坐标t，单元编号ID）；
vector<vector<double>> stiffness::B(double r = 0.0, double s = 0.0, double t = 0.0, int elementID = 0)
{
	vector<vector<vector<double>>> Nxi(8,vector<vector<double>>(3,vector<double>(1,0)));
	for (int i = 0; i < 8; i++)
	{
		Nxi[i] = objM.matrixM(objM.matrixN(J(r, s, t, elementID)), Nri(r, s, t)[i]);
	}
	vector<vector<double>> B{
		{ Nxi[0][0][0],0,0,	Nxi[1][0][0],0,0,	Nxi[2][0][0],0,0,	Nxi[3][0][0],0,0,	Nxi[4][0][0],0,0,	Nxi[5][0][0],0,0,	Nxi[6][0][0],0,0,	Nxi[7][0][0],0,0 },
		{ 0,Nxi[0][1][0],0,	0,Nxi[1][1][0],0,	0,Nxi[2][1][0],0,	0,Nxi[3][1][0],0,	0,Nxi[4][1][0],0,	0,Nxi[5][1][0],0,	0,Nxi[6][1][0],0,	0,Nxi[7][1][0],0 },
		{ 0,0,Nxi[0][2][0],	0,0,Nxi[1][2][0],	0,0,Nxi[2][2][0],	0,0,Nxi[3][2][0],	0,0,Nxi[4][2][0],	0,0,Nxi[5][2][0],	0,0,Nxi[6][2][0],	0,0,Nxi[7][2][0] },
		{ Nxi[0][1][0],Nxi[0][0][0],0,	Nxi[1][1][0],Nxi[1][0][0],0,	Nxi[2][1][0],Nxi[2][0][0],0,	Nxi[3][1][0],Nxi[3][0][0],0,	Nxi[4][1][0],Nxi[4][0][0],0,	Nxi[5][1][0],Nxi[5][0][0],0,	Nxi[6][1][0],Nxi[6][0][0],0,	Nxi[7][1][0],Nxi[7][0][0],0 },
		{ 0,Nxi[0][2][0],Nxi[0][1][0],	0,Nxi[1][2][0],Nxi[1][1][0],	0,Nxi[2][2][0],Nxi[2][1][0],	0,Nxi[3][2][0],Nxi[3][1][0],	0,Nxi[4][2][0],Nxi[4][1][0],	0,Nxi[5][2][0],Nxi[5][1][0],	0,Nxi[6][2][0],Nxi[6][1][0],	0,Nxi[7][2][0],Nxi[7][1][0] },
		{ Nxi[0][2][0],0,Nxi[0][0][0],	Nxi[1][2][0],0,Nxi[1][0][0],	Nxi[2][2][0],0,Nxi[2][0][0],	Nxi[3][2][0],0,Nxi[3][0][0],	Nxi[4][2][0],0,Nxi[4][0][0],	Nxi[5][2][0],0,Nxi[5][0][0],	Nxi[6][2][0],0,Nxi[6][0][0],	Nxi[7][2][0],0,Nxi[7][0][0] }
	};
	return B;
}


//弹性矩阵D，函数（单元编号ID）；
vector<vector<double>> stiffness::D(int elementID = 0)
{
	double E = material[element[elementID][1]][0];
	double mu = material[element[elementID][1]][1];
	double gamma = material[element[elementID][1]][2];
	double lambda = E*mu / ((1 + mu)*(1 - 2 * mu));
	double G = E / (2 * (1 + mu));
	vector<vector<double>> D{
		{lambda+2*G,	lambda,		lambda,		0,	0,	0},
		{lambda,		lambda+2*G,	lambda,		0,	0,	0},
		{lambda,		lambda,		lambda+2*G,	0,	0,	0},
		{0,				0,			0,			G,	0,	0},
		{0,				0,			0,			0,	G,	0},
		{0,				0,			0,			0,	0,	G}
	};
	return D;
}


/********************************************************
//生成单元刚度矩阵elementK;
应用高斯积分（数值积分）。
高斯点取2*2个为(-0.577，+0.577)，权系数为(1.0，1.0)。
*********************************************************/
vector<vector<double>>  stiffness::elementStiffness(int elementID = 0)
{
	double H[2] = { 1.0,1.0 };
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };
	vector<vector<double>>  sonElementStiffness(24, vector<double>(24, 0.0));
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				sonElementStiffness = 
					objM.matrixP(
						objM.matrixM2(
							objM.matrixM2(
								objM.matrixM(
									objM.matrixM(
										objM.matrixT(B(x[i], x[j], x[k], elementID))
									, D(elementID))
								, B(x[i], x[j], x[k], elementID))
							, objM.matrixDet3(J(x[i], x[j], x[k], elementID)))
						, (H[i] * H[j] * H[k]))
					, sonElementStiffness);
			}
		}
	}
	return sonElementStiffness;
}


/*****************************************************************
//组集整体劲度矩阵，函数（）；
循环遍历所有单元；
根据单元局部结点号（单元劲度矩阵行列号）找到全局结点号（整体劲度矩阵行列号）；
将元素值和对应行列号存入COO数组中。
******************************************************************/
/**********************************************************************************
大数法处理：
 -                        -    -        -       -                -
|E+30…………              |  |已知位移1 |     |已知位移1*（E+30）|
|……                      |  |……      |     |……              |
|…………                  |  |……      |  =  |……              |
|………………E+30          |  |已知位移2 |     |已知位移2*（E+30）|
|…………………………      |  |……      |     |……              |
|………………………………  |  |……      |     |……              |
 -                        -    -        -       -                -
**********************************************************************************/
void stiffness::globalStiffness()
{
	for (int k = 0; k < elementNumber; k++)
	{
		vector<vector<double>>  value = elementStiffness(k);
		for (int j = 0; j < 24; j++)
		{
			for (int i = j; i < 24; i++)//只需要下三角元素，包括对角线元素；
			{
				int ii = match(k, i);//全局对应局部；
				int jj = match(k, j);
				if (ii < jj)
				{
					int temp = ii;
					ii = jj;
					jj = temp;
				}//将整体刚度矩阵的上三角元素换到下三角去；
				if ((abs(value[i][j]) > 1E-30))
				{
					COO(value[i][j], ii, jj);//只储存非零元素；
				}
			}
		}
	}
	//大数法处理；
	for (int i = 0; i < boundNumber; i++)
	{
		for (unsigned int j = 0; j < globalStiffnessR.size(); j++)
		{
			if ((globalStiffnessR.at(j) == int(bound[i][0] * 3 + bound[i][1])) && (globalStiffnessR.at(j) == globalStiffnessC.at(j)))
			{
				globalStiffnessV.at(j) = 1E+30;
			}
		}
	}
}


/*************************************************************
//局部结点序号与全局结点序号匹配,函数（单元号，局部结点序号）；
//返回全局的结点序号；
给出整体劲度矩阵的两个结点，则对应有四个元素:
										j结点
									   3j  3j+1	3j+2  列号
						i结点	3i     xx   xy	xz
								3i+1   yx   yy	yz
								3i+2   zx	zy	zz
						行号
**************************************************************/
int stiffness::match(int elementID = 0, int i = 0)
{
	int ii;
	if (0 == (i % 3))//取余；
	{
		ii = 3 * element[elementID][i / 3 + 2];//取商；
	}
	else if (1 == (i % 3))
	{
		ii = 3 * element[elementID][i / 3 + 2] + 1;
	}
	else
	{
		ii = 3 * element[elementID][i / 3 + 2] + 2;
	}
	return ii;
}


/*************************************************************
//COO存储，函数（元素值，行号，列号）；
存储稀疏矩阵。
存储内容：1.矩阵非0元素值globalStiffnessV（value）
2.矩阵非0元素的行号globalStiffnessR（rowNumber）
3.矩阵非0元素的列号globalStiffnessC（columnNumber）
且只存下三角非零元素。
**************************************************************/
void stiffness::COO(double value = 0, int i = 0, int j = 0)
{
	static int VRCi = 0;//用于记录非0元素个数；
	bool panju = false;
	for (int l = 0; l < VRCi; l++)
	{
		if ((i == globalStiffnessR.at(l)) && (j == globalStiffnessC.at(l)))
		{
			globalStiffnessV.at(l) += value;
			panju = true;
		}
	}
	if (false == panju)
	{
		globalStiffnessV.push_back(value);
		globalStiffnessR.push_back(i);
		globalStiffnessC.push_back(j);
		VRCi++;
	}
}


// xiezhuoyu
// mechanics_xzy@163.com