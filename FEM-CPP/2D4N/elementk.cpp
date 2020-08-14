#include <iostream>
#include <vector>
#include <math.h>//sqrt()；
#include "often.h"
#include "elementk.h"

using namespace std;


elementK::elementK(string file)
{
	filename = file;
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
void elementK::doGlobalK()
{
	often runOft;
	for (int k = 0; k < elementNumber; k++)
	{
		vector<vector<double>>  value = getElementK(k);
		for (int j = 0; j < 8; j++)
		{
			for (int i = j; i < 8; i++)//只需要下三角元素；
			{
				int ii = runOft.match(k, i, filename);//全局对应局部；
				int jj = runOft.match(k, j, filename);
				if (ii < jj)
				{
					int t = ii;
					ii = jj;
					jj = t;
				}//将整体刚度矩阵的上三角元素换到下三角去；
				if ((abs(value[i][j]) > 1E-30))
				{
					useCOO(value[i][j], ii, jj);//只储存非零元素；
				}
			}
		}
	}
	//大数法处理；
	for (int i = 0; i < boundNumber; i++)
	{
		for (unsigned int j = 0; j < globalKR.size(); j++)
		{
			if ((globalKR.at(j) == int(bound[i][1] * 2 + bound[i][2])) && (globalKR.at(j) == globalKC.at(j)))
			{
				globalKV.at(j) = 1E+30;
			}
		}
	}
}


/*************************************************************
//COO存储，函数（元素值，行号，列号）；
	存储稀疏矩阵。
	存储内容：1.矩阵非0元素值globalKV（value）
			  2.矩阵非0元素的行号globalKR（rowNumber）
			  3.矩阵非0元素的列号globalKC（columnNumber）
	且只存下三角非零元素。
**************************************************************/
void elementK::useCOO(double value = 0, int i = 0, int j = 0)
{
	static int VRCi = 0;//用于记录非0元素个数；
	bool panju = false;
	for (int l = 0; l < VRCi; l++)
	{
		if ((i == globalKR.at(l)) && (j == globalKC.at(l)))
		{
			globalKV.at(l) += value;
			panju = true;
		}
	}
	if (false == panju)
	{
		globalKV.push_back(value);
		globalKR.push_back(i);
		globalKC.push_back(j);
		VRCi++;
	}
}


/********************************************************
//生成单元刚度矩阵elementK;
	应用高斯积分（数值积分）。
	高斯点取3*3个为(-0.774596,0.0,+0.774596)，权系数为(0.555555,0.888888,0.555555)；
	改为取2*2个,满足精度要求。
*********************************************************/
vector<vector<double>>  elementK::getElementK(int elementID = 0)
{
	double H[2] = { 1.0,1.0 };
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };
	vector<vector<double>>  Ke(8, vector<double>(8, 0));
	often run2;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Ke = run2.matrixP(
				run2.matrixM2(
					run2.matrixM2(
						run2.matrixM(
							run2.matrixM(
								run2.matrixT(getStrB(x[i], x[j], elementID))
								, getElaD(elementID))
							, getStrB(x[i], x[j], elementID))
						, run2.matrixDet2(getJacJ(x[i], x[j], elementID)))
					, (H[i] * H[j]))
				, Ke);
		}
	}
	return Ke;
}


//雅各比矩阵J，函数（局部坐标r，局部坐标s，单元编号ID）；
vector<vector<double>> elementK::getJacJ(double r = 0.0, double s = 0.0, int elementID = 0)
{

	vector<vector<double>> P{ { (-(1 - s) / 4) ,((1 - s) / 4),((1 + s) / 4),(-(1 + s) / 4)},
								{(-(1 - r) / 4),(-(1 + r) / 4),((1 + r) / 4),((1 - r) / 4)} };
	vector<vector<double>> C{ { node[element[elementID][2]][1] ,node[element[elementID][2]][2] },
								{ node[element[elementID][3]][1] ,node[element[elementID][3]][2] },
								{ node[element[elementID][4]][1] ,node[element[elementID][4]][2] },
								{ node[element[elementID][5]][1] ,node[element[elementID][5]][2] } };
	often run2;
	vector<vector<double>> J = run2.matrixM(P, C);
	return J;
}


//应变矩阵B，函数（局部坐标r，局部坐标s，单元编号ID）；
vector<vector<double>> elementK::getStrB(double r = 0.0, double s = 0.0, int elementID = 0)
{
	double N1x, N2x, N3x, N4x, N1y, N2y, N3y, N4y;
	vector<vector<double>> P{ {-(1 - s) / 4,(1 - s) / 4,(1 + s) / 4,-(1 + s) / 4},
								{-(1 - r) / 4,-(1 + r) / 4,(1 + r) / 4,(1 - r) / 4} };
	often run2;
	vector<vector<double>> newdel01 = getJacJ(r, s, elementID);
	vector<vector<double>> JN = run2.matrixN2(newdel01);
	N1x = JN[0][0] * P[0][0] + JN[0][1] * P[1][0]; N2x = JN[0][0] * P[0][1] + JN[0][1] * P[1][1];
	N3x = JN[0][0] * P[0][2] + JN[0][1] * P[1][2]; N4x = JN[0][0] * P[0][3] + JN[0][1] * P[1][3];
	N1y = JN[1][0] * P[0][0] + JN[1][1] * P[1][0]; N2y = JN[1][0] * P[0][1] + JN[1][1] * P[1][1];
	N3y = JN[1][0] * P[0][2] + JN[1][1] * P[1][2]; N4y = JN[1][0] * P[0][3] + JN[1][1] * P[1][3];
	vector<vector<double>> B{ {N1x,0,N2x,0,N3x,0,N4x,0},
								{0,N1y,0,N2y,0,N3y,0,N4y},
								{N1y,N1x,N2y,N2x,N3y,N3x,N4y,N4x} };
	return B;
}


//（平面应力问题）弹性矩阵D，函数（单元编号ID）；
vector<vector<double>> elementK::getElaD(int elementID = 0)
{
	double E = material[element[elementID][1]][1];
	double u = material[element[elementID][1]][2];
	double y = material[element[elementID][1]][3];
	vector<vector<double>> D{ {1,u,0} ,{u,1,0}, {0,0,((1 - u) / 2) } };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			D[i][j] = D[i][j] * (E / (1 - u*u));
		}
	}
	return D;
}

elementK::~elementK()
{
}


// xiezhuoyu
// mechanics_xzy@163.com