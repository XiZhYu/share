#include "load.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

load::load()
{
}


load::~load()
{
}


/**************************************************************************************
//处理集中力形成的结点荷载列阵{NLR}，函数nodeloadR()；
将作用在结点上的集中力，组集到整体结点荷载列阵；
***************************************************************************************/
void load::nodeloadR()
{
	for (int k = 0; k < nodeloadNumber; k++)//遍历所有集中荷载；
	{
		vector<double> P{ nodeload[k][1] , nodeload[k][2] ,nodeload[k][3] };//集中荷载大小；
		int ii = 3 * int(nodeload[k][0]);//全局对应局部；
		for (int i = 0; i < 3; i++)
		{
			loadR[ii + i] += P[i];//叠加到整体结点荷载列阵loadRi；
		}
	}
}


//处理面力形成的结点荷载列阵{FLR}，函数faceloadR()；
void load::faceloadR()
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//坐标；
	double H[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };//系数；
	for (int k = 0; k < faceloadNumber; k++)//遍历所有分布荷载；
	{
		vector<vector<double>> FLR(24, vector<double>(1, 0));
		vector<int> elementID_FaceID = faceload_ElementID_FaceID(k);//确定分布力位于几号相关单元以及第几号面；
		int faceID = elementID_FaceID[1];
		vector<vector<double>> GPP = GaussPiontP(k);//高斯点处的荷载值；
		for (int i = 0; i < 3; i++)//高斯积分；
		{
			for (int j = 0; j < 3; j++)
			{
				switch (elementID_FaceID[1])
				{
				case 1://r=+1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(1.0, x[i], x[j]))
								, A(faceID, i, j, elementID_FaceID[0]))
							, (GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				case 2://r=-1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(-1.0, x[i], x[j]))
								, A(faceID, i, j, elementID_FaceID[0])),
								(GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				case 3://s=+1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(x[j], 1.0, x[i]))//注意ij顺序
								, A(faceID, i, j, elementID_FaceID[0])),
								(GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				case 4://s=-1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(x[j], -1.0, x[i]))
								, A(faceID, i, j, elementID_FaceID[0])),
								(GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				case 5://t=+1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(x[i], x[j], 1.0))
								, A(faceID, i, j, elementID_FaceID[0])),
								(GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				case 6://t=-1
					FLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(x[i], x[j], -1.0))
								, A(faceID, i, j, elementID_FaceID[0])),
								(GPP[i][j] * H[i] * H[j]))
						, FLR);
					break;
				default://正确则只可能出现上述情况；
					cout << "load::faceloadR error!!!" << endl;
					system("pause"); exit(0);
					break;
				}
			}
		}
		for (int l = 0; l < 24; l++)
		{
			int ll = objS.match(elementID_FaceID[0], l);//全局对应局部；
			loadR[ll] += FLR[l][0];//叠加到整体等效结点荷载列阵loadRi；
		}
	}
}


/**************************************************************************
//处理体力形成的结点荷载列阵{BLR}，函数bodyloadR()；
应用高斯积分（数值积分）。
高斯点取2*2个，位置为(-0.577350,0.577350)，权系数为(1.000000,1.000000)；
对单元进行遍历，将体力分给单元结点，再组集到整体结点荷载列阵。
***************************************************************************/
void load::bodyloadR()
{
	double x[2] = { -sqrt(1.0 / 3.0), +sqrt(1.0 / 3.0) };//坐标；
	double H[2] = { 1.0,1.0 };//系数；
	for (int q = 0; q < elementNumber; q++)//遍历所有单元；
	{
		vector<vector<double>> BLR(24, vector<double>(1, 0));
		vector<vector<double>> P{ { 0 },{0}, { -material[element[q][1]][2] } };//荷载列阵；
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					//公式：R=[N]t*[P]*|J|*Hi*Hj*Hk并分别求和两次；
					//N(r,s,t),P(x,y),J(r,s,t,elementID)；
					BLR = objM.matrixP(
						objM.matrixM2(
							objM.matrixM(
								objM.matrixT(N(x[i], x[j], x[k]))
								, P)
							, (objM.matrixDet3(objS.J(x[i], x[j], x[k], q))*H[i] * H[j] * H[k]))
						, BLR);
				}
			}
		}
		for (int l = 0; l < 24; l++)
		{
			int ll = objS.match(q, l);//全局对应局部；
			loadR[ll] += BLR[l][0];//叠加到整体等效结点荷载列阵loadRi；
		}
	}
}


/********************************************************************
//组集整体结点荷载列阵loadR，函数getLoadR()；
大数法处理：
 -                        -    -        -       -                -
|E+30…………              |  |已知位移1 |     |已知位移1*（E+30）|		需要考虑体力、集中力、面力、支座反力。
|……                      |  |……      |     |……              |
|…………                  |  |……      |  =  |……              |
|………………E+30          |  |已知位移2 |     |已知位移2*（E+30）|
|…………………………      |  |……      |     |……              |
|………………………………  |  |……      |     |……              |
 -                        -    -        -       -                -
*********************************************************************/
void load::getLoadR()
{
	loadR.clear();
	loadR.resize(3 * nodeNumber);//整体结点荷载列阵大小为：节点数*空间维数；
	nodeloadR();//考虑集中力；
	faceloadR();//考虑面力；
	bodyloadR();//考虑体力；
	//考虑到支座反力的存在，用大数法处理进行处理，将含支座反力这个未知条件不纳入解方程步骤中；
	for (int i = 0; i < boundNumber; i++)
	{
		loadR[int(bound[i][0] * 3 + bound[i][1])] = bound[i][2] * (1E+30);
	}
}


//形函数矩阵N，函数（r，s，t）
vector<vector<double>> load::N(double r = 0, double s = 0, double t = 0)
{
	double	N1 = (1 - r)*(1 - s)*(1 - t) / 8, N2 = (1 + r)*(1 - s)*(1 - t) / 8, N3 = (1 + r)*(1 + s)*(1 - t) / 8, N4 = (1 - r)*(1 + s)*(1 - t) / 8,
		N5 = (1 - r)*(1 - s)*(1 + t) / 8, N6 = (1 + r)*(1 - s)*(1 + t) / 8, N7 = (1 + r)*(1 + s)*(1 + t) / 8, N8 = (1 - r)*(1 + s)*(1 + t) / 8;
	vector<vector<double>> N{
		{ N1,	0,	0,	N2,	0,	0,	N3,	0,	0,	N4,	0,	0,	N5,	0,	0,	N6,	0,	0,	N7,	0,	0,	N8,	0,	0 },
		{ 0,	N1,	0,	0,	N2,	0,	0,	N3,	0,	0,	N4,	0,	0,	N5,	0,	0,	N6,	0,	0,	N7,	0,	0,	N8,	0 },
		{ 0,	0,	N1,	0,	0,	N2,	0,	0,	N3,	0,	0,	N4,	0,	0,	N5,	0,	0,	N6,	0,	0,	N7,	0,	0,	N8 },
	};
	return N;
}


//确定分布力与单元以及面的关系，函数（faceload_ElementID_FaceID）
//返回分布力所在的单元ID和面ID
vector<int> load::faceload_ElementID_FaceID(int faceloadID = 0)
{
	vector<int> result(2, 0);
	for (int i = 0; i < elementNumber; i++)//重复对比单元的结点信息和分布力的结点信息，直至吻合；
	{
		vector<int> f{ int(faceload[faceloadID][0]) ,int(faceload[faceloadID][1]),int(faceload[faceloadID][2]),int(faceload[faceloadID][3]) };
		vector<int> e{ { element[i][2] ,element[i][3],element[i][4], element[i][5], element[i][6] ,element[i][7],element[i][8], element[i][9] } };
		vector<bool> panju{
				(f[0] == e[1]) && (f[1] == e[2]) && (f[2] == e[6]) && (f[3] == e[5])
			,	(f[0] == e[0]) && (f[1] == e[3]) && (f[2] == e[7]) && (f[3] == e[4])
			,	(f[0] == e[3]) && (f[1] == e[7]) && (f[2] == e[6]) && (f[3] == e[2])
			,	(f[0] == e[0]) && (f[1] == e[4]) && (f[2] == e[5]) && (f[3] == e[1])
			,	(f[0] == e[4]) && (f[1] == e[5]) && (f[2] == e[6]) && (f[3] == e[7])
			,	(f[0] == e[0]) && (f[1] == e[1]) && (f[2] == e[2]) && (f[3] == e[3]) };
		if (panju[0]) { result[0] = i; result[1] = 1; break; }
		if (panju[1]) { result[0] = i; result[1] = 2; break; }
		if (panju[2]) { result[0] = i; result[1] = 3; break; }
		if (panju[3]) { result[0] = i; result[1] = 4; break; }
		if (panju[4]) { result[0] = i; result[1] = 5; break; }
		if (panju[5]) { result[0] = i; result[1] = 6; break; }
	}
	if (0 == result[1])
	{
		cout << "load::faceload_ElementID_FaceID no result!!!" << endl;
		system("pause"); exit(0);
	}
	return result;
}


//计算3*3情况下的高斯点处的荷载值
vector<vector<double>> load::GaussPiontP(int faceloadID = 0)
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//坐标；
	vector<vector<double>> GPP(3, vector<double>(3, 0.0));
	for (int i = 0; i < 3; i++)//利用形函数的特性
	{
		for (int j = 0; j < 3; j++)
		{
			GPP[i][j] = faceload[faceloadID][4] * (1 - x[i])*(1 - x[j]) / 4 + faceload[faceloadID][5] * (1 + x[i])*(1 - x[j]) / 4
				+ faceload[faceloadID][6] * (1 + x[i])*(1 + x[j]) / 4 + faceload[faceloadID][7] * (1 - x[i])*(1 + x[j]) / 4;
		}
	}
	return GPP;
}


//计算公式中的列阵项，函数（）
//返回{A}
vector<vector<double>> load::A(int faceID = 1, int ii = 0, int jj = 0, int elementID = 0)
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//坐标；
	vector<double> xr(3, 0.0);
	vector<double> xs(3, 0.0);
	vector<double> xt(3, 0.0);
	vector<vector<double>> result(3, vector<double>(1, 0.0));
	switch (faceID)
	{
	case 1://r=+1
		for (int i = 0; i < 3; i++)
		{
			xs[i] = (	(1 - x[jj])*(node[element[elementID][1 + 3]][i] - node[element[elementID][1 + 2]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 6]][i])) / 4;
			xt[i] = (	(1 - x[ii])*(node[element[elementID][1 + 6]][i] - node[element[elementID][1 + 2]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 3]][i])) / 4;
		}
		result = {
					{ xs[1] * xt[2] - xs[2] * xt[1] }
			,		{ xs[2] * xt[0] - xs[0] * xt[2] }
			,		{ xs[0] * xt[1] - xs[1] * xt[0] }
		};
		break;
	case 2://r=-1
		for (int i = 0; i < 3; i++)
		{
			xs[i] = (	(1 - x[jj])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 4]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 5]][i] - node[element[elementID][1 + 8]][i])) / 4;
			xt[i] = (	(1 - x[ii])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 5]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 4]][i] - node[element[elementID][1 + 8]][i])) / 4;
		}
		result = {
					{ xs[1] * xt[2] - xs[2] * xt[1] }
			,		{ xs[2] * xt[0] - xs[0] * xt[2] }
			,		{ xs[0] * xt[1] - xs[1] * xt[0] }
		};
		break;
	case 3://s=+1
		for (int i = 0; i < 3; i++)
		{
			xr[i] = (	(1 - x[ii])*(node[element[elementID][1 + 3]][i] - node[element[elementID][1 + 4]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 8]][i])) / 4;
			xt[i] = (	(1 - x[jj])*(node[element[elementID][1 + 8]][i] - node[element[elementID][1 + 4]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 3]][i])) / 4;
		}
		result = {
					{ xt[1] * xr[2] - xt[2] * xr[1] }
			,		{ xt[2] * xr[0] - xt[0] * xr[2] }
			,		{ xt[0] * xr[1] - xt[1] * xr[0] }
		};
		break;
	case 4://s=-1
		for (int i = 0; i < 3; i++)
		{
			xr[i] = (	(1 - x[ii])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 2]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 5]][i] - node[element[elementID][1 + 6]][i])) / 4;
			xt[i] = (	(1 - x[jj])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 5]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 2]][i] - node[element[elementID][1 + 6]][i])) / 4;
		}
		result = {
					{ xt[1] * xr[2] - xt[2] * xr[1] }
			,		{ xt[2] * xr[0] - xt[0] * xr[2] }
			,		{ xt[0] * xr[1] - xt[1] * xr[0] }
		};
		break;
	case 5://t=+1
		for (int i = 0; i < 3; i++)
		{
			xr[i] = (	(1 - x[jj])*(node[element[elementID][1 + 6]][i] - node[element[elementID][1 + 5]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 8]][i])) / 4;
			xs[i] = (	(1 - x[ii])*(node[element[elementID][1 + 8]][i] - node[element[elementID][1 + 5]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 7]][i] - node[element[elementID][1 + 6]][i])) / 4;
		}
		result = {
					{ xr[1] * xs[2] - xr[2] * xs[1] }
					,{ xr[2] * xs[0] - xr[0] * xs[2] }
					,{ xr[0] * xs[1] - xr[1] * xs[0] }
		};
		break;
	case 6://t=-1
		for (int i = 0; i < 3; i++)
		{
			xr[i] = (	(1 - x[jj])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 2]][i])
				+		(1 + x[jj])*(node[element[elementID][1 + 4]][i] - node[element[elementID][1 + 3]][i])) / 4;
			xs[i] = (	(1 - x[ii])*(node[element[elementID][1 + 1]][i] - node[element[elementID][1 + 4]][i])
				+		(1 + x[ii])*(node[element[elementID][1 + 2]][i] - node[element[elementID][1 + 3]][i])) / 4;
		}
		result = {
					{ xr[1] * xs[2] - xr[2] * xs[1] }
					,{ xr[2] * xs[0] - xr[0] * xs[2] }
					,{ xr[0] * xs[1] - xr[1] * xs[0] }
		};
		break;
	default://正确则只可能出现上述情况；
		cout << "load::A error!!!" << endl;
		system("pause"); exit(0);
		break;
	}
	return result;
}

// xiezhuoyu
// mechanics_xzy@163.com