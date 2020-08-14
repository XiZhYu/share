#include "stress.h"
#include <iostream>
#include <vector>
#include "stress.h"
using namespace std;

stress::stress()
{
}


stress::~stress()
{
}


//计算单元中心点处的应力值，函数doStress（）；
//并由单元中心点（应力佳点）的应力值，运用绕点平均法得到结点的应力值
void stress::doStress()
{
	centerStress.clear();
	centerStress.resize(elementNumber, vector<double>(6, 0.0));
	for (int i = 0; i < elementNumber; i++)
	{
		//集成单元的结点位移列阵;
		vector<vector<double>> f(24, vector<double>(1, 0.0));
		for (int j = 0; j < 24; j++)
		{
			f[j][0] = displacement[objS.match(i, j)];
		}
		vector<vector<double>> DBf = objM.matrixM(
										objM.matrixM(
											objS.D(i)
											, objS.B(0, 0, 0, i))
										, f);
		for (int j = 0; j < 6; j++)
		{
			centerStress[i][j] = DBf[j][0];
		}
	}
	nodeLinkElement();
	nodeStress.clear();
	nodeStress.resize(nodeNumber , vector<double>(6, 0.0));
	for (int i = 0; i < nodeNumber; i++)
	{
		vector<vector<double>> mean(1, vector<double>(6, 0.0));
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < v_nodeLinkElement[i].size(); k++)
			{
				mean[0][j] += centerStress[v_nodeLinkElement[i][k]][j] ;
			}
		}
		mean = objM.matrixM2(mean, (1.0 / v_nodeLinkElement[i].size()));
		nodeStress[i] = mean[0];
	}
}


//将结点对应的单元存入数组，函数nodeLinkElement()
//nodeLinkElement[i][j]为i号结点的第j个相关单元的ID
void stress::nodeLinkElement()
{
	vector<vector<int>> result(nodeNumber,vector<int>(0,0));
	for (int i = 0; i < elementNumber; i++)
	{
		for (int j = 2; j < 10; j++)
		{
			result[element[i][j]].push_back(i);
		}
	}
	v_nodeLinkElement = result;
}


//得到各单元中心处应力值，getCenterStress();
vector<vector<double>> stress::getCenterStress()
{
	return centerStress;
}

//得到各结点应力值，getNodeStress();
vector<vector<double>> stress::getNodeStress()
{
	return nodeStress;
}


// xiezhuoyu
// mechanics_xzy@163.com