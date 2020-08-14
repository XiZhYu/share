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


//���㵥Ԫ���ĵ㴦��Ӧ��ֵ������doStress������
//���ɵ�Ԫ���ĵ㣨Ӧ���ѵ㣩��Ӧ��ֵ�������Ƶ�ƽ�����õ�����Ӧ��ֵ
void stress::doStress()
{
	centerStress.clear();
	centerStress.resize(elementNumber, vector<double>(6, 0.0));
	for (int i = 0; i < elementNumber; i++)
	{
		//���ɵ�Ԫ�Ľ��λ������;
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


//������Ӧ�ĵ�Ԫ�������飬����nodeLinkElement()
//nodeLinkElement[i][j]Ϊi�Ž��ĵ�j����ص�Ԫ��ID
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


//�õ�����Ԫ���Ĵ�Ӧ��ֵ��getCenterStress();
vector<vector<double>> stress::getCenterStress()
{
	return centerStress;
}

//�õ������Ӧ��ֵ��getNodeStress();
vector<vector<double>> stress::getNodeStress()
{
	return nodeStress;
}


// xiezhuoyu
// mechanics_xzy@163.com