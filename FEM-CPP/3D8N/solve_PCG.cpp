#include "solve_PCG.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <algorithm>
using namespace std;

solve_PCG::solve_PCG()
{
}


solve_PCG::~solve_PCG()
{
}


//�õ����ڵ�λ��ֵ������getDisplacement()��
//����	λ������displacement[];
vector<double> solve_PCG::getDisplacement()
{
	displacement.clear();//�����ظ�ʹ��ʱdisplacement�д���֮ǰ������
	displacement.resize(nodeNumber * 3,0.0);//�����ռ�
	/*************************************************************************************************/
	//ѡ��ⷨ��
	if ("EBE" == PCG_choice)
	{
		doPCG_EBE();//�ʺ������ͼ���
	}
	else
	{
		doPCG_Normal();
	}
	/*************************************************************************************************/
	return displacement;
}


//PCGʵ�֣�
void solve_PCG::doPCG_Normal()
{
	/********************************************************************/
	//��ʼ��displacementΪ0����vector�Զ���ʼ��Ϊ0��
	//{R}={loadR}-[K]{displcement}
	vector<double> R = loadR;
	//{Z}=[M](-1){R}
	vector<double> MN = getMN_N();//��[M]����ķ���Ԫ�ؼ��Խ�Ԫ��ȡ��������vectorM�У�
	vector<double> Z = getZ(MN, R);
	//{P}={Z}
	vector<double> P = Z;
	/********************************************************************/
	//iteration begin
	vector<double> Q;
	vector<double> aP((3 * nodeNumber), 0);
	vector<double>::iterator max = max_element(begin(aP), end(aP));//�ҵ����������ֵ�ġ�λ�á���
	do
	{
		/********************************/
		//{Q}=[K]{P},��������նȾ���
		Q = getQ_N(P);
		/********************************/
		//��aϵ��
		double a1 = 0, a2 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			a1 += (Z[i] * R[i]);
			a2 += (P[i] * Q[i]);
		}
		if (abs(a2) < 1E-10)
		{
			cout << "solve_PCG:: doPCG_Normal error!!!" << endl;
			system("pause"); exit(0);
		}
		double a = a1 / a2;
		/********************************/
		//���ڼ���ϵ��b
		double b2 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			b2 += (R[i] * Z[i]);
		}
		if (abs(b2) < 1E-10)
		{
			cout << "solve_PCG:: doPCG_Normal error!!!" << endl;
			system("pause"); exit(0);
		}
		/********************************/
		//{displacement}={displacement}+a*{P}
		//{R}={R}-a*{Q}
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			aP[i] = abs(a*P[i]);
			displacement[i] += a*P[i];
			R[i] = R[i] - a*Q[i];
		}
		/********************************/
		//{Z}=[M](-1){R}
		Z = getZ(MN, R);
		/********************************/
		//��bϵ��
		double b1 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			b1 += (R[i] * Z[i]);
		}
		double b = b1 / b2;
		/********************************/
		//{P}={Z}+b*{P}
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			P[i] = Z[i] + b*P[i];
		}
		/********************************/
		max = max_element(begin(aP), end(aP));
	} while (aP[distance(begin(aP), max)] >(1E-8));
	/********************************************************************/
}


//PCGʵ�֣�
void solve_PCG::doPCG_EBE()
{
	/********************************************************************/
	//��ʼ��displacementΪ0����vector�Զ���ʼ��Ϊ0��
	//{R}={loadR}-[K]{displcement}
	vector<double> R = loadR;
	//{Z}=[M](-1){R}
	vector<double> MN = getMN_EBE();
	vector<double> Z = getZ(MN, R);
	//{P}={Z}
	vector<double> P = Z;
	/********************************************************************/
	//iteration begin
	vector<double> Q;
	vector<double> aP((3 * nodeNumber), 0);
	vector<double>::iterator max = max_element(begin(aP), end(aP));//�ҵ����������ֵ�ġ�λ�á���
	do
	{
		/********************************/
		//{Q}=[K]{P},ʹ��element-by-element method;
		Q = getQ_EBE(P);
		/********************************/
		//��aϵ��
		double a1 = 0, a2 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			a1 += (Z[i] * R[i]);
			a2 += (P[i] * Q[i]);
		}
		if (abs(a2) < 1E-10)
		{
			cout << "solve_PCG:: doPCG_EBE error!!!" << endl;
			system("pause"); exit(0);
		}
		double a = a1 / a2;
		/********************************/
		//���ڼ���ϵ��b
		double b2 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			b2 += (R[i] * Z[i]);
		}
		if (abs(b2) < 1E-10)
		{
			cout << "solve_PCG:: doPCG_EBE error!!!" << endl;
			system("pause"); exit(0);
		}
		/********************************/
		//{displacement}={displacement}+a*{P}
		//{R}={R}-a*{Q}
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			aP[i] = abs(a*P[i]);
			displacement[i] += a*P[i];
			R[i] = R[i] - a*Q[i];
		}
		/********************************/
		//{Z}=[M](-1){R}
		Z = getZ(MN, R);
		/********************************/
		//��bϵ��
		double b1 = 0;
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			b1 += (R[i] * Z[i]);
		}
		double b = b1 / b2;
		/********************************/
		//{P}={Z}+b*{P}
		for (int i = 0; i < (3 * nodeNumber); i++)
		{
			P[i] = Z[i] + b*P[i];
		}
		/********************************/
		max = max_element(begin(aP), end(aP));
	} while (aP[distance(begin(aP), max)] >(1E-8));
	/********************************************************************/
}


//�õ�Ԥ������󣬺���������նȾ���ֵ������նȾ���ֵ�кţ�����նȾ���ֵ�кţ���
//����Ԥ������󣬼����Խ�Ԫ�أ�
vector<double> solve_PCG::getMN_N()
{
	vector<double> MN((3 * nodeNumber), 0);
	for (unsigned int i = 0; i < (globalStiffnessV.size()); i++)
	{
		if (globalStiffnessR[i] == globalStiffnessC[i])
		{
			MN[globalStiffnessR[i]] = 1.0 / globalStiffnessV[i];
		}
	}
	return MN;
}


//�õ�Ԥ������󣬺���������
//����Ԥ������󣬼����Խ�Ԫ�أ�
vector<double> solve_PCG::getMN_EBE()
{
	vector<double> MN((3 * nodeNumber), 0);
	for (int j = 0; j < elementNumber; j++)
	{
		for (int k = 2; k < 10; k++)
		{
			MN[element[j][k] * 3] += objS.elementStiffness(j)[3 * (k - 2)][3 * (k - 2)];
			MN[element[j][k] * 3 + 1] += objS.elementStiffness(j)[3 * (k - 2) + 1][3 * (k - 2) + 1];
			MN[element[j][k] * 3 + 2] += objS.elementStiffness(j)[3 * (k - 2) + 2][3 * (k - 2) + 2];
		}
	}
	//δ���Ǵ��������������������Ĵ���
	for (int i = 0; i < boundNumber; i++)
	{
		MN[int(bound[i][0] * 3 + bound[i][1])] = 1E+30;
	}
	for (int i = 0; i < (3 * nodeNumber); i++)
	{
		MN[i] = 1.0 / MN[i];
	}
	return MN;
}


//���ڼ���Z����
vector<double> solve_PCG::getZ(vector<double> MN, vector<double> R)
{
	vector<double> Z((3 * nodeNumber), 0);
	for (int i = 0; i < (3 * nodeNumber); i++)
	{
		Z[i] = MN[i] * R[i];
	}
	return Z;
}


//���ڼ���Q����
vector<double> solve_PCG::getQ_N(vector<double> P)
{
	vector<double> Q((3 * nodeNumber), 0);
	for (unsigned int i = 0; i < (globalStiffnessV.size()); i++)
	{
		Q[globalStiffnessR[i]] += globalStiffnessV[i] * P[globalStiffnessC[i]];
		if (globalStiffnessR[i] > globalStiffnessC[i])
		{
			Q[globalStiffnessC[i]] += globalStiffnessV[i] * P[globalStiffnessR[i]];//��������������㣬ͬʱע�ⲻҪ�ѶԽ�Ԫ���ظ����㣻
		}
	}
	return Q;
}


//���ڼ���Q����
vector<double> solve_PCG::getQ_EBE(vector<double> P)
{
	vector<double> Q((3 * nodeNumber), 0);
	for (int i = 0; i < elementNumber; i++)
	{
		vector<vector<double>> Pson(24, vector<double>(1, 0));
		for (int j = 2; j < 10; j++)
		{
			Pson[(j - 2) * 3][0] = P[element[i][j] * 3];
			Pson[(j - 2) * 3 + 1][0] = P[element[i][j] * 3 + 1];
			Pson[(j - 2) * 3 + 2][0] = P[element[i][j] * 3 + 2];
		}
		vector<vector<double>> Qson = objM.matrixM(objS.elementStiffness(i), Pson);
		for (int j = 2; j < 10; j++)
		{
			Q[element[i][j] * 3] += Qson[(j - 2) * 3][0];
			Q[element[i][j] * 3 + 1] += Qson[(j - 2) * 3 + 1][0];
			Q[element[i][j] * 3 + 2] += Qson[(j - 2) * 3 + 2][0];
		}
	}
	//���������Ĵ���
	for (int k = 0; k < boundNumber; k++)
	{
		Q[int(bound[k][0] * 3 + bound[k][1])] += P[int(bound[k][0] * 3 + bound[k][1])] * 1E+30;
	}
	return Q;
}


// xiezhuoyu
// mechanics_xzy@163.com