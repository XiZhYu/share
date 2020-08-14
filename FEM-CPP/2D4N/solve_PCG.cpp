#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "solve_PCG.h"
#include "elementk.h"
#include "load.h"
#include "often.h"
using namespace std;

solve_PCG::solve_PCG(string file)
{
	filename = file;
}


//�õ����ڵ�λ��ֵ������getDisplacement()��
//����	λ������displacement[];
vector<double> solve_PCG::getDisplacement()
{
	displacement.resize(nodeNumber * 2);//�����ռ䲢��ʼ��Ϊ0��
	/*************************************************************************************************/
	//ѡ��ⷨ��
	if ("N" == PCG_choice)
	{
		doPCG_Normal();
	}
	else if ("EBE" == PCG_choice)//�������ش��������
	{
		doPCG_EBE();
	}
	else
	{
		doPCG_Normal();
	}
	/*************************************************************************************************/
	out();
	return displacement;
}

//PCGʵ�֣�
void solve_PCG::doPCG_Normal()
{
	load runloa(filename);		runloa.doLoadR();
	elementK runele(filename);	runele.doGlobalK();
	/********************************************************************/
	//��ʼ��displacementΪ0����vector�Զ���ʼ��Ϊ0��
	//{R}={loadR}-[K]{displcement}
	vector<double> R = runloa.loadR;
	//{Z}=[M](-1){R}
	vector<double> MN = getMN(runele.globalKV, runele.globalKR, runele.globalKC);//��[M]����ķ���Ԫ�ؼ��Խ�Ԫ��ȡ��������vectorM�У�
	vector<double> Z = getZ(MN, R);
	//{P}={Z}
	vector<double> P = Z;
	/********************************************************************/
	//iteration begin
	vector<double> Q;
	vector<double> aP((2 * nodeNumber), 0);
	vector<double>::iterator max = max_element(begin(aP), end(aP));//�ҵ����������ֵ�ġ�λ�á���
	do
	{
		/********************************/
		//{Q}=[K]{P},��������նȾ���
		Q = getQ(runele.globalKV, runele.globalKR, runele.globalKC, P);
		/********************************/
		//��aϵ��
		double a1 = 0, a2 = 0;
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
		{
			b1 += (R[i] * Z[i]);
		}
		double b = b1 / b2;
		/********************************/
		//{P}={Z}+b*{P}
		for (int i = 0; i < (2 * nodeNumber); i++)
		{
			P[i] = Z[i] + b*P[i];
		}
		/********************************/
		max = max_element(begin(aP), end(aP));
	} while (aP[distance(begin(aP), max)] > (1E-8));
	/********************************************************************/
}
//PCGʵ�֣�
void solve_PCG::doPCG_EBE()
{
	load runloa(filename);		runloa.doLoadR();
	/********************************************************************/
	//��ʼ��displacementΪ0����vector�Զ���ʼ��Ϊ0��
	//{R}={loadR}-[K]{displcement}
	vector<double> R = runloa.loadR;
	//{Z}=[M](-1){R}
	vector<double> MN = getMN();
	vector<double> Z = getZ(MN, R);
	//{P}={Z}
	vector<double> P = Z;
	/********************************************************************/
	//iteration begin
	vector<double> Q;
	vector<double> aP((2 * nodeNumber), 0);
	vector<double>::iterator max = max_element(begin(aP), end(aP));//�ҵ����������ֵ�ġ�λ�á���
	do
	{
		/********************************/
		//{Q}=[K]{P},ʹ��element-by-element method;
		Q = getQ(P);
		/********************************/
		//��aϵ��
		double a1 = 0, a2 = 0;
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
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
		for (int i = 0; i < (2 * nodeNumber); i++)
		{
			b1 += (R[i] * Z[i]);
		}
		double b = b1 / b2;
		/********************************/
		//{P}={Z}+b*{P}
		for (int i = 0; i < (2 * nodeNumber); i++)
		{
			P[i] = Z[i] + b*P[i];
		}
		/********************************/
		max = max_element(begin(aP), end(aP));
	} while (aP[distance(begin(aP), max)] > (1E-8));
	/********************************************************************/
}


//�õ�Ԥ������󣬺���������նȾ���ֵ������նȾ���ֵ�кţ�����նȾ���ֵ�кţ���
//����Ԥ������󣬼����Խ�Ԫ�أ�
vector<double> solve_PCG::getMN(vector<double> globalKV, vector<int> globalKR, vector<int> globalKC)
{
	vector<double> MN((2 * nodeNumber), 0);
	for (unsigned int i = 0; i < (globalKV.size()); i++)
	{
		if (globalKR[i] == globalKC[i])
		{
			MN[globalKR[i]] = 1.0 / globalKV[i];
		}
	}
	return MN;
}
//���غ�����
//�õ�Ԥ������󣬺���������
//����Ԥ������󣬼����Խ�Ԫ�أ�
vector<double> solve_PCG::getMN()
{
	elementK runele(filename);
	vector<double> MN((2 * nodeNumber), 0);
	for (int j = 0; j < elementNumber; j++)
	{
		for (int k = 2; k < 6; k++)
		{
			MN[element[j][k] * 2] += runele.getElementK(j)[2 * (k - 2)][2 * (k - 2)];
			MN[element[j][k] * 2 + 1] += runele.getElementK(j)[2 * (k - 2) + 1][2 * (k - 2) + 1];
		}
	}
	//���������Ĵ���
	for (int i = 0; i < boundNumber; i++)
	{
		MN[int(bound[i][1] * 2 + bound[i][2])] = 1E+30;
	}
	for (int i = 0; i < (2 * nodeNumber); i++)
	{
		MN[i] = 1.0 / MN[i];
	}
	return MN;
}


//���ڼ���Z����
vector<double> solve_PCG::getZ(vector<double> MN, vector<double> R)
{
	vector<double> Z((2 * nodeNumber), 0);
	for (int i = 0; i < (2 * nodeNumber); i++)
	{
		Z[i] = MN[i] * R[i];
	}
	return Z;
}


//���ڼ���Q����
vector<double> solve_PCG::getQ(vector<double>globalKV, vector<int>globalKR, vector<int>globalKC, vector<double> P)
{
	vector<double> Q((2 * nodeNumber), 0);
	for (unsigned int i = 0; i < (globalKV.size()); i++)
	{
		Q[globalKR[i]] += globalKV[i] * P[globalKC[i]];
		if (globalKR[i] > globalKC[i])
		{
			Q[globalKC[i]] += globalKV[i] * P[globalKR[i]];//��������������㣬ͬʱע�ⲻҪ�ѶԽ�Ԫ���ظ����㣻
		}
	}
	return Q;
}
//���غ�����
//���ڼ���Q����
vector<double> solve_PCG::getQ(vector<double> P)
{
	elementK runele(filename);
	often runoft;
	vector<double> Q((2 * nodeNumber), 0);
	for (int i = 0; i < elementNumber; i++)
	{
		vector<vector<double>> Pson(8, vector<double>(1, 0));
		for (int j = 2; j < 6; j++)
		{
			Pson[(j - 2) * 2][0] = P[element[i][j] * 2];
			Pson[(j - 2) * 2 + 1][0] = P[element[i][j] * 2 + 1];
		}
		vector<vector<double>> Qson = runoft.matrixM(runele.getElementK(i), Pson);
		for (int j = 2; j < 6; j++)
		{
			Q[element[i][j] * 2] += Qson[(j - 2) * 2][0];
			Q[element[i][j] * 2 + 1] += Qson[(j - 2) * 2 + 1][0];
		}
	}
	//���������Ĵ���
	for (int k = 0; k < boundNumber; k++)
	{
		Q[int(bound[k][1] * 2 + bound[k][2])] += P[int(bound[k][1] * 2 + bound[k][2])] * 1E+30;
	}
	return Q;
}


//���λ��ֵ��FEMout.txt��
void solve_PCG::out()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ�������ļ���
	if (!fileout) {
		cout << "solve_PCG::out File open error!!!" << endl;
		system("pause"); exit(0);
	}
	fileout << "NoodeID" << "\t\t" << "X-Displacement" << "\t\t" << "Y-Displacement" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < nodeNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 2; j++)
		{
			fileout << displacement[2 * i + j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << resetiosflags(ios_base::scientific);
	fileout << "\r\n";
	fileout.close();
}


solve_PCG::~solve_PCG()
{
}


// xiezhuoyu
// mechanics_xzy@163.com