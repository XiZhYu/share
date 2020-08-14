#include <iostream>
#include <stdlib.h>//exit();
#include <fstream>
#include <iomanip>
#include <vector>
#include "solve_LU.h"
#include "elementk.h"
#include "load.h"
using namespace std;

solve_LU::solve_LU(string file)
{
	filename = file;
}


vector<double> solve_LU::getDisplacement()
{
	displacement.resize(nodeNumber * 2);
	//LU�����Ԥ�����˴�û��ʹ��ѹ�������ʽ�����������������ǵ�Ԫ�أ�
	LU.resize(int((nodeNumber * 2)*(nodeNumber * 2 + 1)*0.5));
	for (int i = 0; i < (nodeNumber * 2); i++)
	{
		LU[getLUij(i, i)] = 1;//������Ԫ��Ϊ����
	}
	doLU();
	forward();
	backwards();
	out();
	return displacement;
}


//���ڲ���LU������getLUij(�кţ��к�)��
//��������õ����к�ֵ�����ش洢��ʽ�µ�λ��ָ�ꣻ
int solve_LU::getLUij(int row = 0, int col = 0)
{
	if (row < col)
	{
		cout << "solve::getLUij input error!!!" << endl;
		system("pause"); exit(0);
	}
	int sum = 0;
	for (int i = 0; i < (row + 1); i++)
	{
		sum += i;
	}
	return (sum + col);
}


//�ֽ�õ�L��U�����е�L���󣬺���LUglobalK()��
//������������LU[]�����У�
void solve_LU::doLU()
{
	elementK runele(filename);	runele.doGlobalK();//ֻ������һ�Σ�
	//����L�����0���е�ֵ��
	for (int i = 0; i < (nodeNumber * 2); i++)
	{
		bool panju = true;
		for (unsigned int k = 0; k < (runele.globalKR.size()); k++)
		{
			if ((i == runele.globalKR.at(k)) && (0 == runele.globalKC.at(k)))
			{
				LU[getLUij(i, 0)] = runele.globalKV.at(k);
				panju = false;
			}
		}
		if (true == panju) { LU[getLUij(i, 0)] = 0; }
	}
	for (int j = 1; j < (nodeNumber * 2); j++)//�ӵ�1���п�ʼ���㣻
	{
		for (int i = j; i < (nodeNumber * 2); i++)
		{
			double sum = 0;
			bool panju = true;
			for (int p = 0; p < j; p++)
			{
				if (abs(LU[getLUij(p, p)]) < 1E-15)//��0����ʾ����
				{
					cout << "solve::doLU error!!!" << endl;
					system("pause"); exit(0);
				}
				sum = sum + (LU[getLUij(i, p)] * LU[getLUij(j, p)]) / LU[getLUij(p, p)];
			}
			for (unsigned int k = 0; k < (runele.globalKR.size()); k++)
			{
				if ((i == runele.globalKR.at(k)) && (j == runele.globalKC.at(k)))
				{
					LU[getLUij(i, j)] = runele.globalKV.at(k) - sum;
					panju = false;
				}
			}
			if (true == panju) { LU[getLUij(i, j)] = 0 - sum; }
		}
	}
}


//ǰ�����[L][g]=[R]�е�[g]������forward������
void solve_LU::forward()
{
	load runloa(filename);	runloa.doLoadR();//ֻ������һ�Σ�
	displacement[0] = runloa.loadR[0] / LU[getLUij(0, 0)];
	for (int i = 1; i < (nodeNumber * 2); i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum = sum + LU[getLUij(i, j)] * displacement[j];
		}
		if (abs(LU[getLUij(i, i)]) < 1E-15)//��0����ʾ����
		{
			cout << "solve::forward error!!!" << endl;
			system("pause"); exit(0);
		}
		displacement[i] = (runloa.loadR[i] - sum) / LU[getLUij(i, i)];
	}
}


//�ش����[U][displacement]=[g]������backwards������
void solve_LU::backwards()
{
	for (int i = (nodeNumber * 2 - 1); i >= 0; i--)
	{
		double sum = 0;
		for (int j = (i + 1); j < (nodeNumber * 2); j++)
		{
			sum = sum + (LU[getLUij(j, i)] / LU[getLUij(i, i)])*displacement[j];
		}
		if (abs(LU[getLUij(i, i)]) < 1E-15)//��0����ʾ����
		{
			cout << "solve::backwards error!!!" << endl;
			system("pause"); exit(0);
		}
		displacement[i] = displacement[i] - sum;
	}
}

//���λ��ֵ��FEMout.txt��
void solve_LU::out()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ�������ļ���
	if (!fileout) {
		cout << "stress::out File open error!!!" << endl;
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

solve_LU::~solve_LU(void)
{
}


// xiezhuoyu
// mechanics_xzy@163.com