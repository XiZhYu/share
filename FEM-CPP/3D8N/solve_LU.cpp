#include "solve_LU.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
using namespace std;

solve_LU::solve_LU()
{
}


solve_LU::~solve_LU()
{
}


vector<double> solve_LU::getDisplacement()
{
	displacement.clear();
	displacement.resize(nodeNumber * 3);
	//LU�����Ԥ�����˴�û��ʹ��ѹ�������ʽ�����������������ǵ�Ԫ�أ�
	LU.clear();
	LU.resize(int((nodeNumber * 3)*(nodeNumber * 3 + 1)*0.5));
	for (int i = 0; i < (nodeNumber * 3); i++)
	{
		LU[getLUij(i, i)] = 1;//������Ԫ��Ϊ����
	}
	doLU();
	forward();
	backwards();
	return displacement;
}


//���ڲ���LU������getLUij(�кţ��к�)��
//��������õ����к�ֵ�����ش洢��ʽ�µ�λ��ָ�ꣻ
int solve_LU::getLUij(int row = 0, int col = 0)
{
	if (row < col)
	{
		cout << "solve_LU::getLUij input error!!!" << endl;
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
	//����L�����0���е�ֵ��
	for (int i = 0; i < (nodeNumber * 3); i++)
	{
		bool panju = true;
		for (unsigned int k = 0; k < (globalStiffnessR.size()); k++)
		{
			if ((i == globalStiffnessR.at(k)) && (0 == globalStiffnessC.at(k)))
			{
				LU[getLUij(i, 0)] = globalStiffnessV.at(k);
				panju = false;
			}
		}
		if (true == panju) { LU[getLUij(i, 0)] = 0; }
	}
	for (int j = 1; j < (nodeNumber * 3); j++)//�ӵ�1���п�ʼ���㣻
	{
		for (int i = j; i < (nodeNumber * 3); i++)
		{
			double sum = 0;
			bool panju = true;
			for (int p = 0; p < j; p++)
			{
				if (abs(LU[getLUij(p, p)]) < 1E-15)//��0����ʾ����
				{
					cout << "solve_LU::doLU error!!!" << endl;
					system("pause"); exit(0);
				}
				sum = sum + (LU[getLUij(i, p)] * LU[getLUij(j, p)]) / LU[getLUij(p, p)];
			}
			for (unsigned int k = 0; k < (globalStiffnessR.size()); k++)
			{
				if ((i == globalStiffnessR.at(k)) && (j == globalStiffnessC.at(k)))
				{
					LU[getLUij(i, j)] = globalStiffnessV.at(k) - sum;
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
	displacement[0] = loadR[0] / LU[getLUij(0, 0)];
	for (int i = 1; i < (nodeNumber * 3); i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum = sum + LU[getLUij(i, j)] * displacement[j];
		}
		if (abs(LU[getLUij(i, i)]) < 1E-15)//��0����ʾ����
		{
			cout << "solve_LU::forward error!!!" << endl;
			system("pause"); exit(0);
		}
		displacement[i] = (loadR[i] - sum) / LU[getLUij(i, i)];
	}
}


//�ش����[U][displacement]=[g]������backwards������
void solve_LU::backwards()
{
	for (int i = (nodeNumber * 3 - 1); i >= 0; i--)
	{
		double sum = 0;
		for (int j = (i + 1); j < (nodeNumber * 3); j++)
		{
			sum = sum + (LU[getLUij(j, i)] / LU[getLUij(i, i)])*displacement[j];
		}
		if (abs(LU[getLUij(i, i)]) < 1E-15)//��0����ʾ����
		{
			cout << "solve_LU::backwards error!!!" << endl;
			system("pause"); exit(0);
		}
		displacement[i] = displacement[i] - sum;
	}
}


// xiezhuoyu
// mechanics_xzy@163.com