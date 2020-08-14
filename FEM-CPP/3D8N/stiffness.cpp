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


//�κ����Ծֲ������ƫ�������������ֲ�����r���ֲ�����s���ֲ�����t��
vector<vector<vector<double>>> stiffness::Nri(double r = 0.0, double s = 0.0, double t = 0.0)
{
	vector<vector<vector<double>>> Nri(8, vector<vector<double>>(3, vector<double>(1, 0)));
	Nri[0][0][0] = -(1 - s)*(1 - t) / 8, Nri[1][0][0] = (1 - s)*(1 - t) / 8, Nri[2][0][0] = (1 + s)*(1 - t) / 8, Nri[3][0][0] = -(1 + s)*(1 - t) / 8, Nri[4][0][0] = -(1 - s)*(1 + t) / 8, Nri[5][0][0] = (1 - s)*(1 + t) / 8, Nri[6][0][0] = (1 + s)*(1 + t) / 8, Nri[7][0][0] = -(1 + s)*(1 + t) / 8;
	Nri[0][1][0] = -(1 - r)*(1 - t) / 8, Nri[1][1][0] = -(1 + r)*(1 - t) / 8, Nri[2][1][0] = (1 + r)*(1 - t) / 8, Nri[3][1][0] = (1 - r)*(1 - t) / 8, Nri[4][1][0] = -(1 - r)*(1 + t) / 8, Nri[5][1][0] = -(1 + r)*(1 + t) / 8, Nri[6][1][0] = (1 + r)*(1 + t) / 8, Nri[7][1][0] = (1 - r)*(1 + t) / 8;
	Nri[0][2][0] = -(1 - r)*(1 - s) / 8, Nri[1][2][0] = -(1 + r)*(1 - s) / 8, Nri[2][2][0] = -(1 + r)*(1 + s) / 8, Nri[3][2][0] = -(1 - r)*(1 + s) / 8, Nri[4][2][0] = (1 - r)*(1 - s) / 8, Nri[5][2][0] = (1 + r)*(1 - s) / 8, Nri[6][2][0] = (1 + r)*(1 + s) / 8, Nri[7][2][0] = (1 - r)*(1 + s) / 8;
	return Nri;
}


//�Ÿ��Ⱦ���J���������ֲ�����r���ֲ�����s���ֲ�����t����Ԫ���ID����
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


//Ӧ�����B���������ֲ�����r���ֲ�����s���ֲ�����t����Ԫ���ID����
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


//���Ծ���D����������Ԫ���ID����
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
//���ɵ�Ԫ�նȾ���elementK;
Ӧ�ø�˹���֣���ֵ���֣���
��˹��ȡ2*2��Ϊ(-0.577��+0.577)��Ȩϵ��Ϊ(1.0��1.0)��
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
//�鼯���徢�Ⱦ��󣬺���������
ѭ���������е�Ԫ��
���ݵ�Ԫ�ֲ����ţ���Ԫ���Ⱦ������кţ��ҵ�ȫ�ֽ��ţ����徢�Ⱦ������кţ���
��Ԫ��ֵ�Ͷ�Ӧ���кŴ���COO�����С�
******************************************************************/
/**********************************************************************************
����������
 -                        -    -        -       -                -
|E+30��������              |  |��֪λ��1 |     |��֪λ��1*��E+30��|
|����                      |  |����      |     |����              |
|��������                  |  |����      |  =  |����              |
|������������E+30          |  |��֪λ��2 |     |��֪λ��2*��E+30��|
|��������������������      |  |����      |     |����              |
|������������������������  |  |����      |     |����              |
 -                        -    -        -       -                -
**********************************************************************************/
void stiffness::globalStiffness()
{
	for (int k = 0; k < elementNumber; k++)
	{
		vector<vector<double>>  value = elementStiffness(k);
		for (int j = 0; j < 24; j++)
		{
			for (int i = j; i < 24; i++)//ֻ��Ҫ������Ԫ�أ������Խ���Ԫ�أ�
			{
				int ii = match(k, i);//ȫ�ֶ�Ӧ�ֲ���
				int jj = match(k, j);
				if (ii < jj)
				{
					int temp = ii;
					ii = jj;
					jj = temp;
				}//������նȾ����������Ԫ�ػ���������ȥ��
				if ((abs(value[i][j]) > 1E-30))
				{
					COO(value[i][j], ii, jj);//ֻ�������Ԫ�أ�
				}
			}
		}
	}
	//����������
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
//�ֲ���������ȫ�ֽ�����ƥ��,��������Ԫ�ţ��ֲ������ţ���
//����ȫ�ֵĽ����ţ�
�������徢�Ⱦ����������㣬���Ӧ���ĸ�Ԫ��:
										j���
									   3j  3j+1	3j+2  �к�
						i���	3i     xx   xy	xz
								3i+1   yx   yy	yz
								3i+2   zx	zy	zz
						�к�
**************************************************************/
int stiffness::match(int elementID = 0, int i = 0)
{
	int ii;
	if (0 == (i % 3))//ȡ�ࣻ
	{
		ii = 3 * element[elementID][i / 3 + 2];//ȡ�̣�
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
//COO�洢��������Ԫ��ֵ���кţ��кţ���
�洢ϡ�����
�洢���ݣ�1.�����0Ԫ��ֵglobalStiffnessV��value��
2.�����0Ԫ�ص��к�globalStiffnessR��rowNumber��
3.�����0Ԫ�ص��к�globalStiffnessC��columnNumber��
��ֻ�������Ƿ���Ԫ�ء�
**************************************************************/
void stiffness::COO(double value = 0, int i = 0, int j = 0)
{
	static int VRCi = 0;//���ڼ�¼��0Ԫ�ظ�����
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