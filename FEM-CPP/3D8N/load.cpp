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
//���������γɵĽ���������{NLR}������nodeloadR()��
�������ڽ���ϵļ��������鼯���������������
***************************************************************************************/
void load::nodeloadR()
{
	for (int k = 0; k < nodeloadNumber; k++)//�������м��к��أ�
	{
		vector<double> P{ nodeload[k][1] , nodeload[k][2] ,nodeload[k][3] };//���к��ش�С��
		int ii = 3 * int(nodeload[k][0]);//ȫ�ֶ�Ӧ�ֲ���
		for (int i = 0; i < 3; i++)
		{
			loadR[ii + i] += P[i];//���ӵ��������������loadRi��
		}
	}
}


//���������γɵĽ���������{FLR}������faceloadR()��
void load::faceloadR()
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//���ꣻ
	double H[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };//ϵ����
	for (int k = 0; k < faceloadNumber; k++)//�������зֲ����أ�
	{
		vector<vector<double>> FLR(24, vector<double>(1, 0));
		vector<int> elementID_FaceID = faceload_ElementID_FaceID(k);//ȷ���ֲ���λ�ڼ�����ص�Ԫ�Լ��ڼ����棻
		int faceID = elementID_FaceID[1];
		vector<vector<double>> GPP = GaussPiontP(k);//��˹�㴦�ĺ���ֵ��
		for (int i = 0; i < 3; i++)//��˹���֣�
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
								objM.matrixT(N(x[j], 1.0, x[i]))//ע��ij˳��
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
				default://��ȷ��ֻ���ܳ������������
					cout << "load::faceloadR error!!!" << endl;
					system("pause"); exit(0);
					break;
				}
			}
		}
		for (int l = 0; l < 24; l++)
		{
			int ll = objS.match(elementID_FaceID[0], l);//ȫ�ֶ�Ӧ�ֲ���
			loadR[ll] += FLR[l][0];//���ӵ������Ч����������loadRi��
		}
	}
}


/**************************************************************************
//���������γɵĽ���������{BLR}������bodyloadR()��
Ӧ�ø�˹���֣���ֵ���֣���
��˹��ȡ2*2����λ��Ϊ(-0.577350,0.577350)��Ȩϵ��Ϊ(1.000000,1.000000)��
�Ե�Ԫ���б������������ָ���Ԫ��㣬���鼯���������������
***************************************************************************/
void load::bodyloadR()
{
	double x[2] = { -sqrt(1.0 / 3.0), +sqrt(1.0 / 3.0) };//���ꣻ
	double H[2] = { 1.0,1.0 };//ϵ����
	for (int q = 0; q < elementNumber; q++)//�������е�Ԫ��
	{
		vector<vector<double>> BLR(24, vector<double>(1, 0));
		vector<vector<double>> P{ { 0 },{0}, { -material[element[q][1]][2] } };//��������
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					//��ʽ��R=[N]t*[P]*|J|*Hi*Hj*Hk���ֱ�������Σ�
					//N(r,s,t),P(x,y),J(r,s,t,elementID)��
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
			int ll = objS.match(q, l);//ȫ�ֶ�Ӧ�ֲ���
			loadR[ll] += BLR[l][0];//���ӵ������Ч����������loadRi��
		}
	}
}


/********************************************************************
//�鼯�������������loadR������getLoadR()��
����������
 -                        -    -        -       -                -
|E+30��������              |  |��֪λ��1 |     |��֪λ��1*��E+30��|		��Ҫ������������������������֧��������
|����                      |  |����      |     |����              |
|��������                  |  |����      |  =  |����              |
|������������E+30          |  |��֪λ��2 |     |��֪λ��2*��E+30��|
|��������������������      |  |����      |     |����              |
|������������������������  |  |����      |     |����              |
 -                        -    -        -       -                -
*********************************************************************/
void load::getLoadR()
{
	loadR.clear();
	loadR.resize(3 * nodeNumber);//��������������СΪ���ڵ���*�ռ�ά����
	nodeloadR();//���Ǽ�������
	faceloadR();//����������
	bodyloadR();//����������
	//���ǵ�֧�������Ĵ��ڣ��ô�����������д�������֧���������δ֪����������ⷽ�̲����У�
	for (int i = 0; i < boundNumber; i++)
	{
		loadR[int(bound[i][0] * 3 + bound[i][1])] = bound[i][2] * (1E+30);
	}
}


//�κ�������N��������r��s��t��
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


//ȷ���ֲ����뵥Ԫ�Լ���Ĺ�ϵ��������faceload_ElementID_FaceID��
//���طֲ������ڵĵ�ԪID����ID
vector<int> load::faceload_ElementID_FaceID(int faceloadID = 0)
{
	vector<int> result(2, 0);
	for (int i = 0; i < elementNumber; i++)//�ظ��Աȵ�Ԫ�Ľ����Ϣ�ͷֲ����Ľ����Ϣ��ֱ���Ǻϣ�
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


//����3*3����µĸ�˹�㴦�ĺ���ֵ
vector<vector<double>> load::GaussPiontP(int faceloadID = 0)
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//���ꣻ
	vector<vector<double>> GPP(3, vector<double>(3, 0.0));
	for (int i = 0; i < 3; i++)//�����κ���������
	{
		for (int j = 0; j < 3; j++)
		{
			GPP[i][j] = faceload[faceloadID][4] * (1 - x[i])*(1 - x[j]) / 4 + faceload[faceloadID][5] * (1 + x[i])*(1 - x[j]) / 4
				+ faceload[faceloadID][6] * (1 + x[i])*(1 + x[j]) / 4 + faceload[faceloadID][7] * (1 - x[i])*(1 + x[j]) / 4;
		}
	}
	return GPP;
}


//���㹫ʽ�е��������������
//����{A}
vector<vector<double>> load::A(int faceID = 1, int ii = 0, int jj = 0, int elementID = 0)
{
	double x[3] = { -sqrt(3.0 / 5.0), 0.0, +sqrt(3.0 / 5.0) };//���ꣻ
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
	default://��ȷ��ֻ���ܳ������������
		cout << "load::A error!!!" << endl;
		system("pause"); exit(0);
		break;
	}
	return result;
}

// xiezhuoyu
// mechanics_xzy@163.com