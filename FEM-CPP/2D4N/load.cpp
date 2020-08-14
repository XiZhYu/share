#include <iostream>
#include <stdlib.h>//exit();
#include <math.h>//pow(),sqrt();
#include <vector>
#include "load.h"
#include "elementk.h"
#include "often.h"
using namespace std;

load::load(string file)
{
	filename = file;
}


/********************************************************************
//�鼯�������������loadR������doLoadR()��
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
void load::doLoadR()
{
	loadR.resize(2 * nodeNumber);//��������������СΪ���ڵ���*�ռ�ά����
	doBodyforce();//����������
	donodeload();//���Ǽ�������
	doFaceload();//����������
	//���ǵ�֧�������Ĵ��ڣ��ô�����������д�������֧���������δ֪����������ⷽ�̲����У�
	for (int i = 0; i < boundNumber; i++)
	{
		loadR[int(bound[i][1] * 2 + bound[i][2])] = bound[i][3] * (1E+30);
	}
}


/**************************************************************************
//���������γɵĵ�Ч������loadRi������doBodyforce()��
	Ӧ�ø�˹���֣���ֵ���֣���
	��˹��ȡ2*2����λ��Ϊ(-0.577350,0.577350)��Ȩϵ��Ϊ(1.000000,1.000000)��
	�Ե�Ԫ���б������������ָ���Ԫ��㣬���鼯�������Ч����������
***************************************************************************/
void load::doBodyforce()
{
	elementK runele(filename);
	often runoft;
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };//���ꣻ
	double H[2] = { 1.0,1.0 };//ϵ����
	for (int k = 0; k < elementNumber; k++)//�������е�Ԫ��
	{
		vector<vector<double>> sonLoadR(8, vector<double>(1, 0));
		vector<vector<double>> P{ {0},{ -material[element[k][1]][3] } };//��������
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				//��ʽ��R=[N]t*[P]*|J|*Hi*Hj���ֱ�������Σ�
				//N(r,s),P(x,y),J(r,s,elementID)��
				sonLoadR = runoft.matrixP(
					runoft.matrixM2(
						runoft.matrixM(
							runoft.matrixT(
								runoft.getN(x[i], x[j]))
							, P)
						, (runoft.matrixDet2(runele.getJacJ(x[i], x[j], k))*H[i] * H[j]))
					, sonLoadR);
			}

		}
		for (int l = 0; l < 8; l++)
		{
			int ll = runoft.match(k, l, filename);//ȫ�ֶ�Ӧ�ֲ���
			loadR[ll] += sonLoadR[l][0];//���ӵ������Ч����������loadRi��
		}
	}
}


/**************************************************************************************
//���������γɵĵ�Ч����������loadRi������donodeload()��
	�Լ��������б������������ָ���������صĵ�Ԫ�Ľ�㣬���鼯�������Ч����������
	����Ҫ�������Ǽ����������ڽ��ͱ�Ե�ϵ�������������롰�����ڵ�Ԫ�ϡ����У�
***************************************************************************************/
void load::donodeload()
{
	often runoft;
	for (int k = 0; k < nodeloadNumber; k++)//�������м��к��أ�
	{
		vector<vector<double>> sonLoadR;
		vector<vector<double>> xy{ { nodeload[k][1] },{ nodeload[k][2] } };//���к���λ�����ꣻ
		vector<vector<double>> P{ { nodeload[k][3] },{ nodeload[k][4] } };//���к��ش�С��
		int elementID = nodeload_Element(k);//�ɼ��������õ�ȷ����ص�Ԫ��
		vector<double> rs = gloToLoc(xy, elementID);//��ȫ������ȷ���ֲ����ꣻ
		//��ʽ��R=[N]t*[P]��
		//N(r,s),P(x,y)��
		sonLoadR = runoft.matrixM(
			runoft.matrixT(
				runoft.getN(rs[0], rs[1]))
			, P);
		for (int l = 0; l < 8; l++)
		{
			int ll = runoft.match(elementID, l, filename);//ȫ�ֶ�Ӧ�ֲ���
			loadR[ll] += sonLoadR[l][0];//���ӵ������Ч����������loadRi��
		}
	}
}


//���������γɵĵ�Ч����������loadRi������doFaceload()��
void load::doFaceload()
{
	often runoft;
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };//���ꣻ
	double H[2] = { 1.0,1.0 };//ϵ����
	double rP = 0, sP = 0;//��˹�㴦�ֲ����ꣻ
	for (int k = 0; k < faceloadNumber; k++)//�������зֲ����أ�
	{
		vector<vector<double>> sonLoadR(8, vector<double>(1, 0));
		vector<int> elementLine = faceload_Line(k);//ȷ���ֲ���λ�ڼ�����ص�Ԫ�Լ��ڼ��űߣ�
		double P1[2] = { faceload[k][3],faceload[k][5] };//�ֲ���1��㴦����ֵ��
		double P2[2] = { faceload[k][4],faceload[k][6] };//�ֲ���2��㴦����ֵ��
		double Px[2] = { 0,0 };//��˹�㴦x�������ֵ��
		double Py[2] = { 0,0 };//��˹�㴦y�������ֵ��
		for (int i = 0; i < 2; i++)//��˹���֣�
		{
			vector<vector<double>> N;
			switch (elementLine[1])//���ݷֲ���λ�ڵ�Ԫ�߱���Լ���������ȷ����˹�㴦�ľֲ����ꡢ�κ�������[N]����˹�㴦�ĺ���ֵ��
			{
			case 1:
				rP = x[i]; sP = -1;
				N = runoft.getN(rP, sP);
				Px[i] = N[0][0] * P1[0] + N[0][2] * P2[0]; Py[i] = N[0][0] * P1[1] + N[0][2] * P2[1];
				break;
			case -1:
				rP = x[i]; sP = -1;
				N = runoft.getN(rP, sP);
				Px[i] = N[0][0] * P2[0] + N[0][2] * P1[0]; Py[i] = N[0][0] * P2[1] + N[0][2] * P1[1];
				break;
			case 2:
				rP = 1; sP = x[i];
				N = runoft.getN(rP, sP);
				Px[i] = N[0][2] * P1[0] + N[0][4] * P2[0]; Py[i] = N[0][2] * P1[1] + N[0][4] * P2[1];
				break;
			case -2:
				rP = 1; sP = x[i];
				N = runoft.getN(rP, sP);
				Px[i] = N[0][2] * P2[0] + N[0][4] * P1[0]; Py[i] = N[0][2] * P2[1] + N[0][4] * P1[1];
				break;
			case 3:
				rP = x[i]; sP = 1;
				N = runoft.getN(rP, sP);
				Px[i] = N[0][4] * P1[0] + N[0][6] * P2[0]; Py[i] = N[0][4] * P1[1] + N[0][6] * P2[1];
				break;
			case -3:
				rP = x[i]; sP = 1;
				N = runoft.getN(rP, sP);
				Px[i] = N[0][4] * P2[0] + N[0][6] * P1[0]; Py[i] = N[0][4] * P2[1] + N[0][6] * P1[1];
				break;
			case 4:
				rP = -1; sP = x[i];
				N = runoft.getN(rP, sP);
				Px[i] = N[0][6] * P1[0] + N[0][0] * P2[0]; Py[i] = N[0][6] * P1[1] + N[0][0] * P2[1];
				break;
			case -4:
				rP = -1; sP = x[i];
				N = runoft.getN(rP, sP);
				Px[i] = N[0][6] * P2[0] + N[0][0] * P1[0]; Py[i] = N[0][6] * P2[1] + N[0][0] * P1[1];
				break;
			default://��ȷ��ֻ���ܳ������������
				cout << "load::doFaceload error!!!" << endl;
				system("pause"); exit(0);
				break;
			}
			vector<vector<double>> P{ { Px[i] },{ Py[i] } };//��˹�㴦�ĺ���ֵ��
			//��ʽ��R=[N]t*[P]*xishu*Hi���2�Σ�
			//N(r,s),P(x,y)��
			double xishu = sqrt(pow(0.5*(node[faceload[k][1]][1] - node[faceload[k][2]][1]), 2)
				+ pow(0.5*(node[faceload[k][1]][2] - node[faceload[k][2]][2]), 2));//��ϵ��Ϊȫ������ת��Ϊ�ֲ�����ʱ�����������ӵģ�
			sonLoadR = runoft.matrixP(
				runoft.matrixM2(
					runoft.matrixM(
						runoft.matrixT(
							runoft.getN(rP, sP))
						, P)
					, (xishu*H[i]))
				, sonLoadR);
		}
		for (int l = 0; l < 8; l++)
		{
			int ll = runoft.match(elementLine[0], l, filename);//ȫ�ֶ�Ӧ�ֲ���
			loadR[ll] += sonLoadR[l][0];//���ӵ������Ч����������loadRi��
		}
	}
}


/***********************************************************************************************
//�ɼ�������ŵõ���صĵ�Ԫ��ţ�����nodeloadElement(���������)��
	���������õ��뵥Ԫ�ĸ���������γ�4�������Σ�������ʽ����ʱ���ƶ���Ϊ˳�򣩼��������������
	�����������õ㲻�ڵ�Ԫ�ڲ����Ե����4���������������1��Ϊ������
	�����������õ��ڵ�Ԫ�ڲ����Ե����4�����������ȫ�����ڵ���0.
************************************************************************************************/
int load::nodeload_Element(int nodeloadID = 0)
{
	often runoft;
	bool panju = true;
	for (int i = 0; i < elementNumber; i++)
	{
		vector<vector<double>> a1{ { 1 ,nodeload[nodeloadID][1],nodeload[nodeloadID][2] },
									{ 1 ,node[element[i][2]][1] ,node[element[i][2]][2] },
									{1,node[element[i][3]][1],node[element[i][3]][2] } };
		vector<vector<double>> a2{ {  1, nodeload[nodeloadID][1],  nodeload[nodeloadID][2] },
									{  1,  node[element[i][3]][1],  node[element[i][3]][2] },
									{  1,  node[element[i][4]][1],  node[element[i][4]][2] } };
		vector<vector<double>> a3{ {  1,  nodeload[nodeloadID][1],  nodeload[nodeloadID][2] },
									{  1, node[element[i][4]][1],  node[element[i][4]][2] },
									{ 1,  node[element[i][5]][1],  node[element[i][5]][2] } };
		vector<vector<double>> a4{ {  1,  nodeload[nodeloadID][1],  nodeload[nodeloadID][2] },
									{ 1,  node[element[i][5]][1],  node[element[i][5]][2] },
									{  1,  node[element[i][2]][1], node[element[i][2]][2] } };
		double A1 = runoft.matrixDet3(a1);
		double A2 = runoft.matrixDet3(a2);
		double A3 = runoft.matrixDet3(a3);
		double A4 = runoft.matrixDet3(a4);
		if ((A1 >= 0) && (A2 >= 0) && (A3 >= 0) && (A4 >= 0))
		{
			panju = false;
			return i;//�ҵ����˳���
		}
	}
	if (true == panju)
	{
		cout << "load::nodeloadElement error!!!" << endl;
		system("pause"); exit(0);
	}
	return 0;
}


/********************************************************************
//��ȫ������ȷ���ֲ����꣬����gloToLoc(ȫ�����꣬���ڵ�Ԫ���)��
//���ؾֲ���������rs[4]={r,s}��
	ȡ(r=0,s=0)Ϊ��ֵ�������Ÿ������ʽ��ȫ�����꣨��ֵ���;ֲ����꣨��ֵ���Ĺ�ϵ��
	������ֱ����������Ҫ��
**********************************************************************/
vector<double> load::gloToLoc(vector<vector<double>> xy, int elementID = 0)
{
	elementK runele(filename);
	often runoft;
	vector<vector<double>> nodexy{ { node[element[elementID][2]][1],node[element[elementID][2]][2] },
									{ node[element[elementID][3]][1],node[element[elementID][3]][2] },
									{ node[element[elementID][4]][1],node[element[elementID][4]][2] },
									{ node[element[elementID][5]][1],node[element[elementID][5]][2] } };//��Ԫ4��������ꣻ
	vector<double> rs(2, 0);//�ֲ����ꣻ
	vector<vector<double>> rscha(2, vector<double>(1, 0));//�ֲ������������Ĳ�ֵ,���ֵΪ1��
	do//������ֱ���ֲ������ֵ���㾫��Ҫ��
	{
		vector<vector<double>> N = runoft.getN(rs[0], rs[1]);
		vector<vector<double>> xynew{ { N[0][0] * nodexy[0][0] + N[0][2] * nodexy[1][0] + N[0][4] * nodexy[2][0] + N[0][6] * nodexy[3][0] },
										{ N[0][0] * nodexy[0][1] + N[0][2] * nodexy[1][1] + N[0][4] * nodexy[2][1] + N[0][6] * nodexy[3][1] } };
		rscha = runoft.matrixM(
			runoft.matrixN2(runele.getJacJ(rs[0], rs[1], elementID))
			, runoft.matrixP2(xy, xynew));
		rs[0] += rscha[0][0]; rs[1] += rscha[1][0];
	} while ((rscha[0][0] > 0.001) || (rscha[1][0] > 0.001));
	return rs;
}


//���������ȷ�����ڵ�Ԫ�ı߱�ţ�����faceloadLine(�ֲ������)��
//���ص�Ԫ��ţ�
//����1��2��3��4�ֱ����ֲ�����ϵ��ʱ����1��2������ڱߣ�2��3��㡭����4��1������ڱߣ���ֵΪ������Ϊ˳ʱ�뷽��ıߣ�
//�磺re[0]=elementID,re[1]=lineID;
vector<int> load::faceload_Line(int faceloadID)
{
	vector<int> result(2, 0);
	for (int i = 0; i < elementNumber; i++)//�ظ��Աȵ�Ԫ�Ľ����Ϣ�ͷֲ����Ľ����Ϣ��ֱ���Ǻϣ�
	{
		vector<double> f{ { faceload[faceloadID][1] ,faceload[faceloadID][2] } };
		vector<int> e{ { element[i][2] ,element[i][3],element[i][4], element[i][5] } };
		if ((abs(f[0] - e[0]) < 1E-10) && (abs(f[1] - e[1]) < 1E-10)) { result[0] = i; result[1] = 1; break; }
		if ((abs(f[1] - e[0]) < 1E-10) && (abs(f[0] - e[1]) < 1E-10)) { result[0] = i; result[1] = -1; break; }
		if ((abs(f[0] - e[1]) < 1E-10) && (abs(f[1] - e[2]) < 1E-10)) { result[0] = i; result[1] = 2; break; }
		if ((abs(f[1] - e[1]) < 1E-10) && (abs(f[0] - e[2]) < 1E-10)) { result[0] = i; result[1] = -2; break; }
		if ((abs(f[0] - e[2]) < 1E-10) && (abs(f[1] - e[3]) < 1E-10)) { result[0] = i; result[1] = 3; break; }
		if ((abs(f[1] - e[2]) < 1E-10) && (abs(f[0] - e[3]) < 1E-10)) { result[0] = i; result[1] = -3; break; }
		if ((abs(f[0] - e[3]) < 1E-10) && (abs(f[1] - e[0]) < 1E-10)) { result[0] = i; result[1] = 4; break; }
		if ((abs(f[1] - e[3]) < 1E-10) && (abs(f[0] - e[0]) < 1E-10)) { result[0] = i; result[1] = -4; break; }
	}
	if (0 == result[1])
	{
		cout << "load::faceloadLine no result!!!" << endl;
		system("pause"); exit(0);
	}
	return result;
}

load::~load(void)
{
}


// xiezhuoyu
// mechanics_xzy@163.com