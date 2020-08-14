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
//组集整体结点荷载列阵loadR，函数doLoadR()；
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
void load::doLoadR()
{
	loadR.resize(2 * nodeNumber);//整体结点荷载列阵大小为：节点数*空间维数；
	doBodyforce();//考虑体力；
	donodeload();//考虑集中力；
	doFaceload();//考虑面力；
	//考虑到支座反力的存在，用大数法处理进行处理，将含支座反力这个未知条件不纳入解方程步骤中；
	for (int i = 0; i < boundNumber; i++)
	{
		loadR[int(bound[i][1] * 2 + bound[i][2])] = bound[i][3] * (1E+30);
	}
}


/**************************************************************************
//处理体力形成的等效结点荷载loadRi，函数doBodyforce()；
	应用高斯积分（数值积分）。
	高斯点取2*2个，位置为(-0.577350,0.577350)，权系数为(1.000000,1.000000)；
	对单元进行遍历，将体力分给单元结点，再组集到整体等效结点荷载列阵。
***************************************************************************/
void load::doBodyforce()
{
	elementK runele(filename);
	often runoft;
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };//坐标；
	double H[2] = { 1.0,1.0 };//系数；
	for (int k = 0; k < elementNumber; k++)//遍历所有单元；
	{
		vector<vector<double>> sonLoadR(8, vector<double>(1, 0));
		vector<vector<double>> P{ {0},{ -material[element[k][1]][3] } };//荷载列阵；
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				//公式：R=[N]t*[P]*|J|*Hi*Hj并分别求和两次；
				//N(r,s),P(x,y),J(r,s,elementID)；
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
			int ll = runoft.match(k, l, filename);//全局对应局部；
			loadR[ll] += sonLoadR[l][0];//叠加到整体等效结点荷载列阵loadRi；
		}
	}
}


/**************************************************************************************
//处理集中力形成的等效结点荷载列阵loadRi，函数donodeload()；
	对集中力进行遍历，将体力分给集中力相关的单元的结点，再组集到整体等效结点荷载列阵；
	不需要单独考虑集中力作用在结点和边缘上的情况，都已纳入“作用在单元上”其中；
***************************************************************************************/
void load::donodeload()
{
	often runoft;
	for (int k = 0; k < nodeloadNumber; k++)//遍历所有集中荷载；
	{
		vector<vector<double>> sonLoadR;
		vector<vector<double>> xy{ { nodeload[k][1] },{ nodeload[k][2] } };//集中荷载位置坐标；
		vector<vector<double>> P{ { nodeload[k][3] },{ nodeload[k][4] } };//集中荷载大小；
		int elementID = nodeload_Element(k);//由集中力作用点确定相关单元；
		vector<double> rs = gloToLoc(xy, elementID);//由全局坐标确定局部坐标；
		//公式：R=[N]t*[P]；
		//N(r,s),P(x,y)；
		sonLoadR = runoft.matrixM(
			runoft.matrixT(
				runoft.getN(rs[0], rs[1]))
			, P);
		for (int l = 0; l < 8; l++)
		{
			int ll = runoft.match(elementID, l, filename);//全局对应局部；
			loadR[ll] += sonLoadR[l][0];//叠加到整体等效结点荷载列阵loadRi；
		}
	}
}


//处理面力形成的等效结点荷载列阵loadRi，函数doFaceload()；
void load::doFaceload()
{
	often runoft;
	double x[2] = { -sqrt(1.0 / 3.0),+sqrt(1.0 / 3.0) };//坐标；
	double H[2] = { 1.0,1.0 };//系数；
	double rP = 0, sP = 0;//高斯点处局部坐标；
	for (int k = 0; k < faceloadNumber; k++)//遍历所有分布荷载；
	{
		vector<vector<double>> sonLoadR(8, vector<double>(1, 0));
		vector<int> elementLine = faceload_Line(k);//确定分布力位于几号相关单元以及第几号边；
		double P1[2] = { faceload[k][3],faceload[k][5] };//分布力1结点处荷载值；
		double P2[2] = { faceload[k][4],faceload[k][6] };//分布力2结点处荷载值；
		double Px[2] = { 0,0 };//高斯点处x方向荷载值；
		double Py[2] = { 0,0 };//高斯点处y方向荷载值；
		for (int i = 0; i < 2; i++)//高斯积分；
		{
			vector<vector<double>> N;
			switch (elementLine[1])//根据分布力位于单元边编号以及正负，可确定高斯点处的局部坐标、形函数矩阵[N]、高斯点处的荷载值；
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
			default://正确则只可能出现上述情况；
				cout << "load::doFaceload error!!!" << endl;
				system("pause"); exit(0);
				break;
			}
			vector<vector<double>> P{ { Px[i] },{ Py[i] } };//高斯点处的荷载值；
			//公式：R=[N]t*[P]*xishu*Hi求和2次；
			//N(r,s),P(x,y)；
			double xishu = sqrt(pow(0.5*(node[faceload[k][1]][1] - node[faceload[k][2]][1]), 2)
				+ pow(0.5*(node[faceload[k][1]][2] - node[faceload[k][2]][2]), 2));//该系数为全局坐标转换为局部坐标时，积分所附加的；
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
			int ll = runoft.match(elementLine[0], l, filename);//全局对应局部；
			loadR[ll] += sonLoadR[l][0];//叠加到整体等效结点荷载列阵loadRi；
		}
	}
}


/***********************************************************************************************
//由集中力编号得到相关的单元编号，函数nodeloadElement(集中力编号)；
	集中力作用点与单元四个结点连线形成4个三角形，用行列式（逆时针绕顶点为顺序）计算三角形面积；
	若集中力作用点不在单元内部或边缘，则4个三角形面积有且1个为负数；
	若集中力作用点在单元内部或边缘，则4个三角形面积全都大于等于0.
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
			return i;//找到即退出；
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
//由全局坐标确定局部坐标，函数gloToLoc(全局坐标，所在单元编号)；
//返回局部坐标列阵rs[4]={r,s}；
	取(r=0,s=0)为初值，利用雅阁比行列式中全局坐标（差值）和局部坐标（差值）的关系；
	迭代，直至精度满足要求。
**********************************************************************/
vector<double> load::gloToLoc(vector<vector<double>> xy, int elementID = 0)
{
	elementK runele(filename);
	often runoft;
	vector<vector<double>> nodexy{ { node[element[elementID][2]][1],node[element[elementID][2]][2] },
									{ node[element[elementID][3]][1],node[element[elementID][3]][2] },
									{ node[element[elementID][4]][1],node[element[elementID][4]][2] },
									{ node[element[elementID][5]][1],node[element[elementID][5]][2] } };//单元4个结点坐标；
	vector<double> rs(2, 0);//局部坐标；
	vector<vector<double>> rscha(2, vector<double>(1, 0));//局部坐标迭代引起的差值,设初值为1；
	do//迭代，直至局部坐标差值满足精度要求；
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


//由面力结点确定所在单元的边编号，函数faceloadLine(分布力编号)；
//返回单元编号；
//返回1、2、3、4分别代表局部坐标系逆时针下1、2结点所在边，2、3结点……，4、1结点所在边；若值为负，则为顺时针方向的边；
//如：re[0]=elementID,re[1]=lineID;
vector<int> load::faceload_Line(int faceloadID)
{
	vector<int> result(2, 0);
	for (int i = 0; i < elementNumber; i++)//重复对比单元的结点信息和分布力的结点信息，直至吻合；
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