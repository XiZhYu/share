#include <iostream>//输入输出；
#include <fstream>//文本输入输出；
#include <string>//字符串；
#include <stdlib.h>//报错，exit();
#include <time.h> //计时；
#include <vector>//数组；
#include "stress.h"
using namespace std;

/*****************************************************************************
读取数据文件：
	结点信息：坐标；
	单元信息：材料ID、相关结点ID；
	材料信息：弹性模量、泊松比、重度；
	约束信息：相关结点ID、结点位移方向、位移大小；
	面力信息：相关结点ID、面力大小；
	集中力信息：力作用点位置、力大小；
*****************************************************************************/
void read(string filename);

/*****************************************************************************
结点信息：结点编号、坐标,结点数；
提取节点坐标node[i][j],i为结点编号，j为0、1、2分别对应结点编号、x坐标、y坐标;
对应格式：
			node	2
			0	0	0
			1	0	1
*****************************************************************************/
vector<vector<double>> node;
int nodeNumber = 0;//结点数目；

/*****************************************************************************
单元信息：单元编号、材料ID、相关结点ID,单元数；
提取单元信息element[i][j]，i为单元编号，j为0到5分别对应单元编号、材料、相关四个节点编号;
对应格式：
			ELEMENT	2
			0	0	0	3	1   2
			1	0	1	3	2   4
*****************************************************************************/
vector<vector<int>> element;
int elementNumber = 0;//单元数目；

/*****************************************************************************
材料信息：材料编号、弹性模量、泊松比、重度,材料数；
提取材料信息,material[i][j];i为材料编号，j为0到3分别对应材料编号、弹性模量、泊松比、重度；
对应格式：
			MATERIAL	1
			0	200000	0.3	77
*****************************************************************************/
vector<vector<double>> material;
int materialNumber = 0;//材料数目；

/*****************************************************************************
约束信息：约束编号、相关结点ID、结点位移方向、位移大小，约束数；
提取约束信息Bound[i][j]，i为约束编号，j为0到3分别对应约束编号、结点编号、xy方向、约束位移值;
对应格式：
			BOUND	2
			0	0	1	0
			1	0	2	0
*****************************************************************************/
vector<vector<double>> bound;
int boundNumber = 0;//约束数目；

/*****************************************************************************
分布力信息：分布力编号、相关结点ID、作用方向、荷载大小，分布荷载数；
提取分布荷载信息faceload[i][j]，i为分布力编号，j为0到6分别为结点1ID、结点2ID、结点1的水平荷载值、结点2的水平荷载值、结点1的竖直荷载值、结点2的竖直荷载值；
对应格式：
			Faceload	1
			0	1	2	0	0	10	20
*****************************************************************************/
vector<vector<double>> faceload;
int faceloadNumber = 0;//面力数目；

/*****************************************************************************
集中力信息：集中力编号、约束编号、力x方向、力y方向、力x分量、力y分量；
提取集中荷载信息nodeload[i][j]，i为约束编号，j为0到4分别集中力编号、作用点x坐标、y坐标、力x分量、力y分量；
对应格式：
			nodeForce	2
			0		3.2		4.8		0		10
			1		1.5		9.7		10		0
*****************************************************************************/
vector<vector<double>> nodeload;
int nodeloadNumber = 0;//集中力数目；

//线性方程组的解法，LU:LU分解-直接法，PCG:预处理共轭梯度法-迭代法；
string solve_choice;
//PCG的解决方法，Normal:正常组成整体刚度矩阵，EBE:逐步单元法；
string PCG_choice;

void main()
{
	cout << "****************************************" << endl << "Program begin." << endl << endl;

	/*****************************************************************************
	主程序控制段；
	*****************************************************************************/
	string filenamein, filenameout;
	cout << "Enter the inputDataFile name:";	cin >> filenamein;
	cout << endl << "Choose your method: PCG or LU" << "\t\t" << "Choice:";		cin >> solve_choice;
	if ("LU" != solve_choice)
	{
		cout << "Choose your PCG_method: Normal(N) or Element-By-Element(EBE)" << "\t\t" << "PCG_Choice:";		cin >> PCG_choice;
	}
	cout << endl << "running……" << endl << endl;
	long time01 = clock();  //计时； 
	//读取数据文件；
	read(filenamein);
	//创建一个输出文件；
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::out | ios_base::binary | ios_base::trunc);//输出|二进制|删除原内容；
	//计算单元中心处应力；
	stress run6(filenamein);	run6.doStress();
	fileout.close();
	/*****************************************************************************
	主程序控制段；
	*****************************************************************************/

	cout << "The results in FEMout.txt!" << endl << endl << "Program end." << endl << "****************************************" << endl << endl;
	long time02 = clock();
	cout << "It took " << double(time02 - time01) / CLOCKS_PER_SEC << "seconds." << endl << endl;
	system("pause");
}

//读取文件；
void read(string filename)
{
	ifstream File1;
	File1.open(filename, ios_base::in | ios_base::binary);

	//若文件打开失败需要提示；
	if (!File1)
	{
		cout << "read File open error!!!" << endl;
		system("pause"); exit(0);//直接终止程序运行，0表示正常退出；
	}

	//输出结点坐标；
	string KeyWord1;
	File1 >> KeyWord1 >> nodeNumber;
	for (int i = 0; i < nodeNumber; i++)
	{
		node.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 3; j++)
		{
			double temp;
			File1 >> temp;
			node[i].push_back(temp);
		}
	}

	//输出单元信息，适用四边形单元；
	string KeyWord2;
	File1 >> KeyWord2 >> elementNumber;
	for (int i = 0; i < elementNumber; i++)
	{
		element.push_back(vector<int>(0, 0));
		for (int j = 0; j < 6; j++)
		{
			double temp;
			File1 >> temp;
			element[i].push_back(temp);
		}
	}

	//输出材料信息；
	string KeyWord3;
	File1 >> KeyWord3 >> materialNumber;
	for (int i = 0; i < materialNumber; i++)
	{
		material.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 4; j++)
		{
			double temp;
			File1 >> temp;
			material[i].push_back(temp);
		}
	}

	//输出约束信息;
	string KeyWord4;
	File1 >> KeyWord4 >> boundNumber;
	for (int i = 0; i < boundNumber; i++)
	{
		bound.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 4; j++)
		{
			double temp;
			File1 >> temp;
			bound[i].push_back(temp);
		}
	}

	//输出分布荷载信息；
	string KeyWord5;
	File1 >> KeyWord4 >> faceloadNumber;
	for (int i = 0; i < faceloadNumber; i++)
	{
		faceload.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 7; j++)
		{
			double temp;
			File1 >> temp;
			faceload[i].push_back(temp);
		}
	}

	//输出集中荷载信息;
	string KeyWord6;
	File1 >> KeyWord6 >> nodeloadNumber;
	for (int i = 0; i < nodeNumber; i++)
	{
		nodeload.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 5; j++)
		{
			double temp;
			File1 >> temp;
			nodeload[i].push_back(temp);
		}
	}

	File1.close();//关闭文件；
}


// xiezhuoyu
// mechanics_xzy@163.com