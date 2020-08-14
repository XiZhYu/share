#include "dataInOut.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
using namespace std;

dataInOut::dataInOut()
{
}


dataInOut::~dataInOut()
{
}

//读取数据文件
void dataInOut::read(string filename)
{
	ifstream File1;
	File1.open(filename, ios_base::in | ios_base::binary);

	//若文件打开失败需要提示；
	if (!File1)
	{
		cout << "dataInOut::read File open error!!!" << endl;
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

	//输出单元信息，适用3D8N单元；
	string KeyWord2;
	File1 >> KeyWord2 >> elementNumber >> eleNodeMax;
	for (int i = 0; i < elementNumber; i++)
	{
		element.push_back(vector<int>(0, 0));
		for (int j = 0; j < 10; j++)
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
		for (int j = 0; j < 3; j++)
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
		for (int j = 0; j < 3; j++)
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
		for (int j = 0; j < 8; j++)
		{
			double temp;
			File1 >> temp;
			faceload[i].push_back(temp);
		}
	}

	//输出集中荷载信息;
	string KeyWord6;
	File1 >> KeyWord6 >> nodeloadNumber;
	for (int i = 0; i < nodeloadNumber; i++)
	{
		nodeload.push_back(vector<double>(0, 0.0));
		for (int j = 0; j < 4; j++)
		{
			double temp;
			File1 >> temp;
			nodeload[i].push_back(temp);
		}
	}

	File1.close();//关闭文件；
}


//输出位移值到FEMout.txt中
void dataInOut::nodeDisplacementOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//打开一个现有文件；
	if (!fileout) {
		cout << "dataInOut::nodeDisplacementOut File open error!!!" << endl;
		system("pause"); exit(0);
	}
	fileout << "NoodeID" << "\t\t" << "X-Displacement" << "\t\t" << "Y-Displacement" << "\t\t" << "Z-Displacement" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < nodeNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 3; j++)
		{
			fileout << displacement[3 * i + j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << resetiosflags(ios_base::scientific);
	fileout << "\r\n";
	fileout.close();
}


//输出单元中心应力值到FEMout.txt中
void dataInOut::centerStressOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//打开一个现存的文件；
	if (!fileout) {
		cout << "dataInOut::centerStressOut File open error!!!" << endl;
		system("pause"); exit(0);
	}
	fileout << "ElementID" << "\t" << "X-NormalStress" << "\t\t" << "Y-NormalStress" << "\t\t" << "Z-NormalStress" <<
		"\t\t" << "XY-ShearStress" << "\t\t" << "YZ-ShearStress" << "\t\t" << "ZX-ShearStress" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < elementNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 6; j++)
		{
			fileout << centerStress[i][j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << resetiosflags(ios_base::scientific);
	fileout << "\r\n";
	fileout.close();
}


//输出结点应力值到FEMout.txt中
void dataInOut::nodeStressOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//打开一个现存的文件；
	if (!fileout) {
		cout << "dataInOut::nodeStressOut File open error!!!" << endl;
		system("pause"); exit(0);
	}
	fileout << "NodeID" << "\t\t" << "X-NormalStress" << "\t\t" << "Y-NormalStress" << "\t\t" << "Z-NormalStress" <<
		"\t\t" << "XY-ShearStress" << "\t\t" << "YZ-ShearStress" << "\t\t" << "ZX-ShearStress" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < nodeNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 6; j++)
		{
			fileout << nodeStress[i][j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << resetiosflags(ios_base::scientific);
	fileout << "\r\n";
	fileout.close();
}

//输出VTK文件
void dataInOut::vtkOut()
{
	ofstream vtk;
	vtk.open("FEMoutVTK.vtk", ios_base::trunc|ios_base::binary);
	if (!vtk) 
	{
		cout << "dataInOut::vtkOut File open error!!!" << endl;
		system("pause"); exit(0);
	}

	vtk << "# vtk DataFile Version 2.0\r\n";// # 说明行（vtk版本）
	vtk << "fem_test\r\n";// 说明行（文件名）
	vtk << "ASCII\r\n"; // 说明行（编码系统）
	vtk << "DATASET UNSTRUCTURED_GRID\r\n";// 说明行（数据设置格式）


	//注：这里往下定义模型的几何形状，包括结点坐标跟单元结点编号
	vtk << "POINTS " << nodeNumber << " float\r\n";// POINTS：结点关键字  322：结点数  float：结点坐标的数据类型
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的三个坐标
		vtk << node[i][0] << " " << node[i][1] << " " << node[i][2] << " \r\n";
	}
	vtk << "CELLS " << elementNumber << " " << (9 * elementNumber) << "\r\n";// CELLS：单元关键字  283：单元数  1415：该部分数据总个数
	for (int i = 0; i < elementNumber; i++)
	{
		// 单元结点数目  该单元的结点编号
		vtk << element[i][0] << " " << element[i][2] << " " << element[i][3] << " " << element[i][4] << " " << element[i][5] << " " << element[i][6] << " " << element[i][7] << " " << element[i][8] << " " << element[i][9] << "\r\n";
	}
	vtk << "CELL_TYPES " << elementNumber << "\r\n";// CELL_TYPES：单元类型关键字  283：单元数
	for (int i = 0; i < elementNumber; i++)
	{
		//单元的单元类型编号，六面体为11
		vtk << 12 << "\r\n";
	}


	//注：以下数据是所要显示的各个变量，用户可根据自己需求添加
	vtk << "POINT_DATA " << nodeNumber << "\r\n";// POINT_DATA：关键字  322：变量数据的行数，一般为结点总数
	vtk << "SCALARS Coor(m) float 3\r\n";// SCALARS：标量数据关键字  Coor(m)：变量名  float：变量的数据类型  3：每一行数据个数
	vtk << "LOOKUP_TABLE default\r\n";// LOOKUP_TABLE：查表格式关键字  default：默认格式
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的三个坐标
		vtk << node[i][0] << " " << node[i][1] << " " << node[i][2] << " \r\n";
	}
	vtk << "VECTORS Displace float\r\n";// VECTORS：矢量数据关键字  Displace：变量名  float：变量的数据类型
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的三个位移
		vtk << displacement[3*i] << " " << displacement[3*i+1] << " " << displacement[3*i+2] << " \r\n";
	}
	vtk << "SCALARS Sx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的x方向正应力
		vtk << nodeStress[i][0]<<"\r\n";
	}
	vtk << "SCALARS Sy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的y方向正应力
		vtk << nodeStress[i][1] << "\r\n";
	}
	vtk << "SCALARS Sz float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的z方向正应力
		vtk << nodeStress[i][2] << "\r\n";
	}
	vtk << "SCALARS Txy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的xy方向切应力
		vtk << nodeStress[i][3] << "\r\n";
	}
	vtk << "SCALARS Tyz float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的yz方向切应力
		vtk << nodeStress[i][4] << "\r\n";
	}
	vtk << "SCALARS Tzx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的zx方向切应力
		vtk << nodeStress[i][5] << "\r\n";
	}
	vtk.close();
}

// xiezhuoyu
// mechanics_xzy@163.com