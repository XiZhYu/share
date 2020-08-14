#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "stress.h"
#include "often.h"
#include "solve_LU.h"
#include "solve_PCG.h"
#include "elementk.h"
using namespace std;

stress::stress(string file)
{
	filename = file;
}


//计算单元中心点处的应力值，函数doStress（）；
void stress::doStress()
{
	often runoft;
	elementK runele(filename);
	/**********************************************************************/
	//选择解法；
	if ("PCG" == solve_choice)
	{
		solve_PCG runPCG(filename);	 dis = runPCG.getDisplacement();
	}
	else if ("LU" == solve_choice)
	{
		solve_LU runLU(filename);	 dis = runLU.getDisplacement();
	}
	else
	{
		solve_PCG runPCG(filename);	 dis = runPCG.getDisplacement();//默认采用PCG法；
	}
	/**********************************************************************/
	p.resize(elementNumber, vector<double>(3, 0));
	for (int i = 0; i < elementNumber; i++)
	{
		//集成单元的结点位移列阵;
		vector<vector<double>> f{ { dis[runoft.match(i,0,filename)] },
									{ dis[runoft.match(i,1,filename)] },
									{ dis[runoft.match(i,2,filename)] },
									{ dis[runoft.match(i,3,filename)] },
									{ dis[runoft.match(i,4,filename)] },
									{ dis[runoft.match(i,5,filename)] },
									{ dis[runoft.match(i,6,filename)] },
									{ dis[runoft.match(i,7,filename)] } };
		vector<vector<double>> DBf = runoft.matrixM(
			runoft.matrixM(
				runele.getElaD(i)
				, runele.getStrB(0, 0, i))
			, f);
		for (int j = 0; j < 3; j++)
		{
			p[i][j] = DBf[j][0];
		}
	}

	nodeStress.clear();
	nodeStress.resize((nodeNumber), vector<double>(3, 0.0));
	often objM;
	for (int i = 0; i < (nodeNumber); i++)
	{
		vector<vector<double>> mean(1, vector<double>(3, 0.0));
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < nodeLinkElement()[i].size(); k++)
			{
				mean[0][j] += p[nodeLinkElement()[i][k]][j];
			}
		}
		mean = objM.matrixM2(mean, (1.0 / nodeLinkElement()[i].size()));
		nodeStress[i] = mean[0];
	}
	out();
	vtkOut();
}

//将结点对应的单元存入数组，函数nodeLinkElement()
//nodeLinkElement[i][j]为i号结点的第j个相关单元的ID
vector<vector<int>> stress::nodeLinkElement()
{
	vector<vector<int>> result(nodeNumber , vector<int>(0, 0));
	for (int i = 0; i < elementNumber; i++)
	{
		for (int j = 2; j < 6; j++)
		{
			result[element[i][j]].push_back(i);
		}
	}
	return result;
}



//输出应力值到FEMout.txt中
void stress::out()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//打开一个现存的文件；
	if (!fileout) {
		cout << "stress::out File open error!!!" << endl;
		system("pause"); exit(0);
	}
	fileout << "ElementID" << "\t" << "X-NormalStress" << "\t\t" << "Y-NormalStress" << "\t\t" << "XY-ShearStress" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < elementNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 3; j++)
		{
			fileout << p[i][j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << "\r\n";

	fileout << "NodeID" << "\t\t" << "X-NormalStress" << "\t\t" << "Y-NormalStress" << "\t\t" << "XY-ShearStress" << "\r\n";
	fileout << setiosflags(ios_base::scientific) << setprecision(3);
	for (int i = 0; i < nodeNumber; i++)
	{
		fileout << i << "\t\t";
		for (int j = 0; j < 3; j++)
		{
			fileout << nodeStress[i][j] << "\t\t";
		}
		fileout << "\r\n";
	}
	fileout << "\r\n";

	fileout << resetiosflags(ios_base::scientific);
	fileout << "\r\n";
	fileout.close();
}


//输出VTK文件
void stress::vtkOut()
{
	ofstream vtk;
	vtk.open("FEMoutVTK.vtk", ios_base::trunc | ios_base::binary);
	if (!vtk)
	{
		cout << "stress::vtkOut File open error!!!" << endl;
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
		vtk << node[i][1] << " " << node[i][2] << " " << 0 << " \r\n";
	}
	vtk << "CELLS " << elementNumber << " " << (5 * elementNumber) << "\r\n";// CELLS：单元关键字  283：单元数  1415：该部分数据总个数
	for (int i = 0; i < elementNumber; i++)
	{
		// 单元结点数目  该单元的结点编号
		vtk << 4 << " " << element[i][2] << " " << element[i][3] << " " << element[i][4] << " " << element[i][5] << "\r\n";
	}
	vtk << "CELL_TYPES " << elementNumber << "\r\n";// CELL_TYPES：单元类型关键字  283：单元数
	for (int i = 0; i < elementNumber; i++)
	{
		//单元的单元类型编号，六面体为11
		vtk << 9 << "\r\n";
	}


	//注：以下数据是所要显示的各个变量，用户可根据自己需求添加
	vtk << "POINT_DATA " << nodeNumber << "\r\n";// POINT_DATA：关键字  322：变量数据的行数，一般为结点总数
	vtk << "SCALARS Coor(m) float 3\r\n";// SCALARS：标量数据关键字  Coor(m)：变量名  float：变量的数据类型  3：每一行数据个数
	vtk << "LOOKUP_TABLE default\r\n";// LOOKUP_TABLE：查表格式关键字  default：默认格式
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的三个坐标
		vtk << node[i][1] << " " << node[i][2] << " " << 0 << " \r\n";
	}
	vtk << "VECTORS Displace float\r\n";// VECTORS：矢量数据关键字  Displace：变量名  float：变量的数据类型
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的三个位移
		vtk << dis[2 * i] << " " << dis[2 * i + 1]  << " " << 0 << " \r\n";
	}
	vtk << "SCALARS Sx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的x方向正应力
		vtk << nodeStress[i][0] << "\r\n";
	}
	vtk << "SCALARS Sy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的y方向正应力
		vtk << nodeStress[i][1] << "\r\n";
	}
	vtk << "SCALARS Txy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// 依次为结点的xy方向切应力
		vtk << nodeStress[i][2] << "\r\n";
	}
	vtk.close();
}

stress::~stress(void)
{
}


// xiezhuoyu
// mechanics_xzy@163.com