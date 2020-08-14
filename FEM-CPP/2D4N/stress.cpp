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


//���㵥Ԫ���ĵ㴦��Ӧ��ֵ������doStress������
void stress::doStress()
{
	often runoft;
	elementK runele(filename);
	/**********************************************************************/
	//ѡ��ⷨ��
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
		solve_PCG runPCG(filename);	 dis = runPCG.getDisplacement();//Ĭ�ϲ���PCG����
	}
	/**********************************************************************/
	p.resize(elementNumber, vector<double>(3, 0));
	for (int i = 0; i < elementNumber; i++)
	{
		//���ɵ�Ԫ�Ľ��λ������;
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

//������Ӧ�ĵ�Ԫ�������飬����nodeLinkElement()
//nodeLinkElement[i][j]Ϊi�Ž��ĵ�j����ص�Ԫ��ID
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



//���Ӧ��ֵ��FEMout.txt��
void stress::out()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ���ִ���ļ���
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


//���VTK�ļ�
void stress::vtkOut()
{
	ofstream vtk;
	vtk.open("FEMoutVTK.vtk", ios_base::trunc | ios_base::binary);
	if (!vtk)
	{
		cout << "stress::vtkOut File open error!!!" << endl;
		system("pause"); exit(0);
	}

	vtk << "# vtk DataFile Version 2.0\r\n";// # ˵���У�vtk�汾��
	vtk << "fem_test\r\n";// ˵���У��ļ�����
	vtk << "ASCII\r\n"; // ˵���У�����ϵͳ��
	vtk << "DATASET UNSTRUCTURED_GRID\r\n";// ˵���У��������ø�ʽ��


	//ע���������¶���ģ�͵ļ�����״����������������Ԫ�����
	vtk << "POINTS " << nodeNumber << " float\r\n";// POINTS�����ؼ���  322�������  float������������������
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ������������
		vtk << node[i][1] << " " << node[i][2] << " " << 0 << " \r\n";
	}
	vtk << "CELLS " << elementNumber << " " << (5 * elementNumber) << "\r\n";// CELLS����Ԫ�ؼ���  283����Ԫ��  1415���ò��������ܸ���
	for (int i = 0; i < elementNumber; i++)
	{
		// ��Ԫ�����Ŀ  �õ�Ԫ�Ľ����
		vtk << 4 << " " << element[i][2] << " " << element[i][3] << " " << element[i][4] << " " << element[i][5] << "\r\n";
	}
	vtk << "CELL_TYPES " << elementNumber << "\r\n";// CELL_TYPES����Ԫ���͹ؼ���  283����Ԫ��
	for (int i = 0; i < elementNumber; i++)
	{
		//��Ԫ�ĵ�Ԫ���ͱ�ţ�������Ϊ11
		vtk << 9 << "\r\n";
	}


	//ע��������������Ҫ��ʾ�ĸ����������û��ɸ����Լ��������
	vtk << "POINT_DATA " << nodeNumber << "\r\n";// POINT_DATA���ؼ���  322���������ݵ�������һ��Ϊ�������
	vtk << "SCALARS Coor(m) float 3\r\n";// SCALARS���������ݹؼ���  Coor(m)��������  float����������������  3��ÿһ�����ݸ���
	vtk << "LOOKUP_TABLE default\r\n";// LOOKUP_TABLE������ʽ�ؼ���  default��Ĭ�ϸ�ʽ
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ������������
		vtk << node[i][1] << " " << node[i][2] << " " << 0 << " \r\n";
	}
	vtk << "VECTORS Displace float\r\n";// VECTORS��ʸ�����ݹؼ���  Displace��������  float����������������
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ��������λ��
		vtk << dis[2 * i] << " " << dis[2 * i + 1]  << " " << 0 << " \r\n";
	}
	vtk << "SCALARS Sx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����x������Ӧ��
		vtk << nodeStress[i][0] << "\r\n";
	}
	vtk << "SCALARS Sy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����y������Ӧ��
		vtk << nodeStress[i][1] << "\r\n";
	}
	vtk << "SCALARS Txy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����xy������Ӧ��
		vtk << nodeStress[i][2] << "\r\n";
	}
	vtk.close();
}

stress::~stress(void)
{
}


// xiezhuoyu
// mechanics_xzy@163.com