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

//��ȡ�����ļ�
void dataInOut::read(string filename)
{
	ifstream File1;
	File1.open(filename, ios_base::in | ios_base::binary);

	//���ļ���ʧ����Ҫ��ʾ��
	if (!File1)
	{
		cout << "dataInOut::read File open error!!!" << endl;
		system("pause"); exit(0);//ֱ����ֹ�������У�0��ʾ�����˳���
	}

	//���������ꣻ
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

	//�����Ԫ��Ϣ������3D8N��Ԫ��
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

	//���������Ϣ��
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

	//���Լ����Ϣ;
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

	//����ֲ�������Ϣ��
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

	//������к�����Ϣ;
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

	File1.close();//�ر��ļ���
}


//���λ��ֵ��FEMout.txt��
void dataInOut::nodeDisplacementOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ�������ļ���
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


//�����Ԫ����Ӧ��ֵ��FEMout.txt��
void dataInOut::centerStressOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ���ִ���ļ���
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


//������Ӧ��ֵ��FEMout.txt��
void dataInOut::nodeStressOut()
{
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::ate | ios_base::app);//��һ���ִ���ļ���
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

//���VTK�ļ�
void dataInOut::vtkOut()
{
	ofstream vtk;
	vtk.open("FEMoutVTK.vtk", ios_base::trunc|ios_base::binary);
	if (!vtk) 
	{
		cout << "dataInOut::vtkOut File open error!!!" << endl;
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
		vtk << node[i][0] << " " << node[i][1] << " " << node[i][2] << " \r\n";
	}
	vtk << "CELLS " << elementNumber << " " << (9 * elementNumber) << "\r\n";// CELLS����Ԫ�ؼ���  283����Ԫ��  1415���ò��������ܸ���
	for (int i = 0; i < elementNumber; i++)
	{
		// ��Ԫ�����Ŀ  �õ�Ԫ�Ľ����
		vtk << element[i][0] << " " << element[i][2] << " " << element[i][3] << " " << element[i][4] << " " << element[i][5] << " " << element[i][6] << " " << element[i][7] << " " << element[i][8] << " " << element[i][9] << "\r\n";
	}
	vtk << "CELL_TYPES " << elementNumber << "\r\n";// CELL_TYPES����Ԫ���͹ؼ���  283����Ԫ��
	for (int i = 0; i < elementNumber; i++)
	{
		//��Ԫ�ĵ�Ԫ���ͱ�ţ�������Ϊ11
		vtk << 12 << "\r\n";
	}


	//ע��������������Ҫ��ʾ�ĸ����������û��ɸ����Լ��������
	vtk << "POINT_DATA " << nodeNumber << "\r\n";// POINT_DATA���ؼ���  322���������ݵ�������һ��Ϊ�������
	vtk << "SCALARS Coor(m) float 3\r\n";// SCALARS���������ݹؼ���  Coor(m)��������  float����������������  3��ÿһ�����ݸ���
	vtk << "LOOKUP_TABLE default\r\n";// LOOKUP_TABLE������ʽ�ؼ���  default��Ĭ�ϸ�ʽ
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ������������
		vtk << node[i][0] << " " << node[i][1] << " " << node[i][2] << " \r\n";
	}
	vtk << "VECTORS Displace float\r\n";// VECTORS��ʸ�����ݹؼ���  Displace��������  float����������������
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ��������λ��
		vtk << displacement[3*i] << " " << displacement[3*i+1] << " " << displacement[3*i+2] << " \r\n";
	}
	vtk << "SCALARS Sx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����x������Ӧ��
		vtk << nodeStress[i][0]<<"\r\n";
	}
	vtk << "SCALARS Sy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����y������Ӧ��
		vtk << nodeStress[i][1] << "\r\n";
	}
	vtk << "SCALARS Sz float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����z������Ӧ��
		vtk << nodeStress[i][2] << "\r\n";
	}
	vtk << "SCALARS Txy float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����xy������Ӧ��
		vtk << nodeStress[i][3] << "\r\n";
	}
	vtk << "SCALARS Tyz float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����yz������Ӧ��
		vtk << nodeStress[i][4] << "\r\n";
	}
	vtk << "SCALARS Tzx float 1\r\n";
	vtk << "LOOKUP_TABLE default\r\n";
	for (int i = 0; i < nodeNumber; i++)
	{
		// ����Ϊ����zx������Ӧ��
		vtk << nodeStress[i][5] << "\r\n";
	}
	vtk.close();
}

// xiezhuoyu
// mechanics_xzy@163.com