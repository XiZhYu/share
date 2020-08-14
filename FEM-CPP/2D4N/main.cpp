#include <iostream>//���������
#include <fstream>//�ı����������
#include <string>//�ַ�����
#include <stdlib.h>//����exit();
#include <time.h> //��ʱ��
#include <vector>//���飻
#include "stress.h"
using namespace std;

/*****************************************************************************
��ȡ�����ļ���
	�����Ϣ�����ꣻ
	��Ԫ��Ϣ������ID����ؽ��ID��
	������Ϣ������ģ�������ɱȡ��ضȣ�
	Լ����Ϣ����ؽ��ID�����λ�Ʒ���λ�ƴ�С��
	������Ϣ����ؽ��ID��������С��
	��������Ϣ�������õ�λ�á�����С��
*****************************************************************************/
void read(string filename);

/*****************************************************************************
�����Ϣ������š�����,�������
��ȡ�ڵ�����node[i][j],iΪ����ţ�jΪ0��1��2�ֱ��Ӧ����š�x���ꡢy����;
��Ӧ��ʽ��
			node	2
			0	0	0
			1	0	1
*****************************************************************************/
vector<vector<double>> node;
int nodeNumber = 0;//�����Ŀ��

/*****************************************************************************
��Ԫ��Ϣ����Ԫ��š�����ID����ؽ��ID,��Ԫ����
��ȡ��Ԫ��Ϣelement[i][j]��iΪ��Ԫ��ţ�jΪ0��5�ֱ��Ӧ��Ԫ��š����ϡ�����ĸ��ڵ���;
��Ӧ��ʽ��
			ELEMENT	2
			0	0	0	3	1   2
			1	0	1	3	2   4
*****************************************************************************/
vector<vector<int>> element;
int elementNumber = 0;//��Ԫ��Ŀ��

/*****************************************************************************
������Ϣ�����ϱ�š�����ģ�������ɱȡ��ض�,��������
��ȡ������Ϣ,material[i][j];iΪ���ϱ�ţ�jΪ0��3�ֱ��Ӧ���ϱ�š�����ģ�������ɱȡ��ضȣ�
��Ӧ��ʽ��
			MATERIAL	1
			0	200000	0.3	77
*****************************************************************************/
vector<vector<double>> material;
int materialNumber = 0;//������Ŀ��

/*****************************************************************************
Լ����Ϣ��Լ����š���ؽ��ID�����λ�Ʒ���λ�ƴ�С��Լ������
��ȡԼ����ϢBound[i][j]��iΪԼ����ţ�jΪ0��3�ֱ��ӦԼ����š�����š�xy����Լ��λ��ֵ;
��Ӧ��ʽ��
			BOUND	2
			0	0	1	0
			1	0	2	0
*****************************************************************************/
vector<vector<double>> bound;
int boundNumber = 0;//Լ����Ŀ��

/*****************************************************************************
�ֲ�����Ϣ���ֲ�����š���ؽ��ID�����÷��򡢺��ش�С���ֲ���������
��ȡ�ֲ�������Ϣfaceload[i][j]��iΪ�ֲ�����ţ�jΪ0��6�ֱ�Ϊ���1ID�����2ID�����1��ˮƽ����ֵ�����2��ˮƽ����ֵ�����1����ֱ����ֵ�����2����ֱ����ֵ��
��Ӧ��ʽ��
			Faceload	1
			0	1	2	0	0	10	20
*****************************************************************************/
vector<vector<double>> faceload;
int faceloadNumber = 0;//������Ŀ��

/*****************************************************************************
��������Ϣ����������š�Լ����š���x������y������x��������y������
��ȡ���к�����Ϣnodeload[i][j]��iΪԼ����ţ�jΪ0��4�ֱ�������š����õ�x���ꡢy���ꡢ��x��������y������
��Ӧ��ʽ��
			nodeForce	2
			0		3.2		4.8		0		10
			1		1.5		9.7		10		0
*****************************************************************************/
vector<vector<double>> nodeload;
int nodeloadNumber = 0;//��������Ŀ��

//���Է�����Ľⷨ��LU:LU�ֽ�-ֱ�ӷ���PCG:Ԥ�������ݶȷ�-��������
string solve_choice;
//PCG�Ľ��������Normal:�����������նȾ���EBE:�𲽵�Ԫ����
string PCG_choice;

void main()
{
	cout << "****************************************" << endl << "Program begin." << endl << endl;

	/*****************************************************************************
	��������ƶΣ�
	*****************************************************************************/
	string filenamein, filenameout;
	cout << "Enter the inputDataFile name:";	cin >> filenamein;
	cout << endl << "Choose your method: PCG or LU" << "\t\t" << "Choice:";		cin >> solve_choice;
	if ("LU" != solve_choice)
	{
		cout << "Choose your PCG_method: Normal(N) or Element-By-Element(EBE)" << "\t\t" << "PCG_Choice:";		cin >> PCG_choice;
	}
	cout << endl << "running����" << endl << endl;
	long time01 = clock();  //��ʱ�� 
	//��ȡ�����ļ���
	read(filenamein);
	//����һ������ļ���
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::out | ios_base::binary | ios_base::trunc);//���|������|ɾ��ԭ���ݣ�
	//���㵥Ԫ���Ĵ�Ӧ����
	stress run6(filenamein);	run6.doStress();
	fileout.close();
	/*****************************************************************************
	��������ƶΣ�
	*****************************************************************************/

	cout << "The results in FEMout.txt!" << endl << endl << "Program end." << endl << "****************************************" << endl << endl;
	long time02 = clock();
	cout << "It took " << double(time02 - time01) / CLOCKS_PER_SEC << "seconds." << endl << endl;
	system("pause");
}

//��ȡ�ļ���
void read(string filename)
{
	ifstream File1;
	File1.open(filename, ios_base::in | ios_base::binary);

	//���ļ���ʧ����Ҫ��ʾ��
	if (!File1)
	{
		cout << "read File open error!!!" << endl;
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

	//�����Ԫ��Ϣ�������ı��ε�Ԫ��
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

	//���������Ϣ��
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

	//���Լ����Ϣ;
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

	//����ֲ�������Ϣ��
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

	//������к�����Ϣ;
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

	File1.close();//�ر��ļ���
}


// xiezhuoyu
// mechanics_xzy@163.com