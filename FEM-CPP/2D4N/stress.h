#pragma once
#include <string>
#include <vector>
using namespace std;
/*********************************************************
	�ɹ�����
	�õ���Ԫ���ĵ�Ӧ��ֵ��
	[P]=[D][B][f]
**********************************************************/
extern int elementNumber;
extern int nodeNumber;
extern vector<vector<double>> node;
extern string solve_choice;
class stress
{
public:
	stress(string file);
	~stress(void);
	//���㵥Ԫ���ĵ㴦��Ӧ��ֵ������doStress������
	void doStress();

	//������Ӧ�ĵ�Ԫ�������飬����nodeLinkElement()
	//nodeLinkElement[i][j]Ϊi�Ž��ĵ�j����ص�Ԫ��ID
	vector<vector<int>> nodeLinkElement();

private:
	string filename;
	vector<vector<double>> nodeStress;//�ڵ�Ӧ��
	vector<vector<double>> p;//����ԪӦ��ֵ��x������Ӧ����y������Ӧ������Ӧ����

	vector<double> dis;

	//���Ӧ��ֵ��FEMout.txt��
	void out();

	//���VTK�ļ�
	void vtkOut();
};


// xiezhuoyu
// mechanics_xzy@163.com