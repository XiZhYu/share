#pragma once
#include <string>
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/*********************************************************
�ɹ�����
�õ���Ԫ���ĺͽ�㴦��Ӧ��ֵ��
[P]=[D][B][f]
**********************************************************/
extern int elementNumber;
extern int nodeNumber;
extern string solve_choice;

extern vector<double> displacement;

//����Ԫ���Ĵ�Ӧ��ֵ��x��y��z������Ӧ����xy��yz��zx������Ӧ����
extern vector<vector<double>> centerStress;
//�����Ӧ��ֵ��x��y��z������Ӧ����xy��yz��zx������Ӧ����
extern vector<vector<double>> nodeStress;

class stress
{
public:
	stress();
	~stress();

	//����Ӧ��ֵ������doStress������
	void doStress();

	//�õ�����Ԫ���Ĵ�Ӧ��ֵ��getCenterStress();
	vector<vector<double>> getCenterStress();

	//�õ������Ӧ��ֵ��getNodeStress();
	vector<vector<double>> getNodeStress();

	//������Ӧ�ĵ�Ԫ�������飬����nodeLinkElement()
	//nodeLinkElement[i][j]Ϊi�Ž��ĵ�j����ص�Ԫ��ID
	void nodeLinkElement();

private:
	stiffness objS;
	matrix objM;

	vector<vector<int>> v_nodeLinkElement;
};


// xiezhuoyu
// mechanics_xzy@163.com