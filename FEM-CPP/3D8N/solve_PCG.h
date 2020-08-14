#pragma once
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/*******************************************************
Ԥ�������ݶȷ�precondition conjugate gradient��
����Jacobian Precondition��Ԥ�������ȡ��Ԫ�ضԽ���
********************************************************/
extern int nodeNumber;
extern int elementNumber;
extern vector<vector<int>> element;
extern int boundNumber;
extern vector<vector<double>> bound;
extern string PCG_choice;

extern vector <double> globalStiffnessV;
extern vector <int> globalStiffnessR;
extern vector <int> globalStiffnessC;
extern vector<double> loadR;

//������λ������
extern vector<double> displacement;

class solve_PCG
{
public:
	solve_PCG();
	~solve_PCG();

	//�õ����ڵ�λ��ֵ������getDisplacement()��
	//����	λ������displacement[];
	vector<double> getDisplacement();

	//PCGʵ�֣�
	void doPCG_Normal();
	void doPCG_EBE();

	//�õ�Ԥ������󣬺���������նȾ���ֵ������նȾ���ֵ�кţ�����նȾ���ֵ�кţ���
	//����Ԥ������󣬼����Խ�Ԫ�أ�
	vector<double> getMN_N();
	vector<double> getMN_EBE();

	//���ڼ���Z����
	vector<double> getZ(vector<double> M, vector<double> R);

	//���ڼ���Q����
	vector<double> getQ_N(vector<double> P);
	vector<double> getQ_EBE(vector<double> P);

private:
	matrix objM;
	stiffness objS;
};


// xiezhuoyu
// mechanics_xzy@163.com