#pragma once
#include <string>
#include <vector>
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
class solve_PCG
{
public:
	solve_PCG(string file);
	~solve_PCG();

	//�õ����ڵ�λ��ֵ������getDisplacement()��
	//����	λ������displacement[];
	vector<double> getDisplacement();

private:
	string filename;

	//PCGʵ�֣�
	void doPCG_Normal();
	void doPCG_EBE();

	//�õ�Ԥ������󣬺���������նȾ���ֵ������նȾ���ֵ�кţ�����նȾ���ֵ�кţ���
	//����Ԥ������󣬼����Խ�Ԫ�أ�
	vector<double> getMN(vector<double> globalKV, vector<int> globalKR, vector<int> globalKC);
	vector<double> getMN();

	//���ڼ���Z����
	vector<double> getZ(vector<double> M, vector<double> R);

	//���ڼ���Q����
	vector<double> getQ(vector<double>globalKV, vector<int>globalKR, vector<int>globalKC, vector<double> P);
	vector<double> getQ(vector<double> P);

	vector<double> displacement;


	//���λ��ֵ��FEMout.txt��
	void out();

};


// xiezhuoyu
// mechanics_xzy@163.com