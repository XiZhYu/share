#pragma once
#include <string>
#include <vector>
#include "matrix.h"
using namespace std;
/**********************************************************
���㵥Ԫ�նȾ��󣬲��鼯������նȾ���
***********************************************************/
extern int elementNumber;
extern int boundNumber;
extern vector<vector<double>> bound;
extern vector<vector<int>> element;
extern vector<vector<double>> node;
extern vector<vector<double>> material;

//���徢�Ⱦ����0Ԫ��ֵ�������ǣ���
extern vector <double> globalStiffnessV;
//���徢�Ⱦ����0Ԫ�ص��кţ�
extern vector <int> globalStiffnessR;
//���徢�Ⱦ����0Ԫ�ص��кţ�
extern vector <int> globalStiffnessC;

class stiffness
{
public:
	stiffness();
	~stiffness();

	//�κ����Ծֲ������ƫ�������������ֲ�����r���ֲ�����s���ֲ�����t���ֲ�����t��
	vector<vector<vector<double>>> Nri(double r, double s, double t);

	//�Ÿ��Ⱦ���J���������ֲ�����r���ֲ�����s���ֲ�����t����Ԫ���ID����
	//����3*3���Ÿ��Ⱦ���J��
	vector<vector<double>> J(double r, double s, double t,  int elementID);

	//Ӧ�����B���������ֲ�����r���ֲ�����s���ֲ�����t����Ԫ���ID����
	//����6*24��Ӧ�����B��
	vector<vector<double>> B(double r, double s, double t, int elementID);

	//���Ծ���D����������Ԫ���ID����
	//����6*6�ĵ��Ծ���D��
	vector<vector<double>> D(int elementID);

	//���㵥Ԫ�նȾ���
	//����24*24�ĵ�Ԫ���Ⱦ���
	vector<vector<double>> elementStiffness(int elementID);

	//�鼯���徢�Ⱦ��󣬺���������
	//��COO��ʽ�洢���徢�Ⱦ���
	void globalStiffness();

	//�ֲ���������ȫ�ֽ�����ƥ��,��������Ԫ�ţ��ֲ������ţ���
	//����ȫ�ֵĽ����ţ�
	int match(int elementID, int i);

	//COO�洢��������Ԫ��ֵ���кţ��кţ���
	void COO(double value, int i, int j);

private:
	matrix objM;
};


// xiezhuoyu
// mechanics_xzy@163.com