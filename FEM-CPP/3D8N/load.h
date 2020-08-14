#pragma once
#include <string>
#include <vector>
#include "stiffness.h"
#include "matrix.h"
using namespace std;
/***********************************************************************************
���㼯�������ֲ�����������Ӧ�Ľ���������{nodeLoadR}��{faceLoadR}��{bodyLoadR}
����װ���������������{loadR}
************************************************************************************/
extern vector<vector<double>> node;
extern int nodeNumber;
extern vector<vector<int>> element;
extern int elementNumber;
extern vector<vector<double>> material;
extern vector<vector<double>> bound;
extern int boundNumber;
extern vector<vector<double>> faceload;
extern int faceloadNumber;
extern vector<vector<double>> nodeload;
extern int nodeloadNumber;

//�������������
extern vector<double> loadR;

class load
{
public:
	load();
	~load();

	//�鼯�������������loadR������doLoadR()��
	void getLoadR();
	
	//�κ�������N������(r,s,t)
	vector<vector<double>> N(double r, double s, double t);

	//���������γɵĽ���������{NLR}������nodeloadR()��
	void nodeloadR();

	//���������γɵĽ���������{FLR}������doFaceload()��
	void faceloadR();

	//���������γɵĽ���������{BLR}������doBodyforce()��
	void bodyloadR();

	//ȷ���ֲ����뵥Ԫ�Լ���Ĺ�ϵ��������faceloadID��
	//���طֲ������ڵĵ�ԪID����ID
	vector<int> faceload_ElementID_FaceID(int faceloadID);

	//����3*3����µĸ�˹�㴦�ĺ���ֵ
	vector<vector<double>> GaussPiontP(int faceloadID);

	//���㹫ʽ�е��������������
	//����{A}
	vector<vector<double>> A(int faceID, int i, int j, int elementID);

private:
	stiffness objS;
	matrix objM;
};

// xiezhuoyu
// mechanics_xzy@163.com