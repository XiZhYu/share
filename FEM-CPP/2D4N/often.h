#pragma once
#include <string>
#include <vector>
using namespace std;
/*******************************************
���������о���ʹ��vector���д��档
*******************************************/
extern vector<vector<int>> element;
class often
{
public:
	often();//���캯����
	~often();//����������

	//����˷�������������A������B����
	vector<vector<double>> matrixM(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

	//����ת�ã�����������A����
	vector<vector<double>> matrixT(vector<vector<double>> matrixA);

	//�����ף�������棬����������A����
	vector<vector<double>> matrixN2(vector<vector<double>> matrixA);

	//������棬����������A����
	vector<vector<double>> matrixN(vector<vector<double>> matrixA);

	//�������ˣ�����������A����B����
	vector<vector<double>> matrixM2(vector<vector<double>> matrixA, double B);

	//������ӣ�����������A������B����
	vector<vector<double>> matrixP(vector<vector<double>> matrixA, vector<vector<double>>matrixB);

	//�������������������A������B����
	vector<vector<double>> matrixP2(vector<vector<double>>matrixA, vector<vector<double>>matrixB);

	//�����ף���������ʽ������������A����
	double matrixDet2(vector<vector<double>> matrixA);

	//�����ף���������ʽ������������A����
	double matrixDet3(vector<vector<double>> matrixA);

	//���������ʽ������������A����
	double matrixDet(vector<vector<double>> matrixA);

	//�ֲ��ڵ����ȫ�ֽڵ��ƥ�䣬��������Ԫ��ţ��ֲ�����Ӧ�����кţ������ļ�����
	//����										  ȫ�ֽ���Ӧ�����к�
	int match(int elementID, int i, string filename);

	//�κ�������N���������ֲ�����r���ֲ�����s����
	//����2*8��Ӧ�����B��
	vector<vector<double>> getN(double r, double s);

private:
	//�������ȡ������������A��row��col����
	vector<vector<double>> matrixSon(vector<vector<double>> matrixA, int row, int col);

};


// xiezhuoyu
// mechanics_xzy@163.com