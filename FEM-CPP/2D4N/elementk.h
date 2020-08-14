#pragma once
#include <string>
#include <vector>
using namespace std;
/**********************************************************
	��������նȾ���
***********************************************************/
extern int elementNumber;
extern int boundNumber;
extern vector<vector<double>> bound;
extern vector<vector<int>> element;
extern vector<vector<double>> node;
extern vector<vector<double>> material;
class elementK
{
public:
	elementK(string file);//���캯����
	~elementK();

	//���㵥Ԫ�նȾ���
	//����8*8�ĵ�Ԫ���Ⱦ���
	vector<vector<double>> getElementK(int elementID);

	//�鼯���徢�Ⱦ��󣬺���������
	//��COO��ʽ�洢���徢�Ⱦ���
	void doGlobalK();//createGlobalK��

	//Ӧ�����B���������ֲ�����r���ֲ�����s����Ԫ���ID����
	//����3*8��Ӧ�����B��
	vector<vector<double>> getStrB(double r, double s, int elementID);//getStrainB

	//��ƽ��Ӧ�����⣩���Ծ���D����������Ԫ���ID����
	//����3*3�ĵ��Ծ���D��
	vector<vector<double>> getElaD(int elementID);//getElasticityD

	//�Ÿ��Ⱦ���J���������ֲ�����r���ֲ�����s����Ԫ���ID����
	//����2*2���Ÿ��Ⱦ���J��
	vector<vector<double>> getJacJ(double r, double s, int elementID);

	//���徢�Ⱦ����0Ԫ��ֵ�������ǣ���
	vector <double> globalKV;
	//���徢�Ⱦ����0Ԫ�ص��кţ�
	vector <int> globalKR;
	//���徢�Ⱦ����0Ԫ�ص��кţ�
	vector <int> globalKC;

private:
	string filename;//�����ļ�����

	//COO�洢��������Ԫ��ֵ���кţ��кţ���
	void useCOO(double value, int i, int j);
};


// xiezhuoyu
// mechanics_xzy@163.com