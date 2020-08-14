#pragma once
#include <vector>
using namespace std;
/*********************************************************
�������֧�䷽��[K]*{f}={R}��
����LU�ֽ�������Դ��������飻
LU�ֽⲽ����ҪΪ���ֽ�õ�L��U����ǰ�����ش���
���յõ��Ľ��Ϊ���λ��ֵ��
**********************************************************/
extern int nodeNumber;

extern vector <double> globalStiffnessV;
extern vector <int> globalStiffnessR;
extern vector <int> globalStiffnessC;
extern vector<double> loadR;

//������λ������
extern vector<double> displacement;

class solve_LU
{
public:
	solve_LU();
	~solve_LU(void);

	//�õ����ڵ�λ��ֵ������getDisplacement()��
	//����	λ������displacement[];
	vector<double> getDisplacement();

	//�ֽ�õ�L��U�����е�L���󣬺���doLU()��
	//������������LU[]�����У�
	void doLU();

	//ǰ�����[L][g]=[R]�е�[g]������forward������
	void forward();

	//�ش����[U][displacement]=[g]������backwards������
	void backwards();

	//���ڲ���LU������getLUij(�кţ��к�)��
	int getLUij(int row, int col);

private:
	vector<double> LU;//���ڴ洢LU�ֽ����ֻ�������ǵ�L����
	
};


// xiezhuoyu
// mechanics_xzy@163.com