#pragma once
#include <string>
#include <vector>
using namespace std;
/*********************************************************
	�������֧�䷽��[K]*{f}={R}��
	����LU�ֽ�������Դ��������飻
	LU�ֽⲽ����ҪΪ���ֽ�õ�L��U����ǰ�����ش���
	���յõ��Ľ��Ϊ���λ��ֵ��
**********************************************************/
extern int nodeNumber;

class solve_LU
{
public:
	solve_LU(string file);
	~solve_LU(void);

	//�õ����ڵ�λ��ֵ������getDisplacement()��
	//����	λ������displacement[];
	vector<double> getDisplacement();

private:
	string filename;

	//�ֽ�õ�L��U�����е�L���󣬺���LUglobalK()��
	//������������LU[]�����У�
	void doLU();

	//ǰ�����[L][g]=[R]�е�[g]������forward������
	void forward();

	//�ش����[U][displacement]=[g]������backwards������
	void backwards();

	vector<double> LU;//���ڴ洢LU�ֽ����ֻ�������ǵ�L����

	//���ڲ���LU������getLUij(�кţ��к�)��
	int getLUij(int row, int col);

	vector<double> displacement;

	//���λ��ֵ��FEMout.txt��
	void out();
};


// xiezhuoyu
// mechanics_xzy@163.com