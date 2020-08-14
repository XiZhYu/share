#pragma once
#include <string>
#include <vector>
using namespace std;
/************************************************************************
�������������loadR���γɣ�
�Ⱥ������������������������ֱ��γɸ��Եĵ�Ч����������loadRi��
�鼯��������loadR��
*************************************************************************/
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
class load
{
public:
	load(string file);
	~load(void);
	//�鼯�����Ч����������loadR������doLoadR()��
	void doLoadR();
	vector<double> loadR;//�����Ч����������

private:
	string filename;

	//���������γɵĵ�Ч����������loadRi������doBodyforce()��
	void doBodyforce();

	//���������γɵĵ�Ч����������loadRi������donodeload()��
	void donodeload();

	//���������γɵĵ�Ч����������loadRi������doFaceload()��
	void doFaceload();

	//�ɼ�������ŵõ���صĵ�Ԫ��ţ�����nodeloadElement(���������)��
	//����	��Ԫ���elementID��
	int nodeload_Element(int nodeloadID);

	//��ȫ������ȷ���ֲ����꣬����gloToLoc(ȫ�����꣬���ڵ�Ԫ���)��
	//����	�ֲ���������rs[4]={r,s}��
	vector<double> gloToLoc(vector<vector<double>> xy, int elementID);

	//�ɷֲ������ȷ�����ڵ�Ԫ����Լ��߱�ţ�����faceloadLine(�ֲ������)��
	//����	��Ԫ���re[0]=elementID,�߱��re[1]=lineID��
	vector<int> faceload_Line(int faceloadID);

};


// xiezhuoyu
// mechanics_xzy@163.com