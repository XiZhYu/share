#pragma once
#include <string>
#include <vector>
using namespace std;
/*********************************************************
	成果整理。
	得到单元中心的应力值；
	[P]=[D][B][f]
**********************************************************/
extern int elementNumber;
extern int nodeNumber;
extern vector<vector<double>> node;
extern string solve_choice;
class stress
{
public:
	stress(string file);
	~stress(void);
	//计算单元中心点处的应力值，函数doStress（）；
	void doStress();

	//将结点对应的单元存入数组，函数nodeLinkElement()
	//nodeLinkElement[i][j]为i号结点的第j个相关单元的ID
	vector<vector<int>> nodeLinkElement();

private:
	string filename;
	vector<vector<double>> nodeStress;//节点应力
	vector<vector<double>> p;//各单元应力值：x方向正应力、y方向正应力、切应力；

	vector<double> dis;

	//输出应力值到FEMout.txt中
	void out();

	//输出VTK文件
	void vtkOut();
};


// xiezhuoyu
// mechanics_xzy@163.com