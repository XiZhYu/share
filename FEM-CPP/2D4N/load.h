#pragma once
#include <string>
#include <vector>
using namespace std;
/************************************************************************
整体结点荷载列阵loadR的形成；
先后处理体力、集中力、面力，分别形成各自的等效结点荷载列阵loadRi；
组集成完整的loadR；
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
	//组集整体等效结点荷载列阵loadR，函数doLoadR()；
	void doLoadR();
	vector<double> loadR;//整体等效结点荷载列阵；

private:
	string filename;

	//处理体力形成的等效结点荷载列阵loadRi，函数doBodyforce()；
	void doBodyforce();

	//处理集中力形成的等效结点荷载列阵loadRi，函数donodeload()；
	void donodeload();

	//处理面力形成的等效结点荷载列阵loadRi，函数doFaceload()；
	void doFaceload();

	//由集中力编号得到相关的单元编号，函数nodeloadElement(集中力编号)；
	//返回	单元编号elementID；
	int nodeload_Element(int nodeloadID);

	//由全局坐标确定局部坐标，函数gloToLoc(全局坐标，所在单元编号)；
	//返回	局部坐标列阵rs[4]={r,s}；
	vector<double> gloToLoc(vector<vector<double>> xy, int elementID);

	//由分布力结点确定所在单元编号以及边编号，函数faceloadLine(分布力编号)；
	//返回	单元编号re[0]=elementID,边编号re[1]=lineID；
	vector<int> faceload_Line(int faceloadID);

};


// xiezhuoyu
// mechanics_xzy@163.com