#include <iostream>//输入输出；
#include <fstream>//文本输入输出；
#include <string>//字符串；
#include <stdlib.h>//报错，exit();
#include <time.h> //计时；
#include <vector>//数组；

#include "dataInOut.h"
#include "stiffness.h"
#include "load.h"
#include "solve_LU.h"
#include "solve_PCG.h"
#include "stress.h"

using namespace std;

/*****源程序共享数据************第一次计算得到后，不再被修改，除非进行第二次计算**************************************************************************/
//结点信息
vector<vector<double>> node;
int nodeNumber = 0;

//单元信息
vector<vector<int>> element;
int elementNumber = 0;
int eleNodeMax = 0;

//材料信息
vector<vector<double>> material;
int materialNumber = 0;

//约束信息
vector<vector<double>> bound;
int boundNumber = 0;

//面力信息
vector<vector<double>> faceload;
int faceloadNumber = 0;

//集中力信息
vector<vector<double>> nodeload;
int nodeloadNumber = 0;

//整体刚度矩阵，COO存储格式
vector <double> globalStiffnessV;
vector <int> globalStiffnessR;
vector <int> globalStiffnessC;

//整体结点荷载列阵
vector<double> loadR;

//整体结点位移列阵
vector<double> displacement;

//各单元中心处应力值
vector<vector<double>> centerStress;
//各结点应力值
vector<vector<double>> nodeStress;

//线性方程组的解法，LU:LU分解-直接法，PCG:预处理共轭梯度法-迭代法；
string solve_choice;
//PCG的解决方法，Normal:正常组成整体刚度矩阵，EBE:逐步单元法；
string PCG_choice;
/*****源程序共享数据***************************************************************************************************************************************/

void main()
{
	cout << "****************************************" << endl << "Program begin." << endl << endl;

	string filenamein, filenameout;
	//输入数据文件名
	cout << "Enter the inputDataFile name:";	cin >> filenamein;

	//选择线性代数方程组解法
	cout << endl << "Choose your method: PCG or LU" << "\t\t" << "Choice:";		cin >> solve_choice;
	if ("LU" == solve_choice)
	{
		cout << endl << "using " << solve_choice << " method" << endl << "running……" << endl << endl;
	}
	else
	{
		cout << "Choose your PCG_method: Normal(N) or Element-By-Element(EBE)" << "\t\t" << "PCG_Choice:";		cin >> PCG_choice;
		cout << endl << "using " << solve_choice << "-" << PCG_choice << " method" << endl << "running……" << endl << endl;
	}

	//计时
	long time01 = clock();

	//读取数据文件
	dataInOut objDI;	objDI.read(filenamein);
	//创建一个输出文件
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::out | ios_base::binary | ios_base::trunc);//输出|二进制|删除原内容；

	if ("EBE" != PCG_choice)//EBE无需求解整体刚度矩阵
	{
		stiffness objS;	objS.globalStiffness();//求解整体刚度矩阵
	}
	load objL;	objL.getLoadR();//求解整体结点荷载列阵
	if ("LU" == solve_choice)
	{
		solve_LU objSL;	objSL.getDisplacement();//按LU求解线性代数方程组
	}
	else
	{
		solve_PCG objSP;	objSP.getDisplacement();//按PCG求解线性代数方程组
	}
	stress objSt;	objSt.doStress();//求解单元中心、结点处的应力值
	dataInOut objDO;	objDO.nodeDisplacementOut();	objDO.centerStressOut();	objDO.nodeStressOut();//输出计算结果
	objDO.vtkOut();

	fileout.close();
	cout << "The results in FEMout.txt!" << endl << endl << "Program end." << endl << "****************************************" << endl << endl;
	long time02 = clock();
	cout << "It took " << double(time02 - time01) / CLOCKS_PER_SEC << "seconds." << endl << endl;
	system("pause");
}


// xiezhuoyu
// mechanics_xzy@163.com