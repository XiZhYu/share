#include <iostream>//���������
#include <fstream>//�ı����������
#include <string>//�ַ�����
#include <stdlib.h>//����exit();
#include <time.h> //��ʱ��
#include <vector>//���飻

#include "dataInOut.h"
#include "stiffness.h"
#include "load.h"
#include "solve_LU.h"
#include "solve_PCG.h"
#include "stress.h"

using namespace std;

/*****Դ����������************��һ�μ���õ��󣬲��ٱ��޸ģ����ǽ��еڶ��μ���**************************************************************************/
//�����Ϣ
vector<vector<double>> node;
int nodeNumber = 0;

//��Ԫ��Ϣ
vector<vector<int>> element;
int elementNumber = 0;
int eleNodeMax = 0;

//������Ϣ
vector<vector<double>> material;
int materialNumber = 0;

//Լ����Ϣ
vector<vector<double>> bound;
int boundNumber = 0;

//������Ϣ
vector<vector<double>> faceload;
int faceloadNumber = 0;

//��������Ϣ
vector<vector<double>> nodeload;
int nodeloadNumber = 0;

//����նȾ���COO�洢��ʽ
vector <double> globalStiffnessV;
vector <int> globalStiffnessR;
vector <int> globalStiffnessC;

//�������������
vector<double> loadR;

//������λ������
vector<double> displacement;

//����Ԫ���Ĵ�Ӧ��ֵ
vector<vector<double>> centerStress;
//�����Ӧ��ֵ
vector<vector<double>> nodeStress;

//���Է�����Ľⷨ��LU:LU�ֽ�-ֱ�ӷ���PCG:Ԥ�������ݶȷ�-��������
string solve_choice;
//PCG�Ľ��������Normal:�����������նȾ���EBE:�𲽵�Ԫ����
string PCG_choice;
/*****Դ����������***************************************************************************************************************************************/

void main()
{
	cout << "****************************************" << endl << "Program begin." << endl << endl;

	string filenamein, filenameout;
	//���������ļ���
	cout << "Enter the inputDataFile name:";	cin >> filenamein;

	//ѡ�����Դ���������ⷨ
	cout << endl << "Choose your method: PCG or LU" << "\t\t" << "Choice:";		cin >> solve_choice;
	if ("LU" == solve_choice)
	{
		cout << endl << "using " << solve_choice << " method" << endl << "running����" << endl << endl;
	}
	else
	{
		cout << "Choose your PCG_method: Normal(N) or Element-By-Element(EBE)" << "\t\t" << "PCG_Choice:";		cin >> PCG_choice;
		cout << endl << "using " << solve_choice << "-" << PCG_choice << " method" << endl << "running����" << endl << endl;
	}

	//��ʱ
	long time01 = clock();

	//��ȡ�����ļ�
	dataInOut objDI;	objDI.read(filenamein);
	//����һ������ļ�
	ofstream fileout;
	fileout.open("FEMout.txt", ios_base::out | ios_base::binary | ios_base::trunc);//���|������|ɾ��ԭ���ݣ�

	if ("EBE" != PCG_choice)//EBE�����������նȾ���
	{
		stiffness objS;	objS.globalStiffness();//�������նȾ���
	}
	load objL;	objL.getLoadR();//����������������
	if ("LU" == solve_choice)
	{
		solve_LU objSL;	objSL.getDisplacement();//��LU������Դ���������
	}
	else
	{
		solve_PCG objSP;	objSP.getDisplacement();//��PCG������Դ���������
	}
	stress objSt;	objSt.doStress();//��ⵥԪ���ġ���㴦��Ӧ��ֵ
	dataInOut objDO;	objDO.nodeDisplacementOut();	objDO.centerStressOut();	objDO.nodeStressOut();//���������
	objDO.vtkOut();

	fileout.close();
	cout << "The results in FEMout.txt!" << endl << endl << "Program end." << endl << "****************************************" << endl << endl;
	long time02 = clock();
	cout << "It took " << double(time02 - time01) / CLOCKS_PER_SEC << "seconds." << endl << endl;
	system("pause");
}


// xiezhuoyu
// mechanics_xzy@163.com