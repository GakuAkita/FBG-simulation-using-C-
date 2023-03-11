#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


using namespace std;
using namespace Eigen;

//���f�����g���������p����B
typedef std::complex<double> com;

//�萔�������ň�C�ɒ�`
//�P��
const double m = 1e-3;
const double micro = 1e-6;
const double nano = pow(0.1,9);
const double gai = pow(10, 20);
const double tera = pow(10, 12);

//���s���ܗ�
const double neff = 1.45431806;
//�~����
const double pi = 3.141592653589;
//����
const double c = 299792458;

//�\���̂̒�`�̂���������������낤���B

void Sum(vector<vector<double>> *array) {//���܂��s���Ă��邩�m�F���邽�߂ɍ�����B
	double sum = 0;
	for (int i = 0; i < (* array).at(0).size(); i++) {
		sum+=(* array).at(0).at(i);
	}
	cout << sum << endl;
}

void SumoneD(vector<double>* array) {//���܂��s���Ă��邩�m�F���邽�߂ɍ�����B
	double sum = 0;
	for (int i = 0; i < (*array).size(); i++) {
		sum += (*array).at(i);
	}
	cout << sum << endl;
}


//�P��g���ɔ��˗����v�Z���邽�߂̊֐��Bwave���g���Bmain�֐��ł��Ƃ͉񂹂Ηǂ�
double Refl(double* wave, vector<double>* sections, double* q) {
	//for���̊O�ɏo���Ȃ��Ƃ��߂Ȃ̂��B
	Matrix2cd R;
	//�ŏ��͒P�ʍs��
	R << com(1, 0), com(0, 0),
		com(0, 0), com(1, 0);
	for (int i = 0; i < (*sections).size(); i++) {//i�̉񐔂𒲐����Ă݂�B
		//�i�q�萔���؂茈��
		double dk = (*sections).at(i);

		//�����v�Z
		double delbeta = 2 * pi * neff / *wave - pi / dk;
		//s^2���v�Z
		double s2 = pow(*q, 2) - pow(delbeta, 2);

		//���`����s2<0���Ԉ���Ă���킽�Ԃ�Bs2>0�͂����Ƃł��Ă���͂��B
		//cout << s2 << endl;

		if (s2 > 0) {
			Matrix2cd Tk;
			Tk << com(cosh(sqrt(s2) * dk), -delbeta / sqrt(s2) * sinh(sqrt(s2) * dk)), com(*q / sqrt(s2) * sinh(sqrt(s2) * dk), 0),
				com(*q / sqrt(s2) * sinh(sqrt(s2) * dk), 0), com(cosh(sqrt(s2) * dk), delbeta / sqrt(s2) * sinh(sqrt(s2) * dk));
			R *= (Tk*Tk);
		}
		else {
			Matrix2cd Tk;
			Tk << com(cos(sqrt(-s2) * dk), -delbeta / sqrt(-s2) * sin(sqrt(-s2) * dk)), com(*q / sqrt(-s2) * sin(sqrt(-s2) * dk), 0),
				com(*q / sqrt(-s2) * sin(sqrt(-s2)*dk), 0), com(cos(sqrt(-s2)* dk), delbeta / sqrt(-s2) * sin(sqrt(-s2) * dk));
			R *= (Tk*Tk);
		}
	}

	double hansya = pow(abs(R(1, 0) / R(1, 1)), 2);
	return hansya;
}




//2�����z�����荞�ށB�R�[�h�̓��e�͂����Ɨ������Ă��Ȃ��B�ł�double�Ŏ�荞�ނł���B
//https://cplusplus.com/forum/general/233338/
vector<double> readCSV(string filename)//1�s�������ꍇ�B�����𒲐������2�����ɂ��ł���B
{
	vector<vector<double>> M;

	ifstream in(filename);
	string line;
	while (getline(in, line))                   // read a whole line of the file
	{
		stringstream ss(line);                     // put it in a stringstream (internal stream)
		vector<double> row;
		string data;
		while (getline(ss, data, ','))           // read (string) items up to a comma
		{
			row.push_back(stod(data));            // use stod() to convert to double; put in row vector
		}
		if (row.size() > 0) M.push_back(row);    // add non-empty rows to matrix
	}
	return M.at(0);
}


//2�����z���csv�ɏo��
//�Q�l������
//https://qiita.com/sakakibarakiaki/items/710aa9407411a06d1067
void write(vector<vector<double>>* data,const double *length,double *q,double *num) {
	cout << fixed << setprecision(25);//�͂����Ă���ɈӖ������邩�͂킩��Ȃ��B
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/���[�h�������_/CSVs/";// �����o����
	
	//csv�̖��O������B�Ȃ񂩂�����python�݂����Ɋy�ɕ�������Ȃ�����@�Ȃ����ȁH
	path.append("L");
	path.append(to_string((int)(*length * 1000)));//FBG�̒������w��
	path.append("_");
	path.append("q");
	path.append(to_string((int)*q));
	path.append("_");
	path.append("num");
	path.append(to_string((int) * num));
	path.append(".csv");

	std::ofstream ofs(path);// �����o���p�X�g���[��
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			for (size_t j = 0; j < (*data).at(i).size(); j++) ofs << (*data).at(i).at(j) << ",";
			ofs << std::endl;
		}
	}
}

void writeString(vector<vector<double>>* data) {
	//cout << fixed << setprecision(25);//�͂����Ă���ɈӖ������邩�͂킩��Ȃ��B
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/���[�h�������_/sample.csv";// �����o����
	std::ofstream ofs(path);// �����o���p�X�g���[��
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			for (size_t j = 0; j < (*data).at(i).size(); j++) ofs << (double)((*data).at(i).at(j)) << ",";
			ofs << std::endl;
		}
	}
}

void writeoneD(vector<double>* data) {
	//cout << fixed << setprecision(25);//�͂����Ă���ɈӖ������邩�͂킩��Ȃ��B
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/���[�h�������_/write1D.csv";// �����o����
	std::ofstream ofs(path);// �����o���p�X�g���[��
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			ofs << (double)((*data).at(i)) << ",";
		}
	}
}


//�ʑ��}�X�N�̍\����(����ς�߂�)

int main() {
	//cout << fixed << setprecision(25);

	//�ʑ��}�X�N�̒���
	const double L = 30 * m;
	//�ʑ��}�X�N�̕�
	const double d = 585.5 * nano;
	double dmin = 578.6 * nano;
	double dmax = 592.3 * nano;
	double interval = 0.25 * nano;
	//FBG�̃Z�N�V������
	const int num =L*2/d;

	//�g��
	double bragg = neff * d;

	//-----------------------------------------------
	//�x�N�g����FBG�̃Z�N�V���������`(Uniform�̂Ƃ��̂�)
	//vector<double> fbgs(num,d/2);
	//-----------------------------------------------

	//int n = 0;
	//double tempd = 0;
	//while (dmin+n*interval<dmax) {
	//	n++;
	//}

	//vector<double> dunique(n);
	//for (int i = 0; i < dunique.size(); i++) {
	//	dunique.at(i) = dmin + i * interval;
	//}

	////chirped�̏ꍇ�����B
	double divide = 31;//divide��ς����FBG�̒������ς��B

	//vector<double> Cfbgs;
	//�Ƃ肠���������߂Ĉ��python�ŏo�͂����ʑ��}�X�N�̃p�^�[�����������œǂݍ���ł�邱�Ƃɂ����B��
	vector<double> chirpedfbg = readCSV("C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/���[�h�������_/PMhalfsec_divide31half.csv");
	//���divide�ƈ�v������B

	//�g��
	double unitWL = 852.3 * nano;

	//�ŏ��g��
	double minwl = 839.5*nano;
	//�ő�g��
	double mxwl = 863.5*nano;
	//�����v�f��
	double xnum = 20;
	//����(�g��) 2�����̔z�������Ă����B
	vector<vector<double>> xywl(2,vector<double>(xnum));
	//�����ŉ������
	double delta = (mxwl - minwl) / xnum;
	//cout << delta << endl;
	for (int j = 0; j < xywl.at(0).size(); j++) {
		xywl.at(0).at(j) = minwl + j * delta;
	}
	//xywl.at(0)�ŉ���������B

	//���ܗ��ϒ���(���) ���f���ɂ��Ă��悢���ǁB
	double q0 = 3000;

	for (int j = 0; j < xywl.at(0).size(); j++) {
		xywl.at(1).at(j) = Refl(&xywl.at(0).at(j), &chirpedfbg, &q0);
	}

	////csv�ɏ����o��
	write(&xywl,&divide,&q0,&xnum);

	cout << "done" << endl;

	return 0;
}


//�v�f��1000��5���B
//�v�f��2000��5���ŏI����Ă��...�Ȃ�ł��B1000�̂Ƃ��̑���ԈႢ���A���邢�͎��R�ɑ����Ȃ��Ă���̂��H


