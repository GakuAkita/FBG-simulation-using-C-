#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


using namespace std;
using namespace Eigen;

//複素数を使う時これを用いる。
typedef std::complex<double> com;

//定数をここで一気に定義
//単位
const double m = 1e-3;
const double micro = 1e-6;
const double nano = pow(0.1,9);
const double gai = pow(10, 20);
const double tera = pow(10, 12);

//実行屈折率
const double neff = 1.45431806;
//円周率
const double pi = 3.141592653589;
//光速
const double c = 299792458;

//構造体の定義のしかたをもう一回やろうか。

void Sum(vector<vector<double>> *array) {//うまく行っているか確認するために作った。
	double sum = 0;
	for (int i = 0; i < (* array).at(0).size(); i++) {
		sum+=(* array).at(0).at(i);
	}
	cout << sum << endl;
}

void SumoneD(vector<double>* array) {//うまく行っているか確認するために作った。
	double sum = 0;
	for (int i = 0; i < (*array).size(); i++) {
		sum += (*array).at(i);
	}
	cout << sum << endl;
}


//単一波長に反射率を計算するための関数。waveが波長。main関数であとは回せば良い
double Refl(double* wave, vector<double>* sections, double* q) {
	//for文の外に出さないとだめなのか。
	Matrix2cd R;
	//最初は単位行列
	R << com(1, 0), com(0, 0),
		com(0, 0), com(1, 0);
	for (int i = 0; i < (*sections).size(); i++) {//iの回数を調整してみる。
		//格子定数を借り決め
		double dk = (*sections).at(i);

		//δβを計算
		double delbeta = 2 * pi * neff / *wave - pi / dk;
		//s^2を計算
		double s2 = pow(*q, 2) - pow(delbeta, 2);

		//あ～これs2<0が間違っているわたぶん。s2>0はちゃんとできているはず。
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




//2次元配列を取り込む。コードの内容はちゃんと理解していない。でもdoubleで取り込むできる。
//https://cplusplus.com/forum/general/233338/
vector<double> readCSV(string filename)//1行だけ取る場合。引数を調整すれば2次元にもできる。
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


//2次元配列をcsvに出力
//参考文献↓
//https://qiita.com/sakakibarakiaki/items/710aa9407411a06d1067
void write(vector<vector<double>>* data,const double *length,double *q,double *num) {
	cout << fixed << setprecision(25);//はたしてこれに意味があるかはわからない。
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/モード結合理論/CSVs/";// 書き出し先
	
	//csvの名前を決定。なんかもっとpythonみたいに楽に文字列をつなげる方法ないかな？
	path.append("L");
	path.append(to_string((int)(*length * 1000)));//FBGの長さを指定
	path.append("_");
	path.append("q");
	path.append(to_string((int)*q));
	path.append("_");
	path.append("num");
	path.append(to_string((int) * num));
	path.append(".csv");

	std::ofstream ofs(path);// 書き出し用ストリーム
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			for (size_t j = 0; j < (*data).at(i).size(); j++) ofs << (*data).at(i).at(j) << ",";
			ofs << std::endl;
		}
	}
}

void writeString(vector<vector<double>>* data) {
	//cout << fixed << setprecision(25);//はたしてこれに意味があるかはわからない。
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/モード結合理論/sample.csv";// 書き出し先
	std::ofstream ofs(path);// 書き出し用ストリーム
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			for (size_t j = 0; j < (*data).at(i).size(); j++) ofs << (double)((*data).at(i).at(j)) << ",";
			ofs << std::endl;
		}
	}
}

void writeoneD(vector<double>* data) {
	//cout << fixed << setprecision(25);//はたしてこれに意味があるかはわからない。
	std::string path = "C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/モード結合理論/write1D.csv";// 書き出し先
	std::ofstream ofs(path);// 書き出し用ストリーム
	if (ofs) {
		for (size_t i = 0; i < (*data).size(); i++) {
			ofs << (double)((*data).at(i)) << ",";
		}
	}
}


//位相マスクの構造体(やっぱやめる)

int main() {
	//cout << fixed << setprecision(25);

	//位相マスクの長さ
	const double L = 30 * m;
	//位相マスクの幅
	const double d = 585.5 * nano;
	double dmin = 578.6 * nano;
	double dmax = 592.3 * nano;
	double interval = 0.25 * nano;
	//FBGのセクション数
	const int num =L*2/d;

	//波長
	double bragg = neff * d;

	//-----------------------------------------------
	//ベクトルでFBGのセクション幅を定義(Uniformのときのみ)
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

	////chirpedの場合を作る。
	double divide = 31;//divideを変えればFBGの長さが変わる。

	//vector<double> Cfbgs;
	//とりあえず一回諦めて一回pythonで出力した位相マスクのパターンをこっちで読み込んでやることにした。↓
	vector<double> chirpedfbg = readCSV("C:/Users/gakum/Box/Aokilab/Team_Folders/FBG2021/Python/モード結合理論/PMhalfsec_divide31half.csv");
	//上のdivideと一致させる。

	//波長
	double unitWL = 852.3 * nano;

	//最小波長
	double minwl = 839.5*nano;
	//最大波長
	double mxwl = 863.5*nano;
	//横軸要素数
	double xnum = 2000;
	//横軸(波長) 2次元の配列を作っていく。
	vector<vector<double>> xywl(2,vector<double>(xnum));
	//ここで横軸作る
	double delta = (mxwl - minwl) / xnum;
	//cout << delta << endl;
	for (int j = 0; j < xywl.at(0).size(); j++) {
		xywl.at(0).at(j) = minwl + j * delta;
	}
	//xywl.at(0)で横軸が取れる。

	//屈折率変調量(一定) 複素数にしてもよいけど。
	double q0 = 3000;

	for (int j = 0; j < xywl.at(0).size(); j++) {
		xywl.at(1).at(j) = Refl(&xywl.at(0).at(j), &chirpedfbg, &q0);
	}

	////csvに書き出し
	write(&xywl,&divide,&q0,&xnum);

	cout << "done" << endl;

	return 0;
}


//要素数1000で5分。
//要素数2000で5分で終わってるな...なんでだ。1000のときの測り間違いか、あるいは自然に速くなっているのか？


