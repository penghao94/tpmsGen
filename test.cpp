#include<iostream>
#include<omp.h>
#include<Eigen/core>
#include<Eigen/sparse>
#include<functional>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <Eigen/Geometry>
#include<TPMSGenerator.h>
#include<queue>
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>

using namespace std;
using namespace Eigen;
bool writeOBJ(
	const std::string str,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F)
{
	using namespace std;
	using namespace Eigen;
	assert(V.cols() == 3 && "V should have 3 columns");
	ofstream s(str);
	if (!s.is_open())
	{
		fprintf(stderr, "IOError: writeOBJ() could not open %s\n", str.c_str());
		return false;
	}
	s <<V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n")) <<
		(F.array() + 1).format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "f ", "", "", "\n"));
	return true;
}
bool writeLOG(const std::string & str, const Eigen::MatrixXd & param, int index=0)
{
	using namespace std;
	using namespace Eigen;

	ofstream csv;

	if (index == 0)
		csv.open(str);
	else
		csv.open(str, ios::app);

	if (!csv.is_open())
	{
		fprintf(stderr, "IOError: writeCSV() could not open %s\n", str.c_str());
		return false;
	}

	csv << param.format(IOFormat(FullPrecision, DontAlignCols, ",", "\n", "", "", "", "\n"));
	csv.close();
}


void func(const Eigen::VectorXd &data) {
		Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(data.data()), 5);
		std::cout << x << std::endl;
}
int main() {
	const auto &f = [](double x, double y, double z)->double {
		return 10 * (cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x))
			- 0.5*(cos(2 * x)*cos(2 * y) + cos(2 * y)*cos(2 * z) + cos(2 * z)*cos(2 * x));
<<<<<<< HEAD
=======
		//return cos(x) + cos(y) + cos(z);
>>>>>>> 5484bc40cf5bddf18dfaed605910d67411e8540d
	};

	//std::vector<double>Vmin = { 0,0,0 };
	//std::vector<double>Vmax = { 2 * PI,2 * PI,2 * PI };

	double x = 2.5, y = 17.5, z =12.5 ;
	std::vector<double>Vmin = { (x - 10)*0.1*PI,(y - 10)*0.1*PI,(z - 10)*0.1*PI };
	std::vector<double>Vmax = { (x + 10)*0.1*PI,(y + 10)*0.1*PI,(z + 10)*0.1*PI };


	double levelset = 9;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	tpmsgen::TPMSGenerator<double, double> tg(f, Vmin, Vmax, 15);
	tg.makeLevelSet(0xFF, levelset, V, F);

	writeOBJ("D:/OneDrive/data/2019-07-01//cell" + std::to_string(x) + "_" + std::to_string(y) + "_" + std::to_string(z) + ".obj", V, F);
	Eigen::SparseMatrix<int> G(4,4);
	G.insert(1, 1) = 9;
	G.insert(2, 2) = 3;
	G.insert(3, 3) = 1;

	std::cout << G << std::endl;

	G.coeffRef(2, 2) = 0;
	G.makeCompressed();
	G.
	for (int k = 0; k < G.outerSize(); ++k)
		for (SparseMatrix<int>::InnerIterator it(G, k); it; ++it)
		{
			std::cout << it.value() << std::endl;; // 元素值
			it.row();   // 行标row index
			it.col();   // 列标（此处等于k）
			it.index(); // 内部索引，此处等于it.row()
		}
	
<<<<<<< HEAD
	double x = 12.5, y = 2.5, z = 17.5;
	std::vector<double>Vmin = { (x - 10)*0.1*PI,(y - 10)*0.1*PI,(z - 10)*0.1*PI };
	std::vector<double>Vmax = { (x + 10)*0.1*PI,(y + 10)*0.1*PI,(z + 10)*0.1*PI };
	double levelset = 5;

	tpmsgen::TPMSGenerator<double,double> tg(f,Vmin, Vmax, 100);
	Eigen::MatrixXd C2W(290, 2);
	for (int i = 0; i <290; i++) {
	std::vector<double> h = {1};
	//for (int k = 0; k < h.size(); k++) {
		std::vector<double> stretch = { 10 * h[0],10 * h[0],10*h[0] };
		double density = 0;
		Eigen::MatrixXd V; Eigen::MatrixXi F;
		tg.makeLevelSet(0x3F, static_cast<double>(14.6-0.1*i), V, F, density);
		C2W.row(i) << static_cast<double>(14.6 - 0.1*i),density;
		std::cout << C2W.row(i) << std::endl;
		//for (int j = 0; j < 3; j++) V.col(j) *= stretch[j] / (2 * PI);

		/*writeOBJ("../data/TEST/TG_level_" + std::to_string(i) + "_density_" + std::to_string(density) + "_X_" + std::to_string(stretch[0])
			+ "_Y_" + std::to_string(stretch[1]) + "_Z_" + std::to_string(stretch[2]) + ".obj", V, F);*/
		//writeOBJ("../data/TEST/TG_level_" + std::to_string(i) + "_density_" + std::to_string(density) + "_X_" + std::to_string(1)
			//+ "_Y_" + std::to_string(1) + "_Z_" + std::to_string(1) + ".obj", V, F);
		//}
		
	}
	writeLOG("../data/levelset2density.csv", C2W);
	system("pause");
	
=======

>>>>>>> 5484bc40cf5bddf18dfaed605910d67411e8540d
}
