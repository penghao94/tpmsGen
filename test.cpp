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
	/*const auto &f = [](double x, double y, double z)->double {
		return 10 * (cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x))
			- 0.5*(cos(2 * x)*cos(2 * y) + cos(2 * y)*cos(2 * z) + cos(2 * z)*cos(2 * x));*/

	const auto &f = [](double x, double y, double z)->double {
		return 10 * (cos(x) + cos(y) + cos(z))-5.1*(cos(x)*cos(y)+cos(y)*cos(z)+cos(z)*cos(x));

	};

	std::vector<double>Vmin = { 0,0,0 };
	std::vector<double>Vmax = { 2 * PI,2 * PI,2 * PI };

	/*double x = 2.5, y = 17.5, z =12.5 ;
	std::vector<double>Vmin = { (x - 10)*0.1*PI,(y - 10)*0.1*PI,(z - 10)*0.1*PI };
	std::vector<double>Vmax = { (x + 10)*0.1*PI,(y + 10)*0.1*PI,(z + 10)*0.1*PI };*/


	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	tpmsgen::TPMSGenerator<double, double> tg(f, Vmin, Vmax, 50);
	
	Eigen::MatrixXd ship(295, 3);
	double density;
	ship.row(0) << 15, 0, 0;
	for (int i = 0; i < 294; i++) {
		double levelset = 14.7 - 0.05*i;
		tg.makeLevelSet(0x3F, levelset, V, F, density);;
		//writeOBJ("../data/TEST/TP_level_" + std::to_string(levelset) + "_density_" + std::to_string(density) + "_X_" + std::to_string(1)
			//+ "_Y_" + std::to_string(1) + "_Z_" + std::to_string(1) + ".obj", V, F);
		ship.row(i+1) << levelset, density, sqrt(density / (3 * PI));
	}	
		//for (int j = 0; j < 3; j++) V.col(j) *= stretch[j] / (2 * PI);

		//writeOBJ("../data/TEST/TG_level_" + std::to_string(levelset) + "_density_" + std::to_string(density) + "_X_" + std::to_string(1)
		//	+ "_Y_" + std::to_string(1) + "_Z_" + std::to_string(1) + ".obj", V, F);
		////writeOBJ("../data/TEST/TG_level_" + std::to_string(i) + "_density_" + std::to_string(density) + "_X_" + std::to_string(1)
		//	//+ "_Y_" + std::to_string(1) + "_Z_" + std::to_string(1) + ".obj", V, F);
		////}
	writeLOG("../data/level_density_radii.csv", ship);
	system("pause");
	

}
