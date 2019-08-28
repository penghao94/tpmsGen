#include<iostream>
#include<Eigen/core>
#include<functional>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cassert>
#include<TPMSGenerator.h>

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
	s <<
		V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n")) <<
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

int main() {
	const auto &f = [](double x, double y, double z)->double {
		return 10 * (cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x))
			- 0.5*(cos(2 * x)*cos(2 * y) + cos(2 * y)*cos(2 * z) + cos(2 * z)*cos(2 * x));
	};
	
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
	
}
