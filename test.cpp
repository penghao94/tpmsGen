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

int main() {
	const auto &f = [](double x, double y, double z)->double {
		return cos(x)+cos(y)+cos(z);
	};
	
	std::vector<double>Vmin = { 0,0,0 };
	std::vector<double>Vmax = { 2*PI,2*PI,2*PI };
	double levelset = 0.5;

	tpmsgen::TPMSGenerator<double,double> tg(f,Vmin, Vmax, 100);
	Eigen::MatrixXd grid, value;
	tg.getDistanceField(grid, value);
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	tg.makeLevelSet(levelset, V, F);
	writeOBJ("D:/project/tpmsGen/tpmsGen/data/p.obj", V, F);
	system("pause");
}
