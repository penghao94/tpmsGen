#pragma once
constexpr auto PI = 3.1415926;
#include<utility>
#include<functional>
#include<vector>
#include<Eigen/core>
#include "MarchingCubes.h"
namespace tpmsgen {

	template<typename DerivedValue, typename DerivedSize>
	class TPMSGenerator
	{
	public:
		std::function<DerivedValue(
			const DerivedSize&,
			const DerivedSize&,
			const DerivedSize&)> func;			//TPMS function
		std::vector<DerivedSize>Vmin,Vmax;	// Bounding box
		size_t solution;			// Resolution for each TPMS cell

		TPMSGenerator(const std::function<DerivedValue(
			const DerivedSize&,
			const DerivedSize&,
			const DerivedSize&)> f,
			const std::vector<DerivedSize>& m,
			const std::vector<DerivedSize>& M,
			const size_t s) {
			func = std::move(f); Vmin = m; Vmax = M; solution = s;
		};
		template<typename DerivedG,typename DerivedV>
		void getDistanceField(Eigen::PlainObjectBase<DerivedG> &grid, Eigen::PlainObjectBase<DerivedV> &value);
		
		template<typename DerivedVec,typename DerivedFct>
		void makeLevelSet(const DerivedValue &level, Eigen::PlainObjectBase<DerivedVec> &vertices, Eigen::PlainObjectBase<DerivedFct>&facets);

		~TPMSGenerator() {};
	};
	
	

	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedG, typename DerivedV>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::getDistanceField(Eigen::PlainObjectBase <DerivedG>& grid, Eigen::PlainObjectBase <DerivedV>& value)
	{
		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution);

		//std::cout << res[0] << res[1] << res[2] << std::endl;
		grid.resize(res[0] * res[1] * res[2], 3);
		value.resize(res[0] * res[1] * res[2], 1);

		const auto lerp = [&](const size_t index, const size_t dim)->DerivedSize {
			return Vmin[dim] + static_cast<DerivedSize>(index) / static_cast<DerivedSize>(res[dim] - 1)
				*(Vmax[dim] - Vmin[dim]);
		};

		for (size_t z = 0; z < res[2]; z++) {
			const DerivedSize pz = lerp(z, 2);
			for (size_t y = 0; y < res[1]; y++) {
				const DerivedSize py = lerp(y, 1);
				for (size_t x = 0; x < res[0]; x++) {
					const DerivedSize px = lerp(x, 0);

					grid.row(x + res[0] * (y + res[1] * z)) << px, py, pz;
					value.row(x + res[0] * (y + res[1] * z)) << func(px, py, pz);
				}
			}
		}
	}

	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedVec, typename DerivedFct>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::makeLevelSet(const DerivedValue & level,
		Eigen::PlainObjectBase<DerivedVec>& vertices, Eigen::PlainObjectBase<DerivedFct>& facets)
	{
		Eigen::Matrix<DerivedSize, Eigen::Dynamic, 3> grid;
		Eigen::Matrix<DerivedValue, Eigen::Dynamic, 1> value;

		getDistanceField(grid, value);
		for (int i = 0; i < value.rows(); i++) value(i) -= level;
		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution);

		tpmsgen::marching_cubes(value, grid, res[0], res[1], res[2],0.0, vertices, facets);
	}

}


