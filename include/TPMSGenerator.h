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
		std::vector<DerivedSize>Vmin, Vmax;	// Bounding box
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
		template<typename DerivedG, typename DerivedV>
		void getDistanceField(Eigen::PlainObjectBase<DerivedG> &grid, Eigen::PlainObjectBase<DerivedV> &value);
		template<typename DerivedG, typename DerivedV>
		void getDistanceField(const size_t label, Eigen::PlainObjectBase<DerivedG> &grid, Eigen::PlainObjectBase<DerivedV> &value);
		template<typename DerivedVec, typename DerivedFct>
		void makeLevelSet(const DerivedValue &level, Eigen::PlainObjectBase<DerivedVec> &vertices, Eigen::PlainObjectBase<DerivedFct>&facets);
		template<typename DerivedVec, typename DerivedFct>
		void makeLevelSet(const size_t label, const DerivedValue &level, Eigen::PlainObjectBase<DerivedVec> &vertices, Eigen::PlainObjectBase<DerivedFct>&facets);
		template<typename DerivedVec, typename DerivedFct>
		void makeLevelSet(const size_t label, const DerivedValue &level, Eigen::PlainObjectBase<DerivedVec> &vertices, Eigen::PlainObjectBase<DerivedFct>&facets, DerivedValue &density);
		~TPMSGenerator() {};
	};



	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedG, typename DerivedV>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::getDistanceField(Eigen::PlainObjectBase <DerivedG>& grid, Eigen::PlainObjectBase <DerivedV>& value)
	{

		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution) + 1;


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

	// | 0/1 0/1 0/1 0/1 0/1 0/1|
	// | xf  xb  yf  yb  zf  zb |


	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedG, typename DerivedV>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::getDistanceField(const size_t label, Eigen::PlainObjectBase<DerivedG>& grid, Eigen::PlainObjectBase<DerivedV>& value)
	{
		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution) + 1;

		std::vector<int> bound(6);
		for (size_t i = 0; i < 6; i++)
			bound[i] = (label >> i) & 0x01;


		int res_0 = bound[0] + res[0] + bound[2];
		int res_1 = bound[1] + res[1] + bound[3];
		int res_2 = bound[4] + res[2] + bound[5];

		grid.resize(res_0 * res_1 * res_2, 3);
		value.resize(res_0 * res_1 * res_2, 1);

		const auto lerp = [&](const DerivedSize index, const size_t dim)->DerivedSize {
			return Vmin[dim] + index / static_cast<DerivedSize>(res[dim] - 1)
				*(Vmax[dim] - Vmin[dim]);
		};
		for (size_t z = 0; z < res_2; z++) {
			const DerivedSize pz = lerp(static_cast<DerivedSize>(z) - bound[4], 2);

			for (size_t y = 0; y < res_1; y++) {
				const DerivedSize py = lerp(static_cast<DerivedSize>(y) - bound[1], 1);

				for (size_t x = 0; x < res_0; x++) {
					const DerivedSize px = lerp(static_cast<DerivedSize>(x) - bound[0], 0);

					grid.row(x + res_0 * (y + res_1 * z)) << px, py, pz;
					value.row(x + res_0 * (y + res_1 * z)) << func(px, py, pz);
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
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution) + 1;

		tpmsgen::marching_cubes(value, grid, res[0], res[1], res[2], 0.0, vertices, facets);
	}

	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedVec, typename DerivedFct>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::makeLevelSet(const size_t label, const DerivedValue & level, Eigen::PlainObjectBase<DerivedVec>& vertices, Eigen::PlainObjectBase<DerivedFct>& facets)
	{
		Eigen::Matrix<DerivedSize, Eigen::Dynamic, 3> grid;
		Eigen::Matrix<DerivedValue, Eigen::Dynamic, 1> value;

		getDistanceField(label, grid, value);
		for (int i = 0; i < value.rows(); i++) value(i) -= level;
		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution) + 1;


		std::vector<size_t> bound(6);
		for (size_t i = 0; i < 6; i++)
			bound[i] = (label >> i) & 0x01;

		int res_0 = bound[0] + res[0] + bound[2];
		int res_1 = bound[1] + res[1] + bound[3];
		int res_2 = bound[4] + res[2] + bound[5];

		for (size_t z = 0; z < res_2; z++) {
			for (size_t y = 0; y < res_1; y++) {
				for (size_t x = 0; x < res_0; x++) {

					if (x == 0 && (label & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (x == res_0 - 1 && ((label >> 2) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (y == 0 && ((label >> 1) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (y == res_1 - 1 && ((label >> 3) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (z == 0 && ((label >> 4) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (z == res_2 - 1 && ((label >> 5) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);

				}
			}
		}
	

		tpmsgen::marching_cubes(value, grid, res_0, res_1, res_2, 0.0, vertices, facets);
		facets.col(1).swap(facets.col(2));
	}

	template<typename DerivedValue, typename DerivedSize>
	template<typename DerivedVec, typename DerivedFct>
	inline void TPMSGenerator<DerivedValue, DerivedSize>::makeLevelSet(const size_t label, const DerivedValue & level, Eigen::PlainObjectBase<DerivedVec>& vertices, Eigen::PlainObjectBase<DerivedFct>& facets, DerivedValue & density)
	{
		Eigen::Matrix<DerivedSize, Eigen::Dynamic, 3> grid;
		Eigen::Matrix<DerivedValue, Eigen::Dynamic, 1> value;

		getDistanceField(label, grid, value);
		for (int i = 0; i < value.rows(); i++) value(i) -= level;

		density = 0;
		for (int i = 0; i < value.rows(); i++)
			if (value(i) >= 0) density += 1;


		std::vector<size_t> res(3);
		for (int i = 0; i < 3; i++) res[i] = static_cast<size_t>((Vmax[i] - Vmin[i]) / PI * solution) + 1;


		std::vector<size_t> bound(6);
		for (size_t i = 0; i < 6; i++)
			bound[i] = (label >> i) & 0x01;

		int res_0 = bound[0] + res[0] + bound[2];
		int res_1 = bound[1] + res[1] + bound[3];
		int res_2 = bound[4] + res[2] + bound[5];
		density = density / res_0 / res_1 / res_2;

		for (size_t z = 0; z < res_2; z++) {
			for (size_t y = 0; y < res_1; y++) {
				for (size_t x = 0; x < res_0; x++) {

					if (x == 0 && (label & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (x == res_0 - 1 && ((label >> 2) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (y == 0 && ((label >> 1) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (y == res_1 - 1 && ((label >> 3) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (z == 0 && ((label >> 4) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);
					else if (z == res_2 - 1 && ((label >> 5) & 0x01))
						value(x + res_0 * (y + res_1 * z), 0) = value(x + res_0 * (y + res_1 * z), 0) > 0 ? 0 : value(x + res_0 * (y + res_1 * z), 0);

				}
			}
		}

		tpmsgen::marching_cubes(value, grid, res_0, res_1, res_2, 0.0, vertices, facets);
		facets.col(1).swap(facets.col(2));
	}

}


