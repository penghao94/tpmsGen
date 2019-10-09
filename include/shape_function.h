#pragma once
#include <Eigen/core>

namespace tpmsgen {

	template<typename DerivedV>
	void ShapeFunction(
		const Eigen::PlainObjectBase<DerivedV> &before,
		const Eigen::PlainObjectBase<DerivedV>& corner,
		Eigen::PlainObjectBase<DerivedV>& after);


	template<typename DerivedV>
	void ShapeFunction(const Eigen::PlainObjectBase<DerivedV>& before, const Eigen::PlainObjectBase<DerivedV>& corner, Eigen::PlainObjectBase<DerivedV>& after)
	{
		typedef Eigen::Matrix<typename DerivedV::Scalar, DerivedV::RowsAtCompileTime, DerivedV::ColsAtCompileTime> Mat;

		Mat Vmax, Vmin;
		Vmax = before.colwise().maxCoeff();
		Vmin = before.colwise().minCoeff();

		Mat coord(before.rows(), 3);
		coord = ((before - Vmin.replicate(before.rows(), 1)).array()
			/ (Vmax - Vmin).replicate(before.rows(), 1).array() * 2 - 1).matrix();
		Mat cn(8, 3);
		cn << -1., -1., -1.,
			1., -1., -1.,
			1., 1., -1.,
			-1., 1., -1.,
			-1., -1., 1.,
			1., -1., 1.,
			1., 1., 1.,
			-1., 1., 1.;

		//The shape function for the eight-node	hexahedral brick element
		Mat N(before.rows(), 8);
		for (int i = 0; i < before.rows(); i++) {

			N.row(i) = (coord.row(i).replicate(8, 1).array()*cn.array() + 1).transpose().colwise().prod() / 8;
		}

		after = N * corner;

	}

}
