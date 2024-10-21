#pragma once

#include <nlopt.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>

#ifndef _MATH_DEFINES_DEFINED
#define _MATH_DEFINES_DEFINED
// Definitions of useful mathematical constants
//
// Define _USE_MATH_DEFINES before including <math.h> to expose these macro
// definitions for common math constants.  These are placed under an #ifdef
// since these commonly-defined names are not part of the C or C++ standards
#define M_E        2.71828182845904523536   // e
#define M_LOG2E    1.44269504088896340736   // log2(e)
#define M_LOG10E   0.434294481903251827651  // log10(e)
#define M_LN2      0.693147180559945309417  // ln(2)
#define M_LN10     2.30258509299404568402   // ln(10)
#define M_PI       3.14159265358979323846   // pi
#define M_PI_2     1.57079632679489661923   // pi/2
#define M_PI_4     0.785398163397448309616  // pi/4
#define M_1_PI     0.318309886183790671538  // 1/pi
#define M_2_PI     0.636619772367581343076  // 2/pi
#define M_2_SQRTPI 1.12837916709551257390   // 2/sqrt(pi)
#define M_SQRT2    1.41421356237309504880   // sqrt(2)
#define M_SQRT1_2  0.707106781186547524401  // 1/sqrt(2)
#endif

//定义局部坐标系下线段的参数化方式: 坐标轴上长为2a, 参数为t的线段
std::vector<double> param_line_point(double a, double t) {
	return { a*t, 0. };
}

std::vector<double> p_to_global(std::vector<double> local, std::vector<double> loc_matrix) {
	if (local.size() < 2 || loc_matrix.size() < 9) {
		return {};
	}

	Eigen::Vector3d lp = { local[0], local[1], 1.0 };
	Eigen::Matrix3d loc_m;
	loc_m << loc_matrix[0], loc_matrix[1], loc_matrix[2],
			 loc_matrix[3], loc_matrix[4], loc_matrix[5],
			 loc_matrix[6], loc_matrix[7], loc_matrix[8];

	Eigen::Vector3d gp = loc_m * lp;

	return { gp[0], gp[1] };
}
class LineFunc {
public:
	static double objective_function(const std::vector<double>& x, std::vector<double>& grad, void* data) {
		//auto & vdata = *reinterpret_cast<std::vector<double> *>(data);
		auto vdata = reinterpret_cast<LineFunc *>(data);

		double theta = vdata->theta;
		double a = vdata->a;
		double t = x.at(0);

		auto loc_p = param_line_point(a, t);

		std::vector<double> loc_matrix = {
			cos(theta), -sin(theta), vdata->u,
			sin(theta),  cos(theta), vdata->v,
				   0.0,         0.0,         1.0
		};

		auto glo_p = p_to_global(loc_p, loc_matrix);

		std::vector<double> des_p = { vdata->x, vdata->y };
		double min_distance_squared = std::pow(des_p.at(0) - glo_p.at(0), 2) + std::pow(des_p.at(1) - glo_p.at(1), 2);
		return min_distance_squared;
	}

public:
	double u = 0.;
	double v = 0.;
	double theta = 0.;
	double a = 1.;
	double x = 0.;
	double y = 0.;
};



void test_getDistance() {

	std::vector<double> x = { 0.0 }; 

	// 创建一个nlopt优化器///LN_COBYLA//LN_BOBYQA
	nlopt::opt optimizer(nlopt::LN_BOBYQA, x.size());

	LineFunc line{ 1.0, 1.0, M_PI_2 / 3, 5.0, 4.0, 10.0 };
	optimizer.set_min_objective(LineFunc::objective_function, &line);

	optimizer.set_lower_bounds({ -1.0 });
	optimizer.set_upper_bounds({ 1.0 });
	optimizer.set_xtol_abs(1e-4);
	optimizer.set_maxeval(1000000);

	double min_distance_squared;
	auto result = optimizer.optimize(x, min_distance_squared);
	if (result > 0) {
	
		auto loc_p = param_line_point(line.a, x.at(0));
		auto theta = line.theta;
		std::vector<double> loc_matrix = {
			cos(theta), -sin(theta),  line.u,
			sin(theta),  cos(theta),  line.v,
				   0.0,         0.0,         1.0
		};

		auto glo_p = p_to_global(loc_p, loc_matrix);

		std::cout << "distance: " << std::sqrt(std::pow(glo_p.at(0) - line.x, 2) + std::pow(glo_p.at(1) - line.y, 2)) << "\n";
	
	}
}
