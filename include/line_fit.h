#pragma once
#include <optional>

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

struct Point {
	double x;
	double y;
};

struct DataBoundingBox {
	double umin, vmin, umax, vmax;
};

std::optional<DataBoundingBox> getPointsBox(std::vector<double>* pdata);

using FuncParamConverter = std::function<std::vector<double>(std::vector<double>)>;

FuncParamConverter ParamConvertFuncFactory(FuncParamConverter param_to_shape_mat);

class FitBase {
public:
	using Ptr = std::shared_ptr<FitBase>;

	virtual ~FitBase() = default;

	virtual void init(std::vector<double>* data) = 0;

	virtual std::optional<double> getFitResult() = 0;

public:
	std::vector<double> param_;
	std::vector<std::vector<double>> param_bounder_;
	std::vector<double>* data_ = nullptr;
};

template<class Derived, class Base = FitBase>
class FitBaseMin : public Base {
public:
	std::optional<double> getFitResult() override {
		if (data_ == nullptr) {
			std::cerr << "please set data to fit.\n";
			return std::nullopt;
		}
		//LN_BOBYQA; LN_COBYLA
		auto opt = std::make_unique<nlopt::opt>(nlopt::LN_BOBYQA, param_.size());

		opt->set_min_objective(Derived::objective_func, this);
		opt->set_lower_bounds(param_bounder_.front());
		opt->set_upper_bounds(param_bounder_.back());
		opt->set_xtol_abs(1e-4);
		opt->set_maxeval(1000000);

		double fit_result;
		auto result = opt->optimize(param_, fit_result);

		if (!(result > 0)) {
			return std::nullopt;
		}

		return fit_result;
	}
};


class FitFactoryBase {
public:
	using Ptr = std::shared_ptr<FitFactoryBase>;
	virtual FitBase::Ptr create() = 0;
	~FitFactoryBase() = default;
};

template<class Fit>
class FitFactoryImpl : public FitFactoryBase {
public:
	FitBase::Ptr create() override {
		return std::make_shared<Fit>();
	}
};



class LineDistanceCalculator : public FitBaseMin<LineDistanceCalculator>{
public:
	LineDistanceCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& t, std::vector<double>& grad, void* data);

};

class EllipticDistanceCalculator : public FitBaseMin<LineDistanceCalculator> {
public:
	EllipticDistanceCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& t, std::vector<double>& grad, void* data);

};

class HalfArcDistanceCalculator : public FitBaseMin<HalfArcDistanceCalculator> {
public:
	HalfArcDistanceCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& t, std::vector<double>& grad, void* data);

};

class OneFourthArcDistanceCalculator : public FitBaseMin<HalfArcDistanceCalculator> {
public:
	OneFourthArcDistanceCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& t, std::vector<double>& grad, void* data);

};

class LineFitErrorCalculator : public FitBaseMin<LineFitErrorCalculator> {
public:

	LineFitErrorCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& param, std::vector<double>& grad, void* p_data);
};

class CircleFitErrorCalculator : public FitBaseMin<CircleFitErrorCalculator> {
public:
	CircleFitErrorCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& param, std::vector<double>& grad, void* p_data);
};

class EllipticFitErrorCalculator : public FitBaseMin<CircleFitErrorCalculator> {
public:
	EllipticFitErrorCalculator();

	void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& param, std::vector<double>& grad, void* p_data);
};

class ShapeCompositeFitErrorCalculator : public FitBaseMin<ShapeCompositeFitErrorCalculator> {
public:
	inline static const std::unordered_map<std::string, FitFactoryBase::Ptr> fit_factory_map = {
		{"line", std::make_shared<FitFactoryImpl<LineDistanceCalculator>>()},
		{"half_arc", std::make_shared<FitFactoryImpl<HalfArcDistanceCalculator>>()},
		{"one_fourth_arc", std::make_shared<FitFactoryImpl<OneFourthArcDistanceCalculator>>()},
	};

	struct InitParam {
		std::string calculator_name_;
		FuncParamConverter func_;
	};

	ShapeCompositeFitErrorCalculator();

	void acceptCalculators(std::vector<InitParam> init_composite);

	//void init(std::vector<double>* data) override;

	static double objective_func(const std::vector<double>& param, std::vector<double>& grad, void* p_data);

	std::vector<InitParam> calculators_;
};

class RoundWaistShapeFitErrorCalculator : public ShapeCompositeFitErrorCalculator {
public:
	RoundWaistShapeFitErrorCalculator();

	void init(std::vector<double>* data) override;
};

class RangleRectangleShapeFitErrorCalculator : public ShapeCompositeFitErrorCalculator {
public:
	RangleRectangleShapeFitErrorCalculator();

	void init(std::vector<double>* data) override;
};