#include "line_fit.h"
#include <exception>

std::optional<DataBoundingBox> getPointsBox(std::vector<double>* pdata)
{
	auto &data = *pdata;
	if (data.size() < 2) {
		return std::nullopt;
	}

	DataBoundingBox box = { data.at(0), data.at(1), data.at(0), data.at(1) };
	for (size_t i = 0; i < data.size(); i += 2) {
		if (data[i] > box.umax) {
			box.umax = data[i];
		}
		if (data[i] < box.umin) {
			box.umin = data[i];
		}
	}

	for (size_t i = 1; i < data.size(); i += 2) {
		if (data[i] > box.vmax) {
			box.vmax = data[i];
		}
		if (data[i] < box.vmin) {
			box.vmin = data[i];
		}
	}
	return box;
}

LineDistanceCalculator::LineDistanceCalculator()
{
	param_ = { 0.0 };
	param_bounder_ = { {-1.0}, {1.0} };
}

void LineDistanceCalculator::init(std::vector<double>* data)
{
	data_ = data;
}

double LineDistanceCalculator::objective_func(const std::vector<double>& t, std::vector<double>& grad, void * data)
{
	auto param_line_point = [](double a, double t) {
		return std::vector<double>{a*t, 0};
	};

	auto locp_to_glop = [](std::vector<double> local, std::vector<double> loc_matrix) -> std::vector<double> {
		Eigen::Vector3d lp = { local[0], local[1], 1.0 };
		Eigen::Matrix3d loc_m;
		loc_m << loc_matrix[0], loc_matrix[1], loc_matrix[2],
			loc_matrix[3], loc_matrix[4], loc_matrix[5],
			loc_matrix[6], loc_matrix[7], loc_matrix[8];

		Eigen::Vector3d gp = loc_m * lp;

		return { gp[0], gp[1] };
	};

	auto self = reinterpret_cast<LineDistanceCalculator*>(data);
	auto vdata = self->data_;

	double u = vdata->at(0);
	double v = vdata->at(1);
	double theta = vdata->at(2);
	double a = vdata->at(3);

	auto loc_p = param_line_point(a, t.at(0));
	std::vector<double> loc_matrix = {
		cos(theta), -sin(theta),     u,
		sin(theta),  cos(theta),     v,
			   0.0,         0.0,   1.0
	};
	auto glo_p = locp_to_glop(loc_p, loc_matrix);

	std::vector<double> des_p = { vdata->at(4), vdata->at(5) };
	return std::pow(des_p.at(0) - glo_p.at(0), 2) + std::pow(des_p.at(1) - glo_p.at(1), 2);
}

EllipticDistanceCalculator::EllipticDistanceCalculator()
{
	param_ = { 0.0 };
	param_bounder_ = { {0.0}, {2*M_PI} };
}

void EllipticDistanceCalculator::init(std::vector<double>* data)
{
	data_ = data;
}

double EllipticDistanceCalculator::objective_func(const std::vector<double>& t, std::vector<double>& grad, void * data)
{
	auto param_elliptic_point = [](double a, double b, double t) {
		return std::vector<double>{a*cos(t), b*sin(t)};
	};

	auto locp_to_glop = [](std::vector<double> local, std::vector<double> loc_matrix) -> std::vector<double> {
		Eigen::Vector3d lp = { local[0], local[1], 1.0 };
		Eigen::Matrix3d loc_m;
		loc_m << loc_matrix[0], loc_matrix[1], loc_matrix[2],
			loc_matrix[3], loc_matrix[4], loc_matrix[5],
			loc_matrix[6], loc_matrix[7], loc_matrix[8];

		Eigen::Vector3d gp = loc_m * lp;

		return { gp[0], gp[1] };
	};

	auto self = reinterpret_cast<EllipticDistanceCalculator*>(data);
	auto vdata = self->data_;

	double u = vdata->at(0);
	double v = vdata->at(1);
	double theta = vdata->at(2);
	double a = vdata->at(3);
	double b = vdata->at(4);

	auto loc_p = param_elliptic_point(a, b, t.at(0));
	std::vector<double> loc_matrix = {
		cos(theta), -sin(theta),     u,
		sin(theta),  cos(theta),     v,
			   0.0,         0.0,   1.0
	};
	auto glo_p = locp_to_glop(loc_p, loc_matrix);

	std::vector<double> des_p = { vdata->at(4), vdata->at(5) };
	return std::pow(des_p.at(0) - glo_p.at(0), 2) + std::pow(des_p.at(1) - glo_p.at(1), 2);
}

HalfArcDistanceCalculator::HalfArcDistanceCalculator()
{
	param_ = { 0.0 };
	param_bounder_ = { {0.0}, {M_PI} };
}

void HalfArcDistanceCalculator::init(std::vector<double>* data)
{
	data_ = data;
}

double HalfArcDistanceCalculator::objective_func(const std::vector<double>& t, std::vector<double>& grad, void * data)
{
	auto param_arc_point = [](double r, double t) {
		return std::vector<double>{r*cos(t), r*sin(t)};
	};

	auto locp_to_glop = [](std::vector<double> local, std::vector<double> loc_matrix) -> std::vector<double> {
		Eigen::Vector3d lp = { local[0], local[1], 1.0 };
		Eigen::Matrix3d loc_m;
		loc_m << loc_matrix[0], loc_matrix[1], loc_matrix[2],
			loc_matrix[3], loc_matrix[4], loc_matrix[5],
			loc_matrix[6], loc_matrix[7], loc_matrix[8];

		Eigen::Vector3d gp = loc_m * lp;

		return { gp[0], gp[1] };
	};

	auto self = reinterpret_cast<LineDistanceCalculator*>(data);
	auto vdata = self->data_;

	double u = vdata->at(0);
	double v = vdata->at(1);
	double theta = vdata->at(2);
	double r = vdata->at(3);

	auto loc_p = param_arc_point(r, t.at(0));
	std::vector<double> loc_matrix = {
		cos(theta), -sin(theta),     u,
		sin(theta),  cos(theta),     v,
			   0.0,         0.0,   1.0
	};
	auto glo_p = locp_to_glop(loc_p, loc_matrix);

	std::vector<double> des_p = { vdata->at(4), vdata->at(5) };
	return std::pow(des_p.at(0) - glo_p.at(0), 2) + std::pow(des_p.at(1) - glo_p.at(1), 2);
}

OneFourthArcDistanceCalculator::OneFourthArcDistanceCalculator()
{
	param_ = { 0.0 };
	param_bounder_ = { {0.0}, {M_PI_2} };
}

void OneFourthArcDistanceCalculator::init(std::vector<double>* data)
{
	data_ = data;
}

double OneFourthArcDistanceCalculator::objective_func(const std::vector<double>& t, std::vector<double>& grad, void * data)
{
	auto param_arc_point = [](double r, double t) {
		return std::vector<double>{r*cos(t), r*sin(t)};
	};

	auto locp_to_glop = [](std::vector<double> local, std::vector<double> loc_matrix) -> std::vector<double> {
		Eigen::Vector3d lp = { local[0], local[1], 1.0 };
		Eigen::Matrix3d loc_m;
		loc_m << loc_matrix[0], loc_matrix[1], loc_matrix[2],
			loc_matrix[3], loc_matrix[4], loc_matrix[5],
			loc_matrix[6], loc_matrix[7], loc_matrix[8];

		Eigen::Vector3d gp = loc_m * lp;

		return { gp[0], gp[1] };
	};

	auto self = reinterpret_cast<LineDistanceCalculator*>(data);
	auto vdata = self->data_;

	double u = vdata->at(0);
	double v = vdata->at(1);
	double theta = vdata->at(2);
	double r = vdata->at(3);

	auto loc_p = param_arc_point(r, t.at(0));
	std::vector<double> loc_matrix = {
		cos(theta), -sin(theta),     u,
		sin(theta),  cos(theta),     v,
			   0.0,         0.0,   1.0
	};
	auto glo_p = locp_to_glop(loc_p, loc_matrix);

	std::vector<double> des_p = { vdata->at(4), vdata->at(5) };
	return std::pow(des_p.at(0) - glo_p.at(0), 2) + std::pow(des_p.at(1) - glo_p.at(1), 2);
}

LineFitErrorCalculator::LineFitErrorCalculator()
{
	param_ = { 0.0, 0.0, 0.0, 1.0 };
	param_bounder_ = { param_, param_ };
}

void LineFitErrorCalculator::init(std::vector<double>* data)
{
	data_ = data;
	auto box = ::getPointsBox(data_);
	double bounder_half_length = sqrt((box->umax - box->umin)*(box->umax - box->umin) + (box->vmax - box->vmin)*(box->vmax - box->vmin));
	param_bounder_ = { {box->umin, box->vmin, 0.0, 0.0}, {box->umax, box->vmax, M_PI, bounder_half_length} };

}

double LineFitErrorCalculator::objective_func(const std::vector<double>& param, std::vector<double>& grad, void * p_data)
{
	auto self = reinterpret_cast<LineFitErrorCalculator*>(p_data);
	auto vdata = self->data_;

	double fit_error_sum = 0.0;
	int fit_size = vdata->size() / 2;

	for (size_t i = 1; i < vdata->size(); i += 2) {
		Point p{ (*vdata)[i-1], (*vdata)[i] };
		LineDistanceCalculator ldc;

		std::vector<double> param_p = param;
		param_p.push_back(p.x);
		param_p.push_back(p.y);

		ldc.init(&param_p);

		if (auto res_opt = ldc.getFitResult()) {
			fit_error_sum += *res_opt;
		}
		else {
			throw std::range_error("out of computable range!");
		}
	}

	return fit_error_sum / fit_size;
}

FuncParamConverter ParamConvertFuncFactory(FuncParamConverter param_to_shape_mat) {
	auto matxmat = [](std::vector<double> mat1, std::vector<double> mat2) -> std::vector<double> {
		Eigen::Matrix3d eimat1, eimat2;
		eimat1 << mat1[0], mat1[1], mat1[2],
			mat1[3], mat1[4], mat1[5],
			mat1[6], mat1[7], mat1[8];

		eimat2 << mat2[0], mat2[1], mat2[2],
			mat2[3], mat2[4], mat2[5],
			mat2[6], mat2[7], mat2[8];

		Eigen::Matrix3d ret_eimat = eimat1 * eimat2;

		return { ret_eimat(0,0), ret_eimat(0,1), ret_eimat(0,2),
			ret_eimat(1,0), ret_eimat(1,1), ret_eimat(1,2),
			ret_eimat(2,0), ret_eimat(2,1), ret_eimat(2,2) };
	};

	auto getLocalParam = [](std::vector<double> mat, double shape_size) {
		double dx = mat.at(2);
		double dy = mat.at(5);
		double theta = 0.0;
		Eigen::Vector2d dir = Eigen::Vector2d{ mat.at(0), mat.at(3) }.normalized();
		if (dir.y() > 1e-4 && abs(dir.x() / dir.y()) < 1e-4) {
			theta = M_PI_2;
		}
		else if (dir.y() < -1e-4 && abs(dir.x() / dir.y()) < 1e-4) {
			theta = -M_PI_2;
		}
		else if (dir.x() > 1e-4 && abs(dir.y() / dir.x()) < 1e-4) {
			theta = 0.0;
		}
		else if (dir.x() < -1e-4 && abs(dir.y() / dir.x()) < 1e-4) {
			theta = M_PI;
		}
		else {
			theta = atan(dir.y() / dir.x());
		}

		return std::vector<double>{dx, dy, theta, shape_size};
	};

	return [=](std::vector<double> param) -> std::vector<double> {

		std::vector<double> shapes_param2 = { param.at(3),param.at(3), param.at(4), param.at(4) };

		double u = param.at(0);
		double v = param.at(1); 
		double theta = param.at(2);
		std::vector<double> param_mat = {
			cos(theta), -sin(theta),     u,
			sin(theta),  cos(theta),     v,
				   0.0,         0.0,   1.0
		};

		auto shape_mat = param_to_shape_mat(param);
		return getLocalParam(matxmat(param_mat, shape_mat), shape_mat.at(9));
	};
}

CircleFitErrorCalculator::CircleFitErrorCalculator()
{
	param_ = { 0.0, 0.0, 1.0 }; //u, v, r
	param_bounder_ = { param_, param_ };
}

void CircleFitErrorCalculator::init(std::vector<double>* data)
{
	data_ = data;
	auto box = ::getPointsBox(data_);
	double bounder_half_length = sqrt((box->umax - box->umin)*(box->umax - box->umin) + (box->vmax - box->vmin)*(box->vmax - box->vmin));
	param_bounder_ = { {box->umin, box->vmin, 0.0}, {box->umax, box->vmax, bounder_half_length} };
}

double CircleFitErrorCalculator::objective_func(const std::vector<double>& param, std::vector<double>& grad, void * p_data)
{
	auto self = reinterpret_cast<CircleFitErrorCalculator*>(p_data);
	auto vdata = self->data_;

	double fit_error_sum = 0.0;
	int fit_size = vdata->size() / 2;

	for (size_t i = 1; i < vdata->size(); i += 2) {
		Point p{ (*vdata)[i - 1], (*vdata)[i] };
		double distance_squ = std::pow(sqrt(std::pow(p.x - param.at(0), 2) + std::pow(p.y - param.at(1), 2)) - param.at(2), 2);
		fit_error_sum += distance_squ;
	}

	return fit_error_sum / fit_size;
}

EllipticFitErrorCalculator::EllipticFitErrorCalculator()
{
	param_ = { 0.0, 0.0, 0.0, 1.0, 1.0}; //u, v, theta, a, b
	param_bounder_ = { param_, param_ };
}

void EllipticFitErrorCalculator::init(std::vector<double>* data)
{
	data_ = data;
	auto box = ::getPointsBox(data_);
	double bounder_half_length = sqrt((box->umax - box->umin)*(box->umax - box->umin) + (box->vmax - box->vmin)*(box->vmax - box->vmin));
	param_bounder_ = { {box->umin, box->vmin, 0.0, 0.0, 0.0}, {box->umax, box->vmax, M_PI, bounder_half_length, bounder_half_length} };
}

double EllipticFitErrorCalculator::objective_func(const std::vector<double>& param, std::vector<double>& grad, void * p_data)
{
	auto self = reinterpret_cast<EllipticFitErrorCalculator*>(p_data);
	auto vdata = self->data_;

	double fit_error_sum = 0.0;
	int fit_size = vdata->size() / 2;

	for (size_t i = 1; i < vdata->size(); i += 2) {
		Point p{ (*vdata)[i - 1], (*vdata)[i] };
		EllipticDistanceCalculator ldc;

		std::vector<double> param_p = param;
		param_p.push_back(p.x);
		param_p.push_back(p.y);

		ldc.init(&param_p);

		if (auto res_opt = ldc.getFitResult()) {
			fit_error_sum += *res_opt;
		}
		else {
			throw std::range_error("out of computable range!");
		}
	}

	return fit_error_sum / fit_size;
}

ShapeCompositeFitErrorCalculator::ShapeCompositeFitErrorCalculator()
{
}

void ShapeCompositeFitErrorCalculator::acceptCalculators(std::vector<InitParam> init_composite)
{
	calculators_ = std::move(init_composite);
}

double ShapeCompositeFitErrorCalculator::objective_func(const std::vector<double>& param, std::vector<double>& grad, void * p_data)
{
	auto self = reinterpret_cast<ShapeCompositeFitErrorCalculator*>(p_data);
	auto vdata = self->data_;
	//param : dx, dy, theta, a, b

	double fit_error_sum = 0.0;
	int fit_size = vdata->size() / 2;

	for (size_t i = 1; i < vdata->size(); i += 2) {
		std::vector<double> distances(self->calculators_.size(), 0.0f);
		Point p{ (*vdata)[i - 1], (*vdata)[i] };
		for (size_t j = 0; j < self->calculators_.size(); j++) {
			FitBase::Ptr lf = fit_factory_map.at(self->calculators_[j].calculator_name_)->create();
			std::vector<double> loc_param = self->calculators_[j].func_(param);
			loc_param.push_back(p.x);
			loc_param.push_back(p.y);
			lf->init(&loc_param);

			if (auto res_opt = lf->getFitResult()) {
				distances[j] = *res_opt;
			}
			else {
				throw std::range_error("out of computable range!");
			}
		}
		double distance_loc_p = *std::min_element(distances.begin(), distances.end());
		fit_error_sum += distance_loc_p;
	}

	return fit_error_sum / fit_size;
}

RoundWaistShapeFitErrorCalculator::RoundWaistShapeFitErrorCalculator()
{
	param_ = { 0.0, 0.0, 0.0, 1.0, 1.0 }; // u, v, theata, a, b
	auto getMat = [](double theta, double du, double dv, double length) {
		//{cos(theta), -sin(theta), du, sin(theta), cos(theta), dv, 0.0, 0.0, 1.0, length}
		return std::vector<double>{cos(theta), -sin(theta), du, sin(theta), cos(theta), dv, 0.0, 0.0, 1.0, length};
	};
	acceptCalculators({
		//{cos(theata), -sin(theata), du, sin(theta), cos(theta), v, 0.0, 0.0, 1.0, length}
		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(0.0, 0.0, -param.at(4), param.at(3)); })},
		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(0.0, 0.0, param.at(4), param.at(3)); })},
		{std::string("half_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(M_PI_2, -param.at(3), 0.0, param.at(4));}) },
		{std::string("half_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(-M_PI_2, param.at(3), 0.0, param.at(4));}) },
	});
}

void RoundWaistShapeFitErrorCalculator::init(std::vector<double>* data)
{
	data_ = data;
	auto box = ::getPointsBox(data_);
	double bounder_half_length = sqrt((box->umax - box->umin)*(box->umax - box->umin) + (box->vmax - box->vmin)*(box->vmax - box->vmin));
	param_bounder_ = { {box->umin, box->vmin, 0.0, 0.0, 0.0}, {box->umax, box->vmax, M_PI, bounder_half_length, bounder_half_length} };
}

RangleRectangleShapeFitErrorCalculator::RangleRectangleShapeFitErrorCalculator()
{
	param_ = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 }; // u, v, theata, a, b, r
	auto getMat = [](double theta, double du, double dv, double length) {
		//{cos(theta), -sin(theta), du, sin(theta), cos(theta), dv, 0.0, 0.0, 1.0, length}
		return std::vector<double>{cos(theta), -sin(theta), du, sin(theta), cos(theta), dv, 0.0, 0.0, 1.0, length};
	};
	acceptCalculators({

		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(0.0, 0.0, -param.at(4)-param.at(5), param.at(3)); })},

		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(0.0, 0.0, param.at(4)+param.at(5), param.at(3)); })},

		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(M_PI_2, -param.at(3)-param.at(5), 0.0, param.at(4)); })},

		{std::string("line"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(-M_PI_2, param.at(3)+param.at(5), 0.0, param.at(4)); })},

		{std::string("one_fourth_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(0.0, param.at(3), param.at(4), param.at(5)); })},

		{std::string("one_fourth_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(M_PI_2, -param.at(3), param.at(4), param.at(5)); }) },

		{std::string("one_fourth_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(M_PI, -param.at(3), -param.at(4), param.at(5)); }) },

		{std::string("one_fourth_arc"), ParamConvertFuncFactory([=](std::vector<double> param) -> std::vector<double> {
			return getMat(-M_PI_2, param.at(3), -param.at(4), param.at(5)); }) },
		});
}

void RangleRectangleShapeFitErrorCalculator::init(std::vector<double>* data)
{
	data_ = data;
	auto box = ::getPointsBox(data_);
	double bounder_half_length = sqrt((box->umax - box->umin)*(box->umax - box->umin) + (box->vmax - box->vmin)*(box->vmax - box->vmin));
	param_bounder_ = { {box->umin, box->vmin, 0.0, 0.0, 0.0, 0.0}, {box->umax, box->vmax, M_PI, bounder_half_length, bounder_half_length, bounder_half_length} };
}


