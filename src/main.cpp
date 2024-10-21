#include <iostream>
#include "line_fit.h"
void test_flann();

void test_fitCurve();

void test_getDistance();

int main() {

	//test_flann();

	//test_fitCurve();

	//test_getDistance();

	/*std::vector<double> points{
		10.0, 0.0,
		10.0, 5.0,
		8.54, 8.54,
		5.0, 10.0,
		0.0, 10.0,
		-5.0, 10.0,
		-8.54, 8.54,
		-10.0, 5.0,
		-10.0, 0.0,
		-10.0, -5.0,
		-8.54, -8.54,
		-5.0, -10.0,
		0.0, -10.0,
		5.0, -10.0,
		8.54, -8.54,
		10.0, -5.0,
	};*/

	std::vector<double> points{
	15.0, 0.0,
	13.54, 3.54,
	10.0, 5.0,
	0.0, 5.0,
	-10.0, 5.0,
	-13.54, 3.54,
	-15.0, 0.0,
	-13.54, -3.54,
	-10.0, -5.0,
	0.0, -5.0,
	10.0, -5.0,
	13.54, -3.54,
	};

	//LineFitErrorCalculator line_fit;
	//line_fit.init(&points);
	//auto error = line_fit.getFitResult();

	auto fit_ptr = std::make_unique<RoundWaistShapeFitErrorCalculator>();
	fit_ptr->init(&points);
	auto error = fit_ptr->getFitResult();

	//RangleRectangleShapeFitErrorCalculator rec_fit;
	//rec_fit.init(&points);
	//auto error = rec_fit.getFitResult();
	std::cout << *error;

	return 0;
}