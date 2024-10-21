#pragma once
#include <vector>
#include <unordered_map>
#include <list>
#include <optional>
#include <Eigen/Dense>

#define N 5

struct PolyLine {

	enum class line_type {
		arc,
		line
	};

	std::list<size_t> vertex_;
	std::list<line_type> vec_type_;

	PolyLine& reverse() {
		vertex_.reverse();
		vec_type_.reverse();
		return *this;
	}

	void addAnothorPline(PolyLine && pl) {
		if (pl.vertex_.back() == vertex_.back()) {
			pl.reverse();
			vertex_.splice(vertex_.end(), pl.vertex_);
			vec_type_.splice(vec_type_.end(), pl.vec_type_);
		}
		else if (pl.vertex_.front() == vertex_.front()) {
			pl.reverse();
			vertex_.splice(vertex_.begin(), pl.vertex_);
			vec_type_.splice(vec_type_.begin(), pl.vec_type_);
		}
		else if (pl.vertex_.back() == vertex_.front()) {
			vertex_.splice(vertex_.begin(), pl.vertex_);
			vec_type_.splice(vec_type_.begin(), pl.vec_type_);
		}
		else {
			vertex_.splice(vertex_.end(), pl.vertex_);
			vec_type_.splice(vec_type_.end(), pl.vec_type_);
		}
	}

};

std::vector<PolyLine> all_line;

static int initPoints = ([]() {
	all_line = { {{1,7}, {PolyLine::line_type::line}},
					{{4,9}, {PolyLine::line_type::line}},
					{{1,4}, {PolyLine::line_type::line}},
					{{9,7}, {PolyLine::line_type::line}},
	};
}(), 0);

struct PointPair {
	size_t pline_id;

	size_t anotherPt(size_t current_pt) {
		size_t fr = all_line.at(pline_id).vertex_.front();
		size_t bc = all_line.at(pline_id).vertex_.back();
		if (current_pt == fr) {
			return bc;
		}
		else if (current_pt == bc) {
			return fr;
		}
		else {
			return 0;
		}
	}
};

void test() {

	std::unordered_map<size_t, PointPair> point_map;

	for (size_t i = 0; i < all_line.size(); i++) {
		size_t front = all_line[i].vertex_.front();
		size_t back = all_line[i].vertex_.back();
		if (point_map.find(front) != point_map.end() && point_map.find(back) != point_map.end()) {
			size_t temp = point_map[front].anotherPt(front);
			point_map[point_map[back].anotherPt(back)] = point_map[temp];
			all_line[point_map[temp].pline_id].addAnothorPline(std::move(all_line[i]));
			all_line[point_map[temp].pline_id].addAnothorPline(std::move(all_line[point_map[back].pline_id]));
			point_map.erase(front);
			point_map.erase(back);
		}
		else if (point_map.find(front) != point_map.end() && point_map.find(back) == point_map.end()) {
			size_t temp = point_map[front].anotherPt(front);
			all_line[point_map[temp].pline_id].addAnothorPline(std::move(all_line[i]));
			point_map.erase(front);
			point_map[back] = point_map[temp];
		}
		else if (point_map.find(front) == point_map.end() && point_map.find(back) != point_map.end()) {
			size_t temp = point_map[front].anotherPt(back);
			all_line[point_map[temp].pline_id].addAnothorPline(std::move(all_line[i]));
			point_map.erase(back);
			point_map[front] = point_map[temp];
		}
		else {
			point_map[front] = PointPair{ i };
			point_map[back] = PointPair{ i };
		}
	}

}