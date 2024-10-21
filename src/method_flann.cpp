#pragma once

#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "flann/flann.h"

struct OutEdge {
	size_t index_;
	//int* edge_num = nullptr;

	std::stack<size_t> * line_ids_;

	bool operator==(const OutEdge& that) const {
		return index_ == that.index_;
	}
};

//使用图的邻接表表示自动连接多段线。
void test_flann() {

	std::vector<float> points;
	int dimension = 3;
	int points_size = 6;

	points = {
		0.f, 0.f, 0.f,
		0.f, 2.f, 0.f, //line1
		0.f, 2.f, 0.f,
		2.f, 2.f, 0.f, //line2
		2.f, 2.f, 0.f,
		0.f, 0.f, 0.f  //line3
	};

	//for (int i = 0; i < dimension*points_size; i++) {
	//	points.push_back(i / 10.f);
	//}

	auto temp_matrix = flann::Matrix<float>(points.data(), points_size, dimension);
	flann::Index<flann::L2<float>> index(temp_matrix, flann::KDTreeIndexParams());

	index.buildIndex();

	std::vector<std::vector<size_t>> indices;
	std::vector<std::vector<float>> dist;
	index.radiusSearch(temp_matrix, indices, dist, 0.01f, flann::SearchParams());

	std::unordered_map<size_t, size_t> map_to_unique_index;
	std::unordered_map<size_t, std::vector<OutEdge>> map_to_outedge;
	std::vector<std::unique_ptr<int>> out_edges;
	std::vector<std::unique_ptr<std::stack<size_t>>> out_edges_;

	for (size_t i = 0; i < indices.size(); i++) {
		if (indices[i].size() < 1) continue;
		size_t unique_index = indices[i][0];
		for (size_t j = 0; j < indices[i].size(); j++) {
			map_to_unique_index.try_emplace(indices[i][j], unique_index);
			//map_to_unique_index[indices[i][j]] = unique_index;
		}
	}

	for (size_t i = 0; i < points_size; i += 2) {
		size_t uindex = map_to_unique_index[i];
		size_t u2index = map_to_unique_index[i + 1];
		if (map_to_outedge.find(uindex) != map_to_outedge.end()) {
			auto & list = map_to_outedge[uindex];
			auto it = std::find(list.begin(), list.end(), OutEdge{ u2index });
			if (it != list.end()) {
				//(*it->edge_num)++;
				it->line_ids_->push(i / 2);
			}
			else {
				out_edges.push_back(std::make_unique<int>(1));
				out_edges_.emplace_back(new std::stack<size_t>({ i / 2 }));
				map_to_outedge[uindex].push_back({ u2index,/* out_edges.back().get(),*/ out_edges_.back().get() });
				map_to_outedge[u2index].push_back({ uindex,/* out_edges.back().get(),*/ out_edges_.back().get() });
			}
		}
		else {
			out_edges.push_back(std::make_unique<int>(1));
			out_edges_.emplace_back(new std::stack<size_t>({ i / 2 }));
			map_to_outedge[uindex].push_back({ u2index,/* out_edges.back().get(),*/ out_edges_.back().get() });
			map_to_outedge[u2index].push_back({ uindex,/* out_edges.back().get(),*/ out_edges_.back().get() });
		}
	}

	std::vector<std::vector<size_t>> lines;

	for (auto it = map_to_outedge.begin(); it != map_to_outedge.end(); it++) {
		size_t key = it->first;
		auto & vec = it->second;
		bool can_add = false;
		for (size_t i = 0; i < vec.size(); i++) {
			if (vec[i].line_ids_->size() > 0) {
				lines.push_back({ key });

				for (bool is_found = true; is_found; ) {
					is_found = false;
					auto & list = map_to_outedge[lines.back().back()];
					for (auto & edge : list) {
						if (edge.line_ids_->size() > 0) {
							lines.back().push_back(edge.line_ids_->top());
							edge.line_ids_->pop();
							lines.back().push_back(edge.index_);
							is_found = true;
							break;
						}
					}
				}
			}
		}
	}

	return;

}