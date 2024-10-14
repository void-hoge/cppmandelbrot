#include "autozoom.hpp"

#include <iostream>
#include <algorithm>

std::pair<double, double> get_most_intricate_place(
	const std::vector<std::vector<int32_t>>& countmap,
	const uint16_t subwidth, const uint16_t subheight,
	const uint16_t col_split, const uint16_t row_split,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max) {

	const uint16_t width = countmap.at(0).size();
	const uint16_t height = countmap.size();

	if (subwidth > width) throw std::logic_error("Target must be smaller or equals to countmap.");
	if (subheight > height) throw std::logic_error("Target must be smaller or equals to countmap.");
	if (subwidth % col_split != 0) throw std::logic_error("Number of row split must be a divisor of subwidth.");
	if (subheight % row_split != 0) throw std::logic_error("Number of col split must be a divisor of subheight.");

	auto sums = std::vector(
		row_split, std::vector<uint64_t>(
			col_split, 0));

	const uint16_t rowbegin = (width - subwidth) / 2;
	const uint16_t colbegin = (height - subheight) / 2;
	for (uint16_t row = 0; row < subheight; row++) {
		auto rowidx = row / (subheight / row_split);
		for (uint16_t col = 0; col < subwidth; col++) {
			auto colidx = col / (subwidth / col_split);
			sums[rowidx][colidx] += countmap[rowbegin + row][colbegin + col];
		}
	}

	auto averages = std::vector(
		row_split, std::vector<double>(
			col_split, 0));
	auto area = (subheight / row_split) * (subwidth / col_split);
	for (uint32_t rowidx = 0; rowidx < row_split; rowidx++) {
		for (uint32_t colidx = 0; colidx < col_split; colidx++) {
			averages[rowidx][colidx] = (double)sums[rowidx][colidx] / area;
		}
	}

	auto variances = std::vector(
		row_split, std::vector<double>(
			col_split, 0));
	for (uint16_t row = 0; row < subheight; row++) {
		auto rowidx = row / (subheight / row_split);
		for (uint16_t col = 0; col < subwidth; col++) {
			auto colidx = col / (subwidth / col_split);
			auto deviation = countmap[rowbegin + row][colbegin + col] * area - averages[rowidx][colidx];
			variances[rowidx][colidx] += deviation * deviation;
		}
	}

	uint32_t maxrow = 0;
	uint32_t maxcol = 0;
	double max_variance = variances[maxrow][maxcol];
	for (uint32_t rowidx = 0; rowidx < row_split; rowidx++) {
		auto max_in_row = std::max_element(variances[rowidx].begin(), variances[rowidx].end());
		if (*max_in_row > max_variance) {
			maxrow = rowidx;
			maxcol = std::distance(variances[rowidx].begin(), max_in_row);
			max_variance = *max_in_row;
		}
	}

	std::cout << maxcol << " " << maxrow << std::endl;

	const double sub_real_max = real_max * ((double)subheight / height);
	const double sub_real_min = real_min * ((double)subheight / height);
	const double sub_real_range_unit = (sub_real_max - sub_real_min) / row_split;
	const double ret_real_center = sub_real_min + sub_real_range_unit * maxrow + sub_real_range_unit / 2;

	const double sub_imag_max = imag_max * ((double)subwidth / width);
	const double sub_imag_min = imag_min * ((double)subwidth / width);
	const double sub_imag_range_unit = (sub_imag_max - sub_imag_min) / col_split;
	const double ret_imag_center = sub_imag_min + sub_imag_range_unit * maxcol + sub_imag_range_unit / 2;

	return {ret_real_center, ret_imag_center};
}
