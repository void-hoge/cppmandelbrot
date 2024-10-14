#pragma once

#include <vector>
#include <cstdint>

std::pair<double, double> get_most_intricate_place(
	const std::vector<std::vector<int32_t>>& countmap,
	const uint16_t subwidth, const uint16_t subheight,
	const uint16_t col_split, const uint16_t row_split,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max);
