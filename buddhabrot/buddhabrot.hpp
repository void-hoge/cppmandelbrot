#pragma once

#include "../mandelbrot/mandelbrot.hpp"

#include <random>

#define DOUBLECOUNT

std::vector<std::vector<std::int32_t>> calc_buddhabrot_countmap(
	const std::uint16_t width, const std::uint16_t height,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max,
	const std::int32_t iter_max, const std::uint32_t samples,
	const std::uint16_t split, region_manager& region);
