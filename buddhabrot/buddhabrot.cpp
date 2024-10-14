#include "buddhabrot.hpp"

uint16_t re2pos(const double real, const double real_min, const double real_unit) {
	return (real - real_min) / real_unit;
}

uint16_t im2pos(const double imag, const double imag_min, const double imag_unit) {
	return (imag - imag_min) / imag_unit;
}

std::pair<uint16_t, uint16_t> complex2pos(
	const double real, const double imag,
	const double real_min, const double real_unit,
	const double imag_min, const double imag_unit) {
	return {re2pos(real, real_min, real_unit), im2pos(imag, imag_min, imag_unit)};
}

int32_t buddhabrot(
	std::vector<double>& xtrace, std::vector<double>& ytrace,
	const double real, const double imag, int32_t iter_max) {

	double x = 0;
	double y = 0;
	double sqx = 0;
	double sqy = 0;
	for (int32_t i = 0; i < iter_max; i++) {
		double tx = sqx - sqy;
		y = x * y * 2.0 + imag;
		x = tx + real;
		xtrace[i] = x;
		ytrace[i] = y;
		sqx = x * x;
		sqy = y * y;
		if (sqx + sqy >= 4.0) {
			return i;
		}
	}
	return iter_max;
}

std::vector<std::vector<int32_t>> buddhabrot_sampling(
	const std::vector<std::vector<int32_t>>& mandelbrot_countmap,
	const uint16_t width, const uint16_t height,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max,
	const int32_t iter_max, const uint32_t samples) {

	const double real_range = real_max - real_min;
	const double real_unit  = real_range / height;
	const double imag_range = imag_max - imag_min;
	const double imag_unit  = imag_range / width;

	std::vector<std::vector<int32_t>> countmap(
		height, std::vector<int32_t>(
			width, 0));

	std::random_device seed;
	std::mt19937 mt(seed());
	std::uniform_real_distribution<double> real_dist(real_min, real_max);
	std::uniform_real_distribution<double> imag_dist(imag_min, imag_max);
	std::vector<double> xtrace(iter_max);
	std::vector<double> ytrace(iter_max);
	for (uint32_t s = 0; s < samples; s++) {
		double re = real_dist(mt);
		double im = imag_dist(mt);
		auto [r, c] = complex2pos(re, im, real_min, real_unit, imag_min, imag_unit);
		if (mandelbrot_countmap[r][c] == iter_max and mandelbrot_countmap[r][c + 1] == iter_max and
			mandelbrot_countmap[r + 1][c] == iter_max and mandelbrot_countmap[r + 1][c + 1] == iter_max) {
			continue;
		}
		int32_t cnt = buddhabrot(xtrace, ytrace, re, im, iter_max);
		if (cnt == iter_max) continue;

		for (int32_t i = 0; i < cnt; i++) {
			uint16_t row = re2pos(xtrace[i], real_min, real_unit);
			uint16_t col = im2pos(ytrace[i], imag_min, imag_unit);
			if (row < height and col < width) {
				countmap[row][col]++;
			}
		}
	}
	return countmap;
}

std::vector<std::vector<int32_t>> calc_buddhabrot_countmap(
	const uint16_t width, const uint16_t height,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max,
	const int32_t iter_max, const uint32_t samples,
	const uint16_t split, region_manager& region) {

	std::cerr << samples << std::endl;

	const double real_range = real_max - real_min;
	const double real_unit  = real_range / height;
	const double imag_range = imag_max - imag_min;
	const double imag_unit  = imag_range / width;

	auto& mandelbrot_countmap = region.countmap;

	std::vector<std::vector<std::vector<int32_t>>> countmaps(split);

#pragma omp parallel for
	for (uint32_t i = 0; i < split; i++) {
		countmaps[i] = buddhabrot_sampling(
			mandelbrot_countmap, width, height,
			real_min, real_max, imag_min, imag_max, iter_max, samples / split);
	}

	std::cerr << "sampling finished" << std::endl;

	auto countmap = std::vector<std::vector<int32_t>>(
		height, std::vector<int32_t>(
			width, 0));
	for (uint32_t i = 0; i < split; i++) {
		for (uint32_t r = 0; r < height; r++) {
			for (uint32_t c = 0; c < width; c++) {
#if defined(DOUBLECOUNT)
				countmap[r][c] += countmaps[i][r][c] + countmaps[i][r][width - c - 1];
#else
				countmap[r][c] += countmaps[i][r][c];
#endif
			}
		}
	}
	return countmap;
}
