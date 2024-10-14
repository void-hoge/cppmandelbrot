#pragma once

#include <vector>
#include <cstdint>
#include <iostream>
#include <queue>
#include <set>
#include <array>
#include <cassert>
#include <bitset>
#include <chrono>
#include <type_traits>

#define ENABLE_AVX

#if defined(ENABLE_AVX) and defined(__AVX2__)
#include <immintrin.h>
#endif

// #define ENABLE_GMP

#if defined(ENABLE_GMP)
#include <gmp.h>
#include <gmpxx.h>
#endif

#define FILL_COUNTMAP

constexpr int32_t INIT = -1;
constexpr int32_t QUEUED = -2;

class region_manager{
public:
	const uint16_t height;
	const uint16_t width;
	const uint16_t col_split;
	const uint16_t row_split;
	std::vector<std::vector<int32_t>> countmap;
	std::vector<std::vector<int32_t *>> subcountmaps;
	std::vector<std::vector<int32_t>> col_perimeters;
	std::vector<std::vector<int32_t>> row_perimeters;
	std::vector<uint16_t> cols;
	std::vector<uint16_t> rows;

	region_manager(
		uint16_t width, uint16_t height,
		uint16_t col_split = 1, uint16_t row_split = 1);

	region_manager(const region_manager& region);

	~region_manager();
};

class nullbuffer : public std::streambuf {
public:
	int overflow(int c) override {
		return c;
	}
};

class nullstream : public std::ostream {
public:
	nullstream() : std::ostream(&nullbuff) {}
private:
	nullbuffer nullbuff;
};

static nullstream nullout;

// naive algorithm
template<typename T>
std::vector<std::vector<int32_t>> calc_mandelbrot_countmap(
	const uint16_t width, const uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const int32_t iter_max, std::ostream& ost = std::cerr);

// iteration boundary trace algorithm
template<typename T>
void calc_mandelbrot_boundary(
	const uint16_t width, const uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const int32_t iter_max, region_manager& region, std::ostream& ost = std::cerr);
