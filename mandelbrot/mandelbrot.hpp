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
#include <cstdio>

#define ENABLE_AVX

#if defined(ENABLE_AVX) and defined(__AVX2__)
#include <immintrin.h>
#endif

#define FILL_COUNTMAP

constexpr std::int32_t init = -1;
constexpr std::int32_t queued = -2;

// naive algorithm
std::vector<std::vector<std::int32_t>> calc_mandelbrot_countmap(
	const std::uint16_t width, const std::uint16_t height,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max,
	const std::int32_t iter_max);

// iteration boundary trace algorithm
std::vector<std::vector<std::int32_t>> calc_mandelbrot_boundary(
	const std::uint16_t width, const std::uint16_t height,
	const double real_min, const double real_max,
	const double imag_min, const double imag_max,
	const std::int32_t iter_max,
	const std::pair<std::uint16_t, std::uint16_t> split = {1, 1});
