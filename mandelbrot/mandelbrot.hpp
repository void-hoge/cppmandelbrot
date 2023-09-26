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

//#define ENABLE_AVX

#if defined(ENABLE_AVX) and defined(__AVX2__)
#include <immintrin.h>
#endif

//#define ENABLE_GMP

#if defined(ENABLE_GMP)
#include <gmp.h>
#include <gmpxx.h>
#endif

#define FILL_COUNTMAP

constexpr std::int32_t init = -1;
constexpr std::int32_t queued = -2;

// naive algorithm
template<typename T>
std::vector<std::vector<std::int32_t>> calc_mandelbrot_countmap(
	const std::uint16_t width, const std::uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const std::int32_t iter_max);

// iteration boundary trace algorithm
template<typename T>
std::vector<std::vector<std::int32_t>> calc_mandelbrot_boundary(
	const std::uint16_t width, const std::uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const std::int32_t iter_max,
	const std::pair<std::uint16_t, std::uint16_t> split = {1, 1});
