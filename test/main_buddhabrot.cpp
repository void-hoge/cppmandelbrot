#include "../buddhabrot/buddhabrot.hpp"

#include <iostream>
#include <chrono>

template<typename T>
std::ostream& operator<<(std::ostream& ost, std::vector<T> vec) {
	for (std::size_t i = 0; i < vec.size(); i++) {
		ost << vec[i] << (i == vec.size()-1 ? "" : " ");
	}
	ost << std::endl;
	return ost;
}

int main() {
	const uint16_t width   = 1920;
	const uint16_t height  = 1080;
	const double real_min       = -1.5;
	const double real_max       = 0.5;
	const double real_range     = real_max - real_min;
	const double imag_range     = real_range * ((double)width / height);
	const double imag_min       = -imag_range / 2;
	const double imag_max       = imag_range / 2;
	const int32_t iter_max = 1000000;
	region_manager region(width, height, 16, 1);
	calc_mandelbrot_boundary(
		width, height, real_min, real_max, imag_min, imag_max, iter_max, region);
	auto countmap = calc_buddhabrot_countmap(
		width, height, real_min, real_max, imag_min, imag_max,
		iter_max, width * height * 16, 16, region);
	std::cout << width << " " << height << std::endl;
	std::cout << countmap << std::endl;
}
