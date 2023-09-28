#include "../mandelbrot/mandelbrot.hpp"

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
#if defined(ENABLE_GMP)
	mpf_set_default_prec(256);
	const std::uint16_t width   = 1920;
	const std::uint16_t height  = 1080;
	const mpf_class real_min       = mpf_class(-1.5);
	const mpf_class real_max       = mpf_class(0.5);
	const mpf_class real_range     = real_max - real_min;
	const mpf_class imag_range     = real_range * ((double)width / height);
	const mpf_class imag_min       = -imag_range / 2;
	const mpf_class imag_max       = imag_range / 2;
	const std::int32_t iter_max = 100000;
#else
	// double zoom;
	// std::cin >> zoom;
	// const std::uint16_t height  = 1080 * 2;
	// const std::uint16_t width   = 1920 * 2;
	// const double real_center    = -1.4855241;
	// const double real_min       = real_center - zoom;
	// const double real_max       = real_center + zoom;
	const std::uint16_t width   = 1920;
	const std::uint16_t height  = 1080;
	const double real_min       = -1.5;
	const double real_max       = 0.5;
	const double real_range     = real_max - real_min;
	const double imag_range     = real_range * ((double)width / height);
	const double imag_min       = -imag_range / 2;
	const double imag_max       = imag_range / 2;
	const std::int32_t iter_max = 100000;
#endif
	region_manager region(width, height, 1, 16);
	calc_mandelbrot_boundary(
		width, height, real_min, real_max, imag_min, imag_max, iter_max, region);
	std::cout << width << " " << height << std::endl;
	std::cout << region.countmap << std::endl;
}
