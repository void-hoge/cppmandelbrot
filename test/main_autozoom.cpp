#include "../mandelbrot/mandelbrot.hpp"
#include "../autozoom/autozoom.hpp"

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
	double zoom = 3.0;
	double real_center = 0.0;
	double imag_center = 0.0;
	for (int i = 0; i < 20; i++) {
		const double real_min = real_center - zoom / 2;
		const double real_max = real_center + zoom / 2;
		const double imag_min = imag_center - zoom / 2 * ((double)width / height);
		const double imag_max = imag_center + zoom / 2 * ((double)width / height);
		region_manager region(width, height, 1, 16);
		calc_mandelbrot_boundary(
			width, height, real_min, real_max, imag_min, imag_max, iter_max, region);

		auto [rc, ic] = get_most_intricate_place(
			region.countmap, width, height, 5, 5,
			real_min, real_max, imag_min, imag_max);
		real_center = rc;
		imag_center = ic;
		zoom *= 0.8;
		std::cout << zoom << " " << real_center << " " << imag_center << std::endl;
	}
}
