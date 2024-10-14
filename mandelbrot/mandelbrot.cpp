#include "mandelbrot.hpp"

// Escape time engine (scalar version)
template<typename T>
uint32_t mandelbrot(
	const T& a, const T& b,
	const int32_t iter_max) {

	T x   = 0;
	T y   = 0;
	T sqx = 0;
	T sqy = 0;
	for (int32_t i = 0; i < iter_max; i++) {
		T tx = sqx - sqy;
		y   = x * y * 2.0 + b;
		x   = tx + a;
		sqx = x * x;
		sqy = y * y;
		if (sqx + sqy >= 4.0) {
			return i;
		}
	}
	return iter_max;
}

#if defined(ENABLE_AVX) and defined(__AVX2__)
// Escape time engine (AVX version)
__m256i mandelbrot_avx(
	const __m256d& a, const __m256d& b,
	const int32_t iter_max) {

	__m256d x    = _mm256_setzero_pd();
	__m256d y    = _mm256_setzero_pd();
	__m256i cnt  = _mm256_setzero_si256();
	for (int32_t i = 0; i < iter_max; i++) {
		__m256d tx = _mm256_fmsub_pd(x, x, _mm256_mul_pd(y, y));
		y = _mm256_fmadd_pd(_mm256_mul_pd(x, y), _mm256_set1_pd(2.0), b);
		x = _mm256_add_pd(tx, a);
		__m256d cmp       = _mm256_cmp_pd(_mm256_fmadd_pd(x, x, _mm256_mul_pd(y, y)), _mm256_set1_pd(4.0), _CMP_LT_OQ);
		__m256i mask      = _mm256_castpd_si256(cmp);
		__m256i increment = _mm256_and_si256(mask, _mm256_set1_epi64x(1));
		cnt = _mm256_add_epi64(cnt, increment);
		if (_mm256_testz_si256(mask, mask)) {
			break;
		}
	}
	return cnt;
}
#endif

// Naive algorithm: compute escape times for each pixel
template<typename T>
std::vector<std::vector<int32_t>> calc_mandelbrot_countmap(
	const uint16_t width, const uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const int32_t iter_max) {

	std::vector<std::vector<int32_t>> countmap(
		height, std::vector<int32_t>(width, 0));
	const T real_range = real_max - real_min;
	const T imag_range = imag_max - imag_min;
	const T real_unit  = real_range / height;
	const T imag_unit  = imag_range / width;
	auto start = std::chrono::system_clock::now();
#if defined(ENABLE_AVX) and defined(__AVX2__) and not defined(ENABLE_GMP)
	static_assert(std::is_same<T, double>::value);
#pragma omp parallel for
	for (uint32_t row = 0; row < height; row++) {
		for (uint32_t col = 0; col < width; col += 4) {
			alignas(32) uint64_t store[4];
			double real = real_min + real_unit * row;
			double imag = imag_min + imag_unit * col;
			__m256d a   = _mm256_set1_pd(real);
			__m256d b   = _mm256_set_pd(
				imag + imag_unit * 3,
				imag + imag_unit * 2,
				imag + imag_unit * 1,
				imag + imag_unit * 0);
			__m256i cnt = mandelbrot_avx(a, b, iter_max);
			_mm256_storeu_si256((__m256i *)store, cnt);
			for (uint32_t i = 0; i < 4; i++) {
				countmap[row][col + i] = store[i];
			}
		}
	}
#else
#pragma omp parallel for
	for (uint16_t row = 0; row < height; row++) {
		for (uint16_t col = 0; col < width; col++) {
			T real = real_min + real_unit * row;
			T imag = imag_min + imag_unit * col;
			countmap[row][col] = mandelbrot(real, imag, iter_max);
		}
	}
#endif
	auto end = std::chrono::system_clock::now();
	double elapsed = (double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() / 1000;
	std::cerr << elapsed << std::endl;
	return countmap;
}

#if defined(ENABLE_AVX) and defined(__AVX2__)
// Pop 4 pixels from QUE, convert to complex numbers,
// and pack them into two __m256d REAL and IMAG.
uint32_t pop4elms(
	std::deque<std::pair<uint16_t, uint16_t>>& que,
	const int32_t *countmap,
	const uint32_t mapwidth, const uint32_t mapheight,
	const double real_min, const double real_unit,
	const double imag_min, const double imag_unit,
	__m256d& real, __m256d& imag,
	std::array<uint16_t, 4>& rpos, std::array<uint16_t, 4>& cpos) {

	alignas(32) double preal[4];
	alignas(32) double pimag[4];
	uint32_t idx = 0;
	while (que.size() and idx < 4) {
		auto&& [row, col] = que.front();
		que.pop_front();
		if (countmap[row * mapwidth + col] >= 0) continue;
		rpos[idx]  = row;
		cpos[idx]  = col;
		preal[idx] = real_min + real_unit * row;
		pimag[idx] = imag_min + imag_unit * col;
		idx++;
	}
	real = _mm256_loadu_pd(preal);
	imag = _mm256_loadu_pd(pimag);
	return idx;
}
#endif

// Compute Mandelbrot set escape time boundary map
// The escape times of perimeter given as UPPER, LOWER, LEFT and RIGHT.
// The result is stored into the COUNRMAP.
// COUNTMAP is allocated in caller side.
template<typename T>
void calc_submandelbrot_boundary(
	const uint16_t width, const uint16_t height,
	const int32_t *upper, const int32_t *lower, const int32_t *left, const int32_t *right,
	const T& real_min, const T& real_unit,
	const T& imag_min, const T& imag_unit,
	const int32_t iter_max, int32_t *countmap) {

	const uint32_t mapwidth  = width + 1;
	const uint32_t mapheight = height + 1;

	// Initialize COUNTMAP
	for (uint32_t row = 0; row < mapheight; row++) {
		for (uint32_t col = 0; col < mapwidth; col++) {
			countmap[row * mapwidth + col] = INIT;
		}
	}

	// Store the perimeter into COUNTMAP
	for (uint32_t row = 0; row < mapheight; row++) {
		countmap[row * mapwidth] = left[row];
		countmap[row * mapwidth + width] = right[row];
	}
	for (uint32_t col = 0; col < mapwidth; col++) {
		countmap[col] = upper[col];
		countmap[height * mapwidth + col] = lower[col];
	}

	// Search boundaries of perimeter, and push the candidates into QUE.
	// Pixels pushed to QUE are marked as QUEUED on the COUNTMAP.
	std::deque<std::pair<uint16_t, uint16_t>> que;
	// horizontal
	for (uint32_t col = 1; col < width; col++) {
		if (countmap[col] != countmap[col+1]) {
			que.push_back({1, col});
			que.push_back({1, col + 1});
			countmap[mapwidth + col] = QUEUED;
			countmap[mapwidth + col + 1] = QUEUED;
		}
		if (countmap[height * mapwidth + col] != countmap[height * mapwidth + col + 1]) {
			que.push_back({height - 1, col});
			que.push_back({height - 1, col + 1});
			countmap[(height - 1) * mapwidth + col]     = QUEUED;
			countmap[(height - 1) * mapwidth + col + 1] = QUEUED;
		}
	}
	// vertical
	for (uint32_t row = 1; row < height; row++) {
		if (countmap[row * mapwidth] != countmap[(row + 1) * mapwidth]) {
			que.push_back({row, 1});
			que.push_back({row + 1, 1});
			countmap[row * mapwidth + 1]       = QUEUED;
			countmap[(row + 1) * mapwidth + 1] = QUEUED;
		}
		if (countmap[row * mapwidth + width] != countmap[(row + 1) * mapwidth + width]) {
			que.push_back({row, width - 1});
			que.push_back({row + 1, width - 1});
			countmap[row * mapwidth + width - 1]       = QUEUED;
			countmap[(row + 1) * mapwidth + width - 1] = QUEUED;
		}
	}
	// offset of eight-neighbours
	// 0 1 2
	// 3 P 4
	// 5 6 7
	alignas(32) const int32_t offset[8] = {
		-(int32_t)mapwidth - 1, // upper left
		-(int32_t)mapwidth,     // upper
		-(int32_t)mapwidth + 1, // upper right
		-1,                          // left
		1,                           // right
		(int32_t)mapwidth - 1,  // lower left
		(int32_t)mapwidth,      // lower
		(int32_t)mapwidth + 1   // lower right
	};
	constexpr int32_t bound_mask[8] = {
		//76543210
		0b00001010, // upper left,  upper or left
		0b00011101, // upper,       left  or right or upper left  or upper right
		0b00010010, // upper right, upper or right
		0b01100011, // left,        upper or lower or upper left  or lower left
		0b11000110, // right,       upper or lower or upper right or lower right
		0b01001000, // lower left,  lower or left
		0b10111000, // lower,       left  or right or lower left  or lower right
		0b01010000  // lower right, lower or right
	};
#if defined(ENABLE_AVX) and defined(__AVX2__) and defined(ENABLE_GMP)
	// index of the gather instruction
	const __m256i vecoffset = _mm256_loadu_si256((__m256i *)offset);
	while (que.size()) {
		auto&& [row, col] = que.front();
		que.pop_front();
		uint32_t pos = (uint32_t)row * mapwidth + col;
		if (countmap[pos] != QUEUED) continue;
		T real = real_min + real_unit * row;
		T imag = imag_min + imag_unit * col;
		countmap[row * mapwidth + col] = mandelbrot(real, imag, iter_max);
		__m256i neighbour = _mm256_i32gather_epi32(
			countmap + pos, vecoffset, 4);
		__m256i is_init = _mm256_cmpeq_epi32(
			neighbour, _mm256_set1_epi32(INIT));
		__m256i is_queued = _mm256_cmpeq_epi32(
			neighbour, _mm256_set1_epi32(QUEUED));
		__m256i is_unknown = _mm256_or_si256(
			is_init, is_queued);
		uint32_t is_init_shrink = _mm256_movemask_ps(
			_mm256_castsi256_ps(is_init));
		__m256i is_same = _mm256_cmpeq_epi32(
			neighbour, _mm256_set1_epi32(countmap[pos]));
		__m256i notbound = _mm256_or_si256(is_same, is_unknown);
		uint32_t bound_shrink = ~_mm256_movemask_ps(
			_mm256_castsi256_ps(notbound));
		for (uint32_t i = 0; i < 8; i++) {
			if (is_init_shrink & (1 << i)) {
				if (bound_shrink & bound_mask[i]) {
					uint32_t p = pos + offset[i];
					uint16_t r = p / mapwidth;
					uint16_t c = p % mapwidth;
					que.push_back({r, c});
					countmap[p] = QUEUED;
				}
			}
		}
	}
#elif defined(ENABLE_AVX) and defined(__AVX2__) and not defined(ENABLE_GMP)
	static_assert(std::is_same<T, double>::value);
	const __m256i vecoffset = _mm256_loadu_si256((__m256i *)offset);
	alignas(32) uint64_t store[4];
	__m256d real, imag;
	std::array<uint16_t, 4> rpos, cpos;

	while (que.size()) {
		uint32_t size = pop4elms(
			que, countmap, mapwidth, mapheight,
			real_min, real_unit, imag_min, imag_unit,
			real, imag, rpos, cpos);
		__m256i cnt = mandelbrot_avx(real, imag, iter_max);
		_mm256_storeu_si256((__m256i *)store, cnt);
		for (uint32_t i = 0; i < size; i++) {
			auto& row = rpos[i];
			auto& col = cpos[i];
			countmap[row * mapwidth + col] = store[i];
		}
		for (uint32_t i = 0; i < size; i++) {
			auto& row = rpos[i];
			auto& col = cpos[i];
			int32_t pos = (int32_t)row * mapwidth + col;
			__m256i neighbour = _mm256_i32gather_epi32(
				countmap + pos, vecoffset, 4);
			__m256i is_init = _mm256_cmpeq_epi32(
				neighbour, _mm256_set1_epi32(INIT));
			__m256i is_queued = _mm256_cmpeq_epi32(
				neighbour, _mm256_set1_epi32(QUEUED));
			__m256i is_unknown = _mm256_or_si256(
				is_init, is_queued);
			uint32_t is_init_shrink = _mm256_movemask_ps(
				_mm256_castsi256_ps(is_init));
			__m256i is_same = _mm256_cmpeq_epi32(
				neighbour, _mm256_set1_epi32(store[i]));
			__m256i notbound = _mm256_or_si256(is_same, is_unknown);
			uint32_t bound_shrink = ~_mm256_movemask_ps(
				_mm256_castsi256_ps(notbound));
			for (uint32_t j = 0; j < 8; j++) {
				if (is_init_shrink & (1 << j)) {
					if (bound_shrink & bound_mask[j]) {
						uint32_t p = pos + offset[j];
						uint16_t r = p / mapwidth;
						uint16_t c = p % mapwidth;
						que.push_back({r, c});
						countmap[p] = QUEUED;
					}
				}
			}
		}
	}
#else
	while (que.size()) {
		auto&& [row, col] = que.front();
		que.pop_front();
		uint32_t pos = (uint32_t)row * mapwidth + col;
		if (countmap[row * mapwidth + col] >= 0) continue;
		T real = real_min + real_unit * row;
		T imag = imag_min + imag_unit * col;
		countmap[row * mapwidth + col] = mandelbrot(real, imag, iter_max);
		int32_t neighbour[8];
		for (int32_t i = 0; i < 8; i++) {
			neighbour[i] = countmap[pos + offset[i]];
		}
		std::uint8_t is_init = 0;
		for (int32_t i = 0; i < 8; i++) {
			is_init |= (neighbour[i] == INIT ? 1 : 0) << i;
		}
		std::uint8_t is_queued = 0;
		for (int32_t i = 0; i < 8; i++) {
			is_queued |= (neighbour[i] == QUEUED ? 1 : 0) << i;
		}
		std::uint8_t is_same = 0;
		for (int32_t i = 0; i < 8; i++) {
			is_same |= (neighbour[i] == countmap[pos] ? 1 : 0) << i;
		}
		std::int8_t is_unknown = is_init | is_queued;
		std::int8_t is_known = ~is_unknown;
		std::int8_t is_bound = is_known & ~is_same;
		for (uint32_t i = 0; i < 8; i++) {
			if (is_init & (1 << i)) {
				if (is_bound & bound_mask[i]) {
					uint32_t p = pos + offset[i];
					uint16_t r = p / mapwidth;
					uint16_t c = p % mapwidth;
					que.push_back({r, c});
					countmap[p] = QUEUED;
				}
			}
		}
	}
#endif
}

// Compute the escape times linearly and return them as a vector
template<typename T>
void calc_mandelbrot_linear(
	const T& real_min, const T& real_unit,
	const T& imag_min, const T& imag_unit,
	const uint16_t len, const int32_t iter_max,
	std::vector<int32_t>& result) {

	for (auto&& val: result) {
		val = INIT;
	}

#if defined(ENABLE_AVX) and defined(__AVX2__) and not defined(ENABLE_GMP)
	static_assert(std::is_same<T, double>::value);
	const uint32_t floored_len = len - (len % 4);
	alignas(32) uint64_t store[4];
	for (uint32_t i = 0; i < floored_len; i += 4) {
		double real = real_min + real_unit * i;
		double imag = imag_min + imag_unit * i;
		__m256d a = _mm256_set_pd(
			real + real_unit * 3,
			real + real_unit * 2,
			real + real_unit * 1,
			real + real_unit * 0);
		__m256d b = _mm256_set_pd(
			imag + imag_unit * 3,
			imag + imag_unit * 2,
			imag + imag_unit * 1,
			imag + imag_unit * 0);
		__m256i cnt = mandelbrot_avx(a, b, iter_max);
		_mm256_storeu_si256((__m256i *)store, cnt);
		for (uint32_t j = 0; j < 4; j++) {
			result[i + j] = store[j];
		}
	}
	for (uint32_t i = floored_len; i < len; i++) {
		double real = real_min + real_unit * i;
		double imag = imag_min + imag_unit * i;
		result[i] = mandelbrot(real, imag, iter_max);
	}
#else
	for (uint32_t i = 0; i < len; i++) {
		T real = real_min + real_unit * i;
		T imag = imag_min + imag_unit * i;
		result[i] = mandelbrot(real, imag, iter_max);
	}
#endif
}

// Compute Mandelbrot set escape time boundary
template<typename T>
void calc_mandelbrot_boundary(
	const uint16_t width, const uint16_t height,
	const T& real_min, const T& real_max,
	const T& imag_min, const T& imag_max,
	const int32_t iter_max, region_manager& region) {

	const T& real_range = real_max - real_min;
	const T& imag_range = imag_max - imag_min;
	const T& real_unit  = real_range / height;
	const T& imag_unit  = imag_range / width;

	const uint16_t row_split = region.row_split;
	const uint16_t col_split = region.col_split;

	auto& rows = region.rows;
	auto& cols = region.cols;
	auto& row_perimeters = region.row_perimeters;
	auto& col_perimeters = region.col_perimeters;
	auto& subcountmaps = region.subcountmaps;
	auto& countmap = region.countmap;
	
	auto start = std::chrono::system_clock::now();

	// Compute the perimeters
	// (real_min, imag_min)---(real_min, imag_max)
	//          |                      |
	//          |                      |
	// (real_max, imag_min)---(real_max, imag_max)
#pragma omp parallel for
	for (uint32_t i = 0; i <= row_split; i++) {
		T real = real_min + real_unit * rows[i];
		calc_mandelbrot_linear(real, T(0.0), imag_min, imag_unit, width + 1, iter_max, row_perimeters[i]);
	}
#pragma omp parallel for
	for (uint32_t i = 0; i <= col_split; i++) {
		T imag = imag_min + imag_unit * cols[i];
		calc_mandelbrot_linear(real_min, real_unit, imag, T(0.0), height + 1, iter_max, col_perimeters[i]);
	}
	// Compute submandelbrot countmaps
#pragma omp parallel for collapse(2)
	for (uint32_t prow = 0; prow < row_split; prow++) {
		for (uint32_t pcol = 0; pcol < col_split; pcol++) {
			uint16_t mapwidth = cols[pcol + 1] - cols[pcol];
			uint16_t mapheight = rows[prow + 1] - rows[prow];
			calc_submandelbrot_boundary(
				mapwidth, mapheight,
				row_perimeters[prow].data() + cols[pcol], row_perimeters[prow + 1].data() + cols[pcol],
				col_perimeters[pcol].data() + rows[prow], col_perimeters[pcol + 1].data() + rows[prow],
				T(real_min + real_unit * rows[prow]), real_unit,
				T(imag_min + real_unit * cols[pcol]), imag_unit, iter_max, subcountmaps[prow][pcol]);
		}
	}

	auto end = std::chrono::system_clock::now();
	double elapsed = (double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() / 1000;
	std::cerr << elapsed << std::endl;

	// Store countmaps into 2-dimentional vector RESULT.
	for (uint32_t i = 0; i < row_split; i++) {
		for (uint32_t j = 0; j < col_split; j++) {
			for (uint32_t row = rows[i]; row < rows[i + 1]; row++) {
				for (uint32_t col = cols[j]; col < cols[j + 1]; col++) {
					uint32_t r = row - rows[i];
					uint32_t c = col - cols[j];
					uint32_t w = cols[j + 1] - cols[j] + 1;
#if defined(FILL_COUNTMAP)
					if (subcountmaps[i][j][r * w + c] == INIT) {
						countmap[row][col] = countmap[row][col - 1];
					}else {
						countmap[row][col] = subcountmaps[i][j][r * w + c];
					}
#else
					countmap[row][col] = subcountmaps[i][j][r * w + c];
#endif
				}
			}
		}
	}
}

region_manager::region_manager(
	uint16_t width, uint16_t height,
	uint16_t col_split, uint16_t row_split) :
	width(width), height(height), col_split(col_split), row_split(row_split) {

	this->countmap = std::vector<std::vector<int32_t>>(
		this->height, std::vector<int32_t>(
			this->width));
	this->rows = std::vector<uint16_t>(this->row_split + 1);
	this->rows[0] = 0;
	for (uint32_t i = 1; i <= this->row_split; i++) {
		this->rows[i] = this->rows[i - 1] + ((i < this->height % this->row_split) ? this->height / this->row_split + 1 : this->height / this->row_split);
	}

	this->cols = std::vector<uint16_t>(this->col_split + 1);
	this->cols[0] = 0;
	for (uint32_t i = 1; i <= this->col_split; i++) {
		this->cols[i] = this->cols[i - 1] + ((i < this->width % this->col_split) ? this->width / this->col_split + 1 : this->width / this->col_split) ;
	}

	this->subcountmaps = std::vector<std::vector<int32_t *>>(
		this->row_split, std::vector<int32_t *>(
			this->col_split));
	for (uint32_t row = 0; row < this->row_split; row++) {
		for (uint32_t col = 0; col < this->col_split; col++) {
			uint16_t mapheight = this->rows[row + 1] - this->rows[row] + 1;
			uint16_t mapwidth = this->cols[col + 1] - this->cols[col] + 1;
			uint32_t maparea = mapwidth * mapheight;
			this->subcountmaps[row][col] = (int32_t *)std::aligned_alloc(32, sizeof(int32_t) * maparea);
		}
	}
	this->row_perimeters = std::vector<std::vector<int32_t>>(
		row_split + 1, std::vector<int32_t>(
			width + 1));
	this->col_perimeters = std::vector<std::vector<int32_t>>(
		col_split + 1, std::vector<int32_t>(
			height + 1));
}

region_manager::region_manager(const region_manager& region) :
	width(region.width), height(region.height), col_split(region.col_split), row_split(region.row_split) {
	this->countmap = region.countmap;
	this->cols = region.cols;
	this->rows = region.rows;
	this->col_perimeters = region.col_perimeters;
	this->row_perimeters = region.row_perimeters;
	this->subcountmaps = std::vector<std::vector<int32_t *>>(
		this->row_split, std::vector<int32_t *>(
			this->col_split));
	for (uint32_t row = 0; row < this->row_split; row++) {
		for (uint32_t col = 0; col < this->col_split; col++) {
			uint16_t mapheight = this->rows[row + 1] - this->rows[row] + 1;
			uint16_t mapwidth = this->cols[col + 1] - this->cols[col] + 1;
			uint32_t maparea = mapwidth * mapheight;
			this->subcountmaps[row][col] = (int32_t *)std::aligned_alloc(32, sizeof(int32_t) * maparea);
			for (uint32_t i = 0; i < maparea; i++) {
				this->subcountmaps[row][col][i] = region.subcountmaps[row][col][i];
			}
		}
	}
}

region_manager::~region_manager() {
	for (uint32_t row = 0; row < this->row_split; row++) {
		for (uint32_t col = 0; col < this->col_split; col++) {
			std::free(this->subcountmaps[row][col]);
		}
	}
}

#if defined(ENABLE_GMP)
template std::vector<std::vector<int32_t>> calc_mandelbrot_countmap (
	const uint16_t, const uint16_t,
	const mpf_class&, const mpf_class&, const mpf_class&, const mpf_class&,
	const int32_t);

template void calc_mandelbrot_boundary<mpf_class> (
	const uint16_t, const uint16_t,
	const mpf_class&, const mpf_class&, const mpf_class&, const mpf_class&,
	const int32_t, region_manager&);
#endif

template std::vector<std::vector<int32_t>> calc_mandelbrot_countmap (
	const uint16_t, const uint16_t,
	const double&, const double&, const double&, const double&,
	const int32_t);

template void calc_mandelbrot_boundary<double> (
	const uint16_t, const uint16_t,
	const double&, const double&, const double&, const double&,
	const int32_t, region_manager&);

