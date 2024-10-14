#include "gui.hpp"

MandelbrotGUI::MandelbrotGUI(const uint32_t width, const uint32_t height, const uint32_t num_threads, std::ostream& ost) : buffersize({width, height}), region(width, height, num_threads, 1), ost(ost) {
	this->buffer = std::vector<unsigned char>(width * height * 3, 0);
	this->update = true;
	this->windowsize = {width, height};

#if defined(ENABLE_GMP)
	this->center = {mpf_class(-0.5), mpf_class(0.0)};
	this->range = {mpf_class(2.0), mpf_class(2.0) * ((double)width/ height)};
#else
	this->center = {-0.5, 0.0};
	this->range = {2.0, 2.0 * ((double)width / height)};
#endif
	this->iter_max = 1<<10;
	this->zoomexp = 1;

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	this->window = SDL_CreateWindow(
		"MandelbrotGUI",
		SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
		this->windowsize.first, this->windowsize.second, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
	if (!this->window) {
		throw std::runtime_error("Failed to init SDL window.");
	}
	SDL_SetWindowResizable(this->window, SDL_TRUE);
	this->context = SDL_GL_CreateContext(this->window);
	if (!this->context) {
		throw std::runtime_error("Failed to init SDL OpenGL contest.");
	}
	this->renderer = SDL_CreateRenderer(this->window, -1, SDL_RENDERER_ACCELERATED);
	if (!this->renderer) {
		throw std::runtime_error("Failed to init SDL renderer.");
	}
	glViewport(0, 0, this->windowsize.first, this->windowsize.second);
}

MandelbrotGUI::~MandelbrotGUI() {
	SDL_QuitSubSystem(SDL_INIT_VIDEO);
	SDL_DestroyRenderer(this->renderer);
	SDL_GL_DeleteContext(this->context);
	SDL_Quit();
}

void MandelbrotGUI::handle_window_events(const SDL_Event& event) {
	if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
		this->update = true;
		this->windowsize = {event.window.data1, event.window.data2};
		glViewport(0, 0, this->windowsize.first, this->windowsize.second);
		glPixelZoom((float)this->windowsize.first / this->buffersize.first,
					(float)this->windowsize.second / this->buffersize.second);
	}else if (event.window.event == SDL_WINDOWEVENT_EXPOSED) {
		this->update = true;
	}
}

void MandelbrotGUI::handle_mouse_press(const SDL_Event& event) {
	if (event.button.button == SDL_BUTTON_LEFT || event.button.button == SDL_BUTTON_RIGHT) {
		this->update = true;
		this->reload = true;
		const double real_relative = (double)event.button.y / this->windowsize.second;
		const double imag_relative = (double)event.button.x / this->windowsize.first;
		const auto& [real_center, imag_center] = this->center;
		const auto& [real_range, imag_range] = this->range;
#if defined(ENABLE_GMP)
		const mpf_class real_min = real_center - real_range / 2.0;
		const mpf_class imag_min = imag_center - imag_range / 2.0;
#else
		const auto real_min = real_center - real_range / 2.0;
		const auto imag_min = imag_center - imag_range / 2.0;
#endif
		this->center = {real_min + real_range * real_relative, imag_min + imag_range * imag_relative};
		if (event.button.button == SDL_BUTTON_LEFT) {
			this->range = {real_range / 2.0, imag_range / 2.0};
			this->iter_max += 1<<10;
			this->zoomexp += 1;
		}else {
			this->range = {real_range * 2.0, imag_range * 2.0};
			this->iter_max = (this->iter_max == (1 << 10)) ? (1 << 10) : (this->iter_max - (1 << 10));
			this->zoomexp -= 1;
		}
	}
}

void MandelbrotGUI::render_mandelbrot() {
	if (this->reload && this->update) {
		const auto& [width, height] = this->buffersize;
		const auto [real_center, imag_center] = this->center;
		const auto [real_range, imag_range] = this->range;
#if defined(ENABLE_GMP)
		const mpf_class real_min = real_center - real_range / 2;
		const mpf_class real_max = real_center + real_range / 2;
		const mpf_class imag_min = imag_center - imag_range / 2;
		const mpf_class imag_max = imag_center + imag_range / 2;
#else
		const auto real_min = real_center - real_range / 2;
		const auto real_max = real_center + real_range / 2;
		const auto imag_min = imag_center - imag_range / 2;
		const auto imag_max = imag_center + imag_range / 2;
#endif
		this->ost << "rendering..." << std::endl;
		calc_mandelbrot_boundary(width, height, real_min, real_max, imag_min, imag_max, iter_max, this->region, this->ost);
		this->ost << std::fixed << std::setprecision(std::ceil(this->zoomexp / std::log2(10)) + 5) << "real, imag = " << this->center.first << ", " << this->center.second << std::endl;
		this->ost << "height (real) = " << this->range.first << std::endl;
		this->ost << "iter_max = " << this->iter_max << std::endl;
		this->ost << std::defaultfloat << std::setprecision(std::cout.precision());
		for (int row = 0; row < this->buffersize.second; row++) {
			for (int col = 0; col < this->buffersize.first; col++) {
				auto grayscale = this->region.countmap[this->buffersize.second - row - 1][col] % 255;
				this->buffer[(row * this->buffersize.first + col) * 3 + 0] = grayscale;
				this->buffer[(row * this->buffersize.first + col) * 3 + 1] = grayscale;
				this->buffer[(row * this->buffersize.first + col) * 3 + 2] = grayscale;
			}
		}
		this->update = false;
		this->reload = false;
	}
}

void MandelbrotGUI::mainloop() {
	this->update = true;
	this->reload = true;
	while (true) {
		auto framestart = SDL_GetTicks();
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_QUIT) {
				return;
			}else if (event.type == SDL_MOUSEBUTTONDOWN) {
				this->handle_mouse_press(event);
			}else if (event.type == SDL_WINDOWEVENT) {
				this->handle_window_events(event);
			}else if (event.type == SDL_KEYDOWN) {
				if (event.key.keysym.sym == SDLK_q) {
					return;
				}
			}
		}
		if (this->update) {
			glClear(GL_COLOR_BUFFER_BIT);
			this->render_mandelbrot();
			glDrawPixels(this->buffersize.first, this->buffersize.second, GL_RGB, GL_UNSIGNED_BYTE, this->buffer.data());
			SDL_GL_SwapWindow(this->window);
			this->update = false;
		}
		auto frametime = SDL_GetTicks() - framestart;
		if (frametime < FRAME_DELAY) {
			SDL_Delay(FRAME_DELAY - frametime);
		}
	}
}
