#pragma once

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>

#include <SDL2/SDL.h>
#include <SDL_opengl.h>

#include "../mandelbrot/mandelbrot.hpp"

#define FRAME_RATE 60
#define FRAME_DELAY (1000 / FRAME_RATE)

class MandelbrotGUI {
private:
	std::ostream& ost;
	
	SDL_Window *window;
	SDL_GLContext context;
	SDL_Renderer *renderer;
	bool update;
	bool reload;

	std::pair<uint32_t, uint32_t> windowsize;
	const std::pair<uint32_t, uint32_t> buffersize;
	std::vector<unsigned char> buffer;

	region_manager region;
	std::pair<double, double> center;
	std::pair<double, double> range;
	uint32_t iter_max;

	void handle_window_events(const SDL_Event& event);
	void handle_mouse_press(const SDL_Event& event);
	void render_mandelbrot();
public:
	MandelbrotGUI(const uint32_t width = 1920, const uint32_t height = 1080, std::ostream& ost = std::cerr);
	MandelbrotGUI(const MandelbrotGUI&) = delete;
	MandelbrotGUI& operator=(const MandelbrotGUI&) = delete;
	~MandelbrotGUI();
	void mainloop();
};
