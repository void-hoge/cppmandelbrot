#include "../gui/gui.hpp"

#include <sstream>

int main(const int argc, const char *argv[]) {
	if (argc < 2) {
		std::stringstream ss;
		ss << "usage: " << argv[0] << " <resolution>" << std::endl;
		throw std::invalid_argument(ss.str());
	}
	std::string resolution = argv[1];
	size_t xpos = resolution.find('x');
	if (xpos == std::string::npos) {
		std::stringstream ss;
		ss << "<resolution> must be formmated as <width>x<height>" << std::endl;
		throw std::invalid_argument(ss.str());
	}
	int width = std::stoi(resolution.substr(0, xpos));
	int height = std::stoi(resolution.substr(xpos + 1));
	auto gui = MandelbrotGUI(width, height, std::cerr);
	gui.mainloop();
}
