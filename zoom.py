#!/usr/bin/env python3

import sys
import os
import numpy as np
from PIL import Image

def render(filename, zoom):
    cmd = f'echo {zoom} | ./{sys.argv[1]}'
    output = os.popen(cmd).read()
    lines = output.strip().split('\n')
    width, height = map(int, lines[0].strip().split())
    grayscale = []
    for line in lines[1:]:
        col = line.strip().split()
        grayscale.extend(map(lambda x: int(x)%256, col))
    arr = np.array(grayscale, dtype=np.uint8)
    arr = np.reshape(arr, (height, width))
    image = Image.fromarray(arr, mode='L')
    image.save(filename)

if __name__ == '__main__':
    zoom = 3
    for i in range(0, 100):
        render(f'mandelbrot{i:02}.png', zoom * 0.8 ** i)
