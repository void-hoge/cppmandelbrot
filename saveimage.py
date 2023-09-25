#!/usr/bin/env python3

import sys
import os
import numpy as np
from PIL import Image

def render():
    cmd = f'cat {sys.argv[1]}'
    output = os.popen(cmd).read()
    print('converting to image')
    lines = output.strip().split('\n')
    width, height = map(int, lines[0].strip().split())
    grayscale = []
    for line in lines[1:]:
        col = line.strip().split()
        grayscale.extend(map(lambda x: int(x)%256, col))
    arr = np.array(grayscale, dtype=np.uint8)
    arr = np.reshape(arr, (height, width))
    image = Image.fromarray(arr, mode='L')
    image.save(f'{sys.argv[2]}')

if __name__ == '__main__':
    render()
