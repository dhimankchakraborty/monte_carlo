import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys



try:
    total_points_count = int(sys.argv[1])
except ValueError:
    print("Input should be an integer.")
    exit()

rng = np.random.default_rng()

points_count_inside_circle = 0
for i in tqdm(range(total_points_count)):
    if (rng.random()**2 + rng.random()**2) <= 1:
        points_count_inside_circle += 1

estimated_pi = 4 * np.float128(points_count_inside_circle) / np.float128(total_points_count)
print(f"Estimated value of Pi: {estimated_pi}")
print(f"Error in the esimated value of Pi: {((estimated_pi - np.pi) / np.pi * 100):.5f} %")