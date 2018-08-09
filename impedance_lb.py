import numpy as np

R = np.array([10.2, 40.1, 200, 990, 40000])
err_r = np.repeat(0.3, len(R))
V = np.array([0.111, 0.420, 1.760, 4.91, 8.72])
err_v = np.array([0.001, 0.001, 0.001, 0.01, 0.01])
I = V/R
err_cur = ((err_v**2 + (I * err_r)**2)**0.5)/R
