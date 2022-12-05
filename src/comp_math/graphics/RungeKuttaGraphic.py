import matplotlib.pyplot as plt
import numpy as np
plt.figure()
fig, ax = plt.subplots(figsize=(12,8))
data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_rayleigh.txt')

ax.plot(data[:, 1], data[:, 0], 'x', label=f'Численное решение,\nКоличество узлов {data.shape[0]}', markersize=2)
plt.legend(fontsize=16)
plt.title('Рунге-Кутта 4')
ax.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('graphic2_yx3.png')
plt.show()