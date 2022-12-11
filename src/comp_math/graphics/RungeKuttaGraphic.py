import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_rayleigh_rk.txt')

plt.scatter(data[1, :], data[0, :], label=f'Численное решение,\nКоличество узлов {data.shape[1]}', s=0.1)
plt.legend(fontsize=16)
plt.title('Рунге-Кутта 4')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/graphic_rayleigh_rk.png')
plt.show()
