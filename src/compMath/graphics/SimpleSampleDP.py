import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_simple_sample_dp.txt')

plt.scatter(data[0, :], data[1, :], label=f'Численное решение y(x),\nКоличество узлов {data.shape[1]}', s=0.1)
plt.scatter(data[0, :], data[2, :], label=f'Численное решение y\'(x),\nКоличество узлов {data.shape[1]}', s=0.1)
plt.legend()
plt.title('Дорманд Принс')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/graphic_simple_sample_dp.png')
plt.show()
