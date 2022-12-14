import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_three_bodies_dp7.txt')
n = 20000
plt.scatter(data[0, :n], data[1, :n], label=f'Численное решение,\nКоличество узлов {n}', s=0.1)

plt.legend()
plt.title('Дорман Принц с контролем шага')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/test_three_bodies_dp7.png')
plt.show()
