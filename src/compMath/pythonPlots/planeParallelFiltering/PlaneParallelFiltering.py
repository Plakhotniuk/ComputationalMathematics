import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/study/CompMath/ComputationalMathematics/tests/data_files/planeParallelFiltering.txt')

plt.scatter(data[1, :], data[0, :], label=f'Численное решение,\nКоличество узлов {data.shape[1]}')
plt.legend()
plt.title('Плоскопараллельная однофазная фильтрация')
plt.grid()
plt.xlabel('x')
plt.ylabel('y(x)')
plt.savefig('planeParallelFiltering.png')
plt.show()