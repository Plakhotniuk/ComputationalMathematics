import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_nonlin_shooting_method.txt')

plt.scatter(data[0, :], data[1, :], label=f'Численное решение методом стрельбы,\nКоличество узлов {data.shape[1]}', s=0.1)

plt.legend()
plt.title('Метод стрельбы')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/test_nonlin_shooting_method.png')
plt.show()