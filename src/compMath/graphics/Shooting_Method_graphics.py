import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_shooting_method.txt')

data2 = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_check_shooting_method.txt')

plt.scatter(data[0, :], data[1, :], label=f'Численное решение методом стрельбы,\nКоличество узлов {data.shape[1]}', s=0.5)

plt.scatter(data2[0, :], data2[1, :], label=f'Численное решение методом прогонки,\nКоличество узлов {data.shape[1]}', s=0.5)
plt.legend()
plt.title('Метод стрельбы')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/test_shooting_method.png')
plt.show()



