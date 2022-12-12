import matplotlib.pyplot as plt
import numpy as np
plt.figure()

data = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_nonlin_bvp.txt')
data2 = np.loadtxt(f'/Users/arseniy/Desktop/CompMath/ComputationalMathematics/cmake-build-debug/tests/test_nonlin_shooting_method.txt')

plt.plot(data[1, :], data[0, :], label=f'Численное решение методом квазилинеаризации,\nКоличество узлов {data.shape[1]}', color="red")
plt.scatter(data[1, :], data[0, :], label=f'Численное решение методом методом стрельб,\nКоличество узлов {data2.shape[1]}', s=3)

plt.legend()
plt.title('Метод стрельбы')
plt.grid()
plt.xlabel('x', fontsize=14)
plt.ylabel('y(x)', fontsize=14)
plt.savefig('data/test_nonlin_shooting_method.png')
plt.show()