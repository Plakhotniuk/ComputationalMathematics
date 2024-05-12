import matplotlib.pyplot as plt
import numpy as np
import imageio

data = np.loadtxt(f'/Users/arseniy/study/CompMath/ComputationalMathematics/tests/data_files/twoDimEqThermalCond.txt')
def create_frame(t):
    fig = plt.figure(figsize=(6, 6))
    plt.imshow(data[data.shape[1]*(t):data.shape[1]*(t+1)])
    plt.colorbar()
    plt.title(f'step {t}',
              fontsize=14)

    plt.savefig(f'img_{t}.png',
                transparent = False,
                facecolor = 'white'
                )
    plt.close()


time = [i for i in range(0, int(data.shape[0]/data.shape[1]))]


for t in time:
    create_frame(t)

frames = []
for t in time:
    image = imageio.v2.imread(f'img_{t}.png')
    frames.append(image)


imageio.mimsave('./example.gif', # output gif
                frames,          # array of input frames
                fps = 5)         # optional: frames per second

print(data.shape)

print(data[0:data.shape[1]] - data[-data.shape[1]:])

# plt.imshow(data)
# plt.show()
