import numpy as np
import sys
import matplotlib.pyplot as plt



if len(sys.argv) < 3:
    print("Usage: python plotter.py file_to_plot [file_to_save_to]")

x = np.fromfile(sys.argv[1], np.float64, -1)
n = int(np.sqrt(x.shape[0]))
x = np.flipud(x.reshape((n,n)))

plt.imshow(x, cmap="gray")
plt.show()
