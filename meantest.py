import os
import numpy as np

data = [[1,2],[2,3],[4,5],[6,7]]
for i in range(len(data[0])):
	print np.mean(data[:, i])
