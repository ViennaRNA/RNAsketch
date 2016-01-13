import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('data.csv', delimiter=';', skip_header=1)

plt.rc('font', family='Liberation Sans', size='18')
plt.grid(True)
x = [ x*100 for x in data[:,3] ]
y1 = [ y1*100 for y1 in data[:,4] ]
y2 = [ y2*100 for y2 in data[:,5] ]
plt.plot(x, y1, 'o-', markersize=12, markeredgewidth=0, linewidth=3.0, label='Fair Sampling')
plt.plot(x, y2, '^-', markersize=12, markeredgewidth=0, linewidth=3.0, label='Unfair Sampling')
plt.legend(loc=0, fontsize=18)
plt.ylabel('Unique solutions [%]')
plt.xlabel('Sample size relative to size of solution space [%]')
plt.savefig('xy.svg')
