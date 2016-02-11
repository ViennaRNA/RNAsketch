import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('data.csv', delimiter=';', skip_header=1)
fig = plt.figure()
fig.set_size_inches(14, 7)
ax = fig.add_subplot(111)
plt.rc('font', family='Liberation Sans', size='18')
ax.grid(True)
x = [ x*100 for x in data[:,3] ]
y1 = [ y1*100 for y1 in data[:,4] ]
y2 = [ y2*100 for y2 in data[:,5] ]
ax.plot(x, y1, 'o-', markersize=12, markeredgewidth=0, linewidth=3.0, color='#006DDB', label='Fair Sampling')
ax.plot(x, y2, '^-', markersize=12, markeredgewidth=0, linewidth=3.0, color='#FF8000', label='Unfair Sampling')
ax.legend(loc=0, fontsize=18)
ax.set_ylabel('Unique solutions [%]')
ax.set_xlabel('Sample size relative to size of solution space [%]')
fig.savefig('xy.svg')
