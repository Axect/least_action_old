import matplotlib.pyplot as plt
import numpy as np

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot
plt.figure(figsize=(10,10), dpi=300)
plt.title(r"2D Paths of Algorithmic Least Action", fontsize=16)
plt.xlabel(r'$x$', fontsize=14)
plt.ylabel(r'$y$', fontsize=14)

x1 = [5*i for i in range(0, 5)]
y2 = [0,15,20,15,0]
x3 = [0,6,12,17,20]
y3 = [20,17,12,6,0]
plt.plot(x1,x1, 'o-', label=r"Free Body $m=3$", alpha=0.5)
plt.plot(x1,y2, 'o-', label=r"Uniform Gravity $m=3$", alpha=0.5)
plt.plot(x3,y3, 'o-', label=r"General Gravity $m=3$", alpha=0.5)

x = np.array([0.1*i for i in range(0, 201)])
ug = 4*x - 0.2*x**2
plt.plot(x, ug, label=r"Uniform Gravity (real)", alpha=0.5)


plt.legend(fontsize=12)
plt.grid()
plt.savefig("motion_2d.png", dpi=300)
