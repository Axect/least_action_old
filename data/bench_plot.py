import csv
import numpy as np
import pylab as plt

x = [50*i for i in range(1, 11)]
fb_1d_bf = []

multiplier = 2**(20)-1
x2 = [i*multiplier for i in range(1, 11)]
fb_1d_dc = []

with open('fb_1d_bf.csv', 'r') as fr:
    rdr = csv.reader(fr)
    next(rdr, None)
    for line in rdr:
        fb_1d_bf.append(float(line[1]))

with open('fb_1d_dc.csv', 'r') as fr2:
    rdr2 = csv.reader(fr2)
    next(rdr2, None)
    for line in rdr2:
        fb_1d_dc.append(float(line[1]))

fb_1d_bf_coef = np.polyfit(x, fb_1d_bf, 3)
fb_1d_bf_poly = np.poly1d(fb_1d_bf_coef)

fb_1d_dc_coef = np.polyfit(x2, fb_1d_dc, 1)
fb_1d_dc_poly = np.poly1d(fb_1d_dc_coef)

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot
plt.figure(figsize=(10,6), dpi=300)

# plt.title(r"Free Body 1D Bruteforce", fontsize=16)
# plt.xlabel(r'$N$', fontsize=14)
# plt.ylabel(r'time (sec)', fontsize=14)
# plt.plot(x, fb_1d_bf, 'o-', label=r'Bruteforce $m=3$', alpha=0.5)
# plt.plot(x, fb_1d_bf_poly(x), label=r'polyfit: order 3', alpha=0.5)

plt.title(r"Free Body 1D Divide and Conquer", fontsize=16)
plt.xlabel(r'$N$', fontsize=14)
plt.ylabel(r'time (sec)', fontsize=14)
plt.plot(x2, fb_1d_dc, 'o-', label=r'Divide and Conquer $m=7$', alpha=0.5)
plt.plot(x2, fb_1d_dc_poly(x2), label=r'polyfit: order 1', alpha=0.5)

plt.legend(fontsize=12)
plt.grid()
plt.savefig("fb_1d_dc.png", dpi=300)
