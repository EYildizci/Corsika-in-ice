#!/usr/bin/env python
import pandas as pd
import matplotlib
import numpy as np
matplotlib.use("agg")
import matplotlib.pyplot as plt

print 'run: echo "id energy px py pz x y" > output2 &&  grep "  [56] " output >> output2'
df = pd.read_csv("output2", sep=r"\s*")
mmin = df.energy.min()
mmean = df.energy.mean()
mmax = df.energy.max()
fig = plt.figure()
ax = fig.gca()
bins = np.logspace(np.log10(mmin*0.9), np.log10(mmax*1.1), 50)
df.energy.hist(bins=bins, ax=ax)
ax.semilogx()
ax.set_xlabel("mu energy [GeV]")
print "mean = ", mmean, mmin , mmax
fig.savefig("output2.png")
