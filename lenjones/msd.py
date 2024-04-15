import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("msd.csv")
T = np.asarray(df['T'])
MSD = np.asarray(df['MSD'])

plt.xlabel(r'time', fontsize=14)
plt.ylabel(r'$\Delta R^2$', fontsize=14)
plt.title(r'$\Delta R^2 (t)$', fontsize=14)
plt.grid(True)
plt.scatter(T, MSD, color="firebrick", label=r'')
# plt.savefig("density"+str(i)+".png")
# plt.clf()
y = MSD
dy = 0
x = T
dx = 0

pol = np.polyfit(x[1500:], y[1500:], deg=1, cov=True)
print(pol)

k = pol[0][0]
b = pol[0][1]
dk = pol[1][0][0]
db = pol[1][1][1]

pol = pol[0]
yfit = np.polyval(pol, x)
plt.plot(x, yfit, color="darkgreen", label=r'k = '+str(k)+' $\pm$ '+str(dk)+'\n'+'b = '+str(b)+' $\pm$ '+str(db))


plt.show()