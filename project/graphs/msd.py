import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("log.csv")
t = df.Time
r1 = df.c_msd1
r2 = df.c_msd2
r3 = df.c_msd3
r4 = df.c_msd4
r5 = df.c_msd5

plt.xlabel(r'Time, seconds', fontsize=14)
plt.ylabel(r'MSD, $cm^2$', fontsize=14)
plt.grid(True)
plt.plot(t, r1, color="firebrick", label=r'1st region')
plt.plot(t, r2, color="darkorange", label=r'2nd region')
plt.plot(t, r3, color="darkgreen", label=r'3rd region')
plt.plot(t, r4, color="darkblue", label=r'4th region')
plt.plot(t, r5, color="darkmagenta", label=r'5th region')
plt.legend(loc='best', fontsize=12)
# plt.savefig("eng.png")
plt.show()
