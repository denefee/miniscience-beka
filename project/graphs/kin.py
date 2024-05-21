import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
df = pd.read_csv("log.csv")
t = df.Time
r1 = df.c_avke1
r2 = df.c_avke2
r3 = df.c_avke3
r4 = df.c_avke4
r5 = df.c_avke5

plt.xlabel(r'Time, seconds', fontsize=14)
plt.ylabel(r'Average kinetic energy, ergs', fontsize=14)
plt.grid(True)
plt.plot(t, r1, color="firebrick", label=r'1st region')
plt.plot(t, r2, color="darkorange", label=r'2nd region')
plt.plot(t, r3, color="darkgreen", label=r'3rd region')
plt.plot(t, r4, color="darkblue", label=r'4th region')
plt.plot(t, r5, color="darkmagenta", label=r'5th region')
plt.legend(loc='best', fontsize=12)
# plt.savefig("eng.png")

r1 = r1.mean()
r2 = r2.mean()
r3 = r3.mean()
r4 = r4.mean()
r5 = r5.mean()

print(r1, r2, r3, r4, r5)
plt.show()
