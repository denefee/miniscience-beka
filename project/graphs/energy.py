import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("log.csv")
Pot_E = df.PotEng
Tot_E = df.TotEng
Kin_E = df.KinEng
T = df.Time

plt.xlabel(r'Time, seconds', fontsize=14)
plt.ylabel(r'Energy, ergs', fontsize=14)
plt.grid(True)
plt.plot(T, Pot_E, color="darkblue", label=r'Potential energy')
plt.plot(T, Tot_E, color="black", label=r'Total energy')
plt.plot(T, Kin_E, color="firebrick", label=r'Kinetic energy')
plt.legend(loc='best', fontsize=12)
# plt.savefig("eng.png")
plt.show()
