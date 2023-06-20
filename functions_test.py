import numpy as np
import AtmosModel as At
import matplotlib.pyplot as plt

QN = np.arange(-3, 3, 0.01)

tau_c = [At.get_scale_time_condesation(1, 0, i, 1) for i in QN]
tau_c = np.array(tau_c)
# plt.plot(1 / tau_c)
# plt.show()

QR = np.arange(-3, 3, 0.01)

VT = [At.get_terminalvelocity(1, np.exp(i), 1) for i in QR]
VT = np.array(VT)
# plt.plot(VT)
# plt.show()

VTN = [At.get_aerosolvelocity(1, 1, i, 1) for i in QR]
VTN = np.array(VTN)
plt.plot(VTN)
plt.show()
