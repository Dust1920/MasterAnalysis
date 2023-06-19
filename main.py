import numpy as np
import AtmosModel as At
import json

with open("parameters.json",'r') as f:
 Parameters = json.load(f)

print(Parameters)
# Escalas
Ls = 10  # Km
Ths = 3  # K
Ts = 15  # min
Qs = 1000  # g/Kg

Vs = (Ls * 1000) / (Ts * 60)

# Condiciones Tiempo-Espacio
z0 = 0.5
z0 = z0 / Ls
dZ0 = 0.1
dZ0 = dZ0 / Ls

t0 = 0
t0 = t0 / Ts
dT0 = 0.15
dT0 = dT0 / Ts

# Condiciones iniciales
W10 = 0  # Velocidad.
W10 = W10 / Vs
W11 = 0
W11 = W11 / Vs

T10 = 300  # Temperatura potencial.
T10 = T10 / Ths
T11 = 300
T11 = T11 / Ths

QV10 = 16  # Proporción de vapror de agua
QV10 = QV10 / Qs
QV11 = 16
QV11 = QV11 / Qs

QR10 = 0  # Proporción de agua líquida
QR10 = QR10 / Qs
QR11 = 0
QR11 = QR11 / Qs

QN10 = 0  # Proporción de núcleos de condensación
QN10 = QN10 / Qs
QN11 = 0
QN11 = QN11 / Qs

# Zona de trabajo
Wblock = np.zeros((2, 2))
Tblock = np.zeros((2, 2))
QVblock = np.zeros((2, 2))
QRblock = np.zeros((2, 2))
QNblock = np.zeros((2, 2))

# Parametros
B = 3  # K
B = B * Ls / Ths

L = 2.5e6
cp = 1e3
cp = cp * Ths
LCP = L / cp

q_star = 1000
q_star = q_star / Qs
qvs0 = 28
qvs0 = qvs0/ Qs
vt0 = 1
vtn0 = 1
epsilon = 0.6
ThetaE = (T11 + B * z0) + LCP * At.approxfqv(z0, QV11)

tau_w = 7.5  # min
tau_w = tau_w / Ts

g = 9.81  # m/s^2
g = g / Vs * 60 * Ts

b_w = 0.491  # m/s sqrt(s)
b_w = b_w / (Ls * 1000) * (Ts * 60) * np.sqrt(Ts * 60)

tau0 = 0.1
gamma = 0.3
taue = 0.2
c = 0.2
c0 = 0.4

vpar = np.array([B, QV11, LCP, epsilon, ThetaE, g])

# El programa recupera el trabajo de la Tesis de Licenciatura
# qv0 = 16
# qv0 = qv0/Qs
# qvs0 = 22
# qvs0 = qvs0/Qs
# Z=np.arange(0,1.5,1e-03)
# Bouyancypoints=[At.bouyancyforcematlab(z, T0, vpar, 0.0138, 28/Qs) for z in Z]
# plt.plot(np.array(Bouyancypoints)*(Vs/(Ts*60)),Z*Ls)
# plt.plot()

dT = dT0  # Desplazamiento en el tiempo
dZ = dZ0  # Desplazamieno en el espacio

Wblock[1, 0] = W10
Wblock[1, 1] = W11

Tblock[1, 0] = T10
Tblock[1, 1] = T11

QVblock[1, 0] = QV10
QVblock[1, 1] = QV11

QRblock[1, 0] = QR10
QRblock[1, 1] = QR11

QNblock[1, 0] = QN10
QNblock[1, 1] = QN11


for i in range(2):
    Wblock[0, 1] = Wblock[1,1] + tau_w * 1 / (1+ tau_w) * (At.getbouyancyforce(z0 + dZ * i, T11, vpar, QVblock[1, 1], QRblock[1, 1]) -
                                                        At.getbouyancyforce(z0 + dZ * (i - 1), T11, vpar, QVblock[1, 1], QRblock[1, 1]))
    Wblock[1, 1] = Wblock[0, 1]
    CFL=dT / dZ
    Tblock[0, 1] = Tblock[1, 1] - CFL *(Wblock[1, 1] * Tblock[1, 1]- Wblock[1, 0]*Tblock[1, 0]) + dT * LCP * (
                        At.auxcdev(T11, QN11, QNblock[1, 1], gamma, QVblock[1, 1], qvs0, QRblock[1, 1],taue, q_star, z0 + i * dZ))
    Tblock[1, 1] = Tblock[0, 1]

    QVblock[0, 1] = QVblock[1, 1] - CFL * (Wblock[1, 1] * QVblock[1, 1]- Wblock[1, 0]*QVblock[1, 0])- dT * (
                        At.auxcdev(T11, QN11, QNblock[1, 1], gamma, QVblock[1, 1], qvs0, QRblock[1, 1],taue, q_star, z0 + i * dZ))
    QVblock[1, 1] = QVblock[0, 1]

    QRblock[0, 1] = QRblock[1, 1] - CFL * (QRblock[1, 1] * (Wblock[1, 1]-At.get_terminalvelocity(vt0, QRblock[1, 1], q_star)) -
                                           QRblock[1, 0] * (Wblock[1, 0]-At.get_terminalvelocity(vt0, QRblock[1, 0], q_star))) - dT * (
                            At.auxcdev(T11, QN11, QNblock[1, 1], gamma, QVblock[1, 1], qvs0, QRblock[1, 1],taue, q_star, z0 + i * dZ))
    QRblock[1, 1] = QRblock[0, 1]

    QNblock[0, 1] = np.exp(-QNblock[1, 1])
    QNblock[1, 1] = QNblock[0, 1]