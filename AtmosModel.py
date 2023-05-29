import numpy as np


def approxfqv(z, c):
    # Dr.Gerardo Hernandez Dueñas
    a = 18.04
    b = 3.27
    c0 = 0.1
    d = 3.48
    pz = (1 - b * np.log(1 + c0 * z)) ** d
    fqv = (c / pz) * np.exp(
        - a * (1 / ((1 - b * np.log(1 + c0 * z)) * (1 + c0 * z)) - 1))
    return fqv


def getbouyancyforce(z, theta0, par, qv, qr):
    B = par[0]
    qv0 = par[1]
    LCP = par[2]
    epsilon = par[3]
    thetae = par[4]
    g = par[5]
    b1 = theta0 + z * B
    b2 = approxfqv(z, qv0)
    b3 = b1 + LCP * b2
    b4 = b3 + theta0 * (epsilon - LCP / theta0) * b2
    b5 = thetae + theta0 * (epsilon - LCP / theta0) * qv - theta0 * qr
    b = g / theta0 * (b5 - b4)
    return b


def plotbouyancyfunctions(z, theta0, par, qv, qr, k):
    B = par[0]
    qv0 = par[1]
    LCP = par[2]
    epsilon = par[3]
    thetae = par[4]
    b1 = theta0 + z * B
    b2 = approxfqv(z, qv0)
    b3 = b1 + LCP * b2
    b4 = b3 + theta0 * (epsilon - LCP / theta0) * b2
    b5 = thetae + theta0 * (epsilon - LCP / theta0) * qv - theta0 * qr
    B = np.zeros(5)
    B[0] = b1
    B[1] = b2
    B[2] = b3
    B[3] = b4
    B[4] = b5
    return B[k]


def bouyancyforcematlab(z, theta0, par, qt, qvs0):
    B = par[0]
    qv0 = par[1]
    LCP = par[2]
    epsilon = par[3]
    thetae = par[4]
    g = par[5]
    b1 = theta0 + z * B
    b2 = approxfqv(z, qv0)
    b3 = approxfqv(z, qvs0)
    b4 = b1 + LCP * b2
    b5 = b4 + theta0 * (epsilon - LCP / theta0) * b2
    b6 = np.min([b3, qt])
    b7 = np.max([qt - b3, 0])
    b8 = thetae + theta0 * (epsilon - LCP / theta0) * b6 - theta0 * b7
    b = g / theta0 * (b8 - b5)
    return b


def get_scale_time_condesation(tau_0, qn0, qn, gamma):
    tau_c = tau_0 * np.exp(((qn-qn0) / gamma) ** 2)
    return tau_c


def get_condensation(tau_0, qn0, qn, gamma, qv, qvs0, z):
    cd = 1 / get_scale_time_condesation(tau_0, qn0, qn, gamma) * np.max([qv - approxfqv(z, qvs0), 0])
    return cd


def get_evaporation(qr, tau_e, q_star, qv, qvs0, z):
    ev = qr / (tau_e * q_star) * np.max([approxfqv(z, qvs0) - qv, 0])
    return ev


def get_terminalvelocity(vt0, qr, q_star):
    vt = vt0 * qr / q_star
    return vt


def get_aerosolvelocity(vtnd, vt0, qr, q_star):
    vtn = vtnd + np.min([qr / q_star, 1]) * np.max([get_terminalvelocity(vt0, qr, q_star), 0])
    return vtn


def auxcdev(cte, qv, qvs0, qr, tau_e, q_star, z):
    cdev = get_condensation(cte, qv, qvs0, z) - get_evaporation(qr, tau_e, q_star, qv, qvs0, z)
    return cdev
