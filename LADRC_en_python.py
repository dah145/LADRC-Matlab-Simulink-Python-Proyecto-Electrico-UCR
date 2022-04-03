from functools import partial
import numpy as np
import math


class LADRC:
    """Clase que contiene métodos para la implementación de un LADRC."""

    def __init__(self, ordenProceso, gananciaNominal, anchoBandaControlador, anchoBandaLESO, condicionInicial,
                 imprimirMat=0):
        """Inicializa las ganancias y construye las matrices que describen el MVE del sistema de control según el
        orden requerido."""

        self.u = 0
        self.h = 0.001  # tamaño del paso de integración del método runkut4
        self.nx = int(ordenProceso)
        self.wc = anchoBandaControlador
        self.wo = anchoBandaLESO
        self.bo = gananciaNominal
        self.zo = condicionInicial
        self.Cg, self.Ac, self.Bc, self.Cc, self.zo, self.L, self.K, self.z = self.ConstruirMatrices(imprimirMat)

    def ConstruirMatrices(self, imprimirMat=0):
        """Función que construye las matrices que componen el LADRC.
        Las matrices son calculadas según el orden escogido para el control.
        Ac, Bc y Cc corresponden a las matrices del MVE del sistema de control LADRC
        en su formulación convencional."""

        n = self.nx + 1

        K = np.zeros([1, self.nx])
        for i in range(self.nx):
            K[0, i] = math.pow(self.wc, n - (i + 1)) * (
                    (math.factorial(self.nx)) / (math.factorial((i + 1) - 1) * math.factorial(n - (i + 1))))

        L = np.zeros([n, 1])
        for i in range(n):
            L[i] = math.pow(self.wo, i + 1) * (
                        (math.factorial(n)) / (math.factorial(i + 1) * math.factorial(n - (i + 1))))

        Cg = 1 / self.bo

        Ac = np.vstack((np.hstack((np.zeros([n - 1, 1]), np.identity(n - 1))), np.zeros([1, n])))
        Bc = np.vstack((np.zeros([self.nx - 1, 1]), self.bo, 0))
        Cc = np.hstack(([[1]], np.zeros([1, n - 1])))
        zo = np.vstack(([[self.zo]], np.zeros([n - 1, 1])))
        z = np.zeros([n, 1])

        if imprimirMat == 1:
            print(f"K =\n {K}\n")
            print(f"L =\n {L}\n")
            print(f"Ac =\n {Ac}\n")
            print(f"Bc =\n {Bc}\n")
            print(f"Cc =\n {Cc}\n")
            print(f"z =\n {z}\n")
            print(f"zo =\n {zo}\n")

        return Cg, Ac, Bc, Cc, zo, L, K, z

    def LESO(self, u, y, z):
        """Calcula los estados del observador lineal de estados extendido."""
        return np.matmul(self.Ac, z) + self.Bc * u + self.L * (y - np.matmul(self.Cc, z))

    def Runkut4(self, F, z, h):
        """Calcula el próximo estado del ADRC usando el método de integración Runge-Kutta de 4to orden."""

        k0 = h * F(z)
        k1 = h * F(z + 0.5 * k0)
        k2 = h * F(z + 0.5 * k1)
        k3 = h * F(z + k2)

        return z + (1 / 6) * (k0 + 2 * k1 + 2 * k2 + k3)

    def SalidaControl(self, r, y):
        """Genera la señal de salida del LADRC."""

        leso = partial(self.LESO, self.u, y)
        self.z = self.Runkut4(leso, self.z, self.h)
        u0 = self.K[0, 0] * (r - self.z[0, 0])

        for i in range(self.nx - 1):
            u0 -= self.K[0, i + 1] * self.z[i + 1, 0]

        self.u = (u0 - self.z[self.nx, 0]) * self.Cg

        return self.u
