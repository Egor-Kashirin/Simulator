import numpy as np


class Pipe:
    def __init__(self, D, roughness):
        self.D = D
        self.roughness = roughness

    def _get_lambda(self, Re):
        if Re < 2300:
            return 64.0 / Re

        rel_rough = self.roughness / self.D
        lam = 0.02
        for _ in range(50):
            rhs = -2.0 * np.log10(rel_rough / 3.7 + 2.51 / (Re * np.sqrt(lam)))
            lam_new = 1.0 / rhs ** 2
            if abs(lam_new - lam) < 1e-8:
                break
            lam = lam_new
        return lam

    def calculate_dP(self, q, L, rho, mu_Pas, Bg, T=310.0):
        Q = q / (24 * 3600) * Bg
        A = np.pi * self.D ** 2 / 4
        v = Q / A

        Re = rho * v * self.D / mu_Pas if mu_Pas > 0 else 1e6

        lam = self._get_lambda(Re)

        dP = lam * (L / self.D) * (rho * v ** 2 / 2) / 101325
        return dP