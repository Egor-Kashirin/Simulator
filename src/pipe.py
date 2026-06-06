import math
import numpy as np


class Pipe:
    def __init__(self, L, D, roughness, H=0.0):
        self.L = L
        self.D = D
        self.roughness = roughness
        self.H = H
        self.g = 9.81

    def _calc_lambda(self, Re, rel_rough):
        if Re < 1e-6:
            return 0.0
        if Re < 2300.0:
            return 64.0 / Re

        lam = 0.02
        for _ in range(50):
            term = rel_rough / 3.7 + 2.51 / (Re * math.sqrt(lam))
            if term <= 0: break
            lam_new = 1.0 / (-2.0 * math.log10(term)) ** 2
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def calculate_dP(self, P, q, interp):
        rho_kgm3 = (P * 101325 * 0.01604) / (interp['Z'].predict(P) * 8.314 * 310.0)
        Bg = interp['Bg'].predict(P)
        mu_cP = interp['visc'].predict(P)
        mu_Pa_s = mu_cP / 1000.0

        v = (4.0 * q * Bg) / (math.pi * self.D ** 2 * 86400.0)

        Re = (rho_kgm3 * v * self.D) / mu_Pa_s if mu_Pa_s > 0 else 1e6

        rel_rough = self.roughness / self.D
        lam = self._calc_lambda(Re, rel_rough)

        delta_P_friction_Pa = lam * (self.L / self.D) * (rho_kgm3 * v ** 2 / 2.0)
        delta_P_hydro_Pa = rho_kgm3 * self.g * self.H
        delta_P_total_Pa = delta_P_friction_Pa + delta_P_hydro_Pa

        return delta_P_total_Pa / 101325.0

    def get_vlp(self, P_man, q_values, interp):
        qs = []
        pbhps = []
        for q in q_values:
            dP = self.calculate_dP(P_man, q, interp)
            P_bhp = P_man + dP
            qs.append(q)
            pbhps.append(P_bhp)
        return qs, pbhps