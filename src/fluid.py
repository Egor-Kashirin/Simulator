import pandas as pd
from src.interpolator import LinearInterpolator

P_STD = 101325.0
T_STD = 293.15
R = 8.314

class Fluid:
    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        self.M = M
        self.rho_c = rho_c
        self.xa = xa
        self.xy = xy
        self.T = T

        try:
            df = pd.read_csv('interp_data.csv')
            self.mu_interp = LinearInterpolator(df['P'].tolist(), df['mu'].tolist())
        except Exception:
            self.mu_interp = LinearInterpolator([1, 100], [0.011, 0.018])

    def z(self, P_atm: float) -> float:
        P_pa = P_atm * 101325
        return max(0.7, min(1.0, 1.0 - 0.000004 * P_pa))

    def ro(self, P_atm: float) -> float:
        P_pa = P_atm * 101325
        Z = self.z(P_atm)
        return (P_pa * self.M) / (Z * R * self.T)

    def bg(self, P_atm: float) -> float:
        Z = self.z(P_atm)
        return (P_STD * Z * self.T) / (P_atm * 101325 * T_STD)

    def mu(self, P_atm: float) -> float:
        return self.mu_interp.predict(P_atm)