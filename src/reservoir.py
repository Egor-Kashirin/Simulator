from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    P: float
    V: float
    T: float

class Reservoir:
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        rho_std = (101325.0 * self.fluid.M) / (1.0 * 8.314 * 293.15)

        P_cur = self.resprops.P
        rho_cur = self.fluid.ro(P_cur)
        Z_cur = self.fluid.z(P_cur)

        delta_P = (Z_cur * rho_std / rho_cur) * (q_total * dt / self.resprops.V)
        return P_cur - delta_P