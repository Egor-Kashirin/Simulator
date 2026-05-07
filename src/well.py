import math
from src.fluid import Fluid
from src.pipe import Pipe

BETA = 0.00852702


class Well:
    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float, pipe: Pipe = None):
        self.fluid = fluid
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.pipe = pipe

    def q(self, P_res: float, P_bhp: float) -> float:

        if P_res <= P_bhp:
            return 0.0
        mu = self.fluid.mu(P_res)
        if mu <= 0 or self.rw <= 0 or self.re <= self.rw:
            return 0.0
        ln_ratio = math.log(self.re / self.rw)
        C = (BETA * self.k * self.h) / (mu * ln_ratio)
        return C * (P_res - P_bhp)

    def get_ipr_curve(self, P_res: float, n_points: int = 20) -> list:
        """Генерация точек кривой IPR: (P_bhp, q)"""
        curve = []
        for i in range(n_points + 1):
            P_bhp = P_res * (1.0 - i / n_points)
            curve.append((P_bhp, self.q(P_res, P_bhp)))
        return curve

    def get_vlp_curve(self, P_man: float, max_q: float = 1500.0, n_points: int = 20) -> list:
        """Генерация точек кривой VLP: (P_bhp, q) через трубу"""
        if self.pipe is None:
            raise ValueError("Pipe object is required for VLP")

        curve = []
        for i in range(n_points + 1):
            q_val = max_q * (i / n_points)
            # Считаем перепад в трубе для данного дебита
            node = self.pipe.dp(P_man, q_val)
            # P_bhp = P_man + потери
            curve.append((P_man + node.dP, q_val))
        return curve