import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from tqdm import tqdm

from src.state import NodeState


class FieldSimulator:
    def __init__(self, reservoir, wells, shlyf, dcs):
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs
        self.P_line = dcs.P_line

    def _system_equations(self, x, P_res):
        q1, q2, q3, P_man = x
        q_list = [q1, q2, q3]
        F = []

        # Уравнения для каждой скважины
        for i, well in enumerate(self.wells):
            q_i = q_list[i]

            # Используем ИНДИВИДУАЛЬНУЮ трубу каждой скважины
            pipe_node = well.pipe.dp(P_man, q_i)
            P_bhp = P_man + pipe_node.dP

            q_ipr = well.q(P_res, P_bhp)
            F.append(q_i - q_ipr)

        # Уравнение для шлейфа
        q_total = q1 + q2 + q3 + self.dcs.q_ext
        P_in_dcs = self.dcs.P_in()

        shlyf_node = self.shlyf.dp(P_in_dcs, q_total)
        P_man_calc = P_in_dcs + shlyf_node.dP

        F.append(P_man - P_man_calc)

        return F

    def solve_step(self, P_res: float, q_init: list = None) -> dict:
        if q_init is None:
            q_init = [1000.0, 1000.0, 1000.0]

        P_man_init = self.dcs.P_in() + 2.0
        x0 = q_init + [P_man_init]

        sol = fsolve(self._system_equations, x0, args=(P_res,), full_output=True)
        x = sol[0]
        info = sol[1]

        q1, q2, q3, P_man = x

        q1 = max(0.0, q1)
        q2 = max(0.0, q2)
        q3 = max(0.0, q3)
        P_man = max(0.0, P_man)

        return {
            'q1': q1, 'q2': q2, 'q3': q3,
            'P_man': P_man,
            'converged': np.max(np.abs(info['fvec'])) < 1.0
        }

    def run(self, N_days: int, dt: float = 1.0, P_res_init: float = None) -> pd.DataFrame:
        if P_res_init is None:
            P_res = self.reservoir.resprops.P
        else:
            P_res = P_res_init

        history = []
        q_prev = [1000.0, 1000.0, 1000.0]
        cumulative_prod = 0.0
        depletion_rate = 5e-5

        for day in tqdm(range(0, N_days, int(dt)), desc="Симуляция"):
            res = self.solve_step(P_res, q_init=q_prev)

            q1, q2, q3 = res['q1'], res['q2'], res['q3']
            P_man = res['P_man']
            q_total = q1 + q2 + q3

            q_prev = [q1, q2, q3]

            prod_step = q_total * dt
            cumulative_prod += prod_step

            P_res_new = P_res - (prod_step * depletion_rate)
            P_res = max(1.0, P_res_new)

            history.append({
                'day': day,
                'P_res': P_res,
                'P_man': P_man,
                'q1': q1, 'q2': q2, 'q3': q3,
                'q_total': q_total,
                'q_cumulative': cumulative_prod
            })

        return pd.DataFrame(history)