import numpy as np
from scipy.optimize import fsolve


class FieldSimulator:
    def __init__(self, reservoir, wells, shlyf, P_line):
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.P_line = P_line

    def _residuals(self, x, P_res, interp):
        q1, q2, q3, P_man = x

        # 1. Защита области определения: давления и дебиты должны быть положительными
        # Это критично для log() в формуле трения и sqrt() в IPR
        qs = [max(1.0, q) for q in [q1, q2, q3]]
        P_man_safe = max(self.P_line + 0.1, P_man)

        res = [0.0] * 4

        for i, well in enumerate(self.wells):
            # Внутренняя итерация для нахождения P_bhp при фиксированных q и P_man
            P_bhp = P_man_safe + 10.0
            for _ in range(5):
                try:
                    # Считаем потери в трубе
                    dP_tube = well.pipe.calculate_dP(P_man_safe, qs[i], interp)
                    P_bhp_new = P_man_safe + dP_tube
                    if abs(P_bhp_new - P_bhp) < 0.01: break
                    P_bhp = P_bhp_new
                except:
                    break

            try:
                # Дебит из пласта (IPR)
                q_ipr = well.get_q(P_bhp)
            except:
                q_ipr = 0.0

            # Невязка: разница между заданным q и тем, что дает пласт
            res[i] = qs[i] - q_ipr

        # Баланс шлейфа
        q_total_sys = sum(qs)
        try:
            dP_shlyf = self.shlyf.calculate_dP(P_man_safe, q_total_sys, interp)
        except:
            dP_shlyf = 10.0  # Штрафное значение, если ошибка

        # Невязка давления: P_man должно равняться P_line + потери
        res[3] = P_man_safe - (self.P_line + dP_shlyf)

        return res

    def solve(self, P_res, interp):
        # Улучшенное начальное приближение
        # Если начать слишком далеко, fsolve может уйти в недопустимую область (отрицательные числа)
        x0 = [100000.0, 100000.0, 100000.0, self.P_line + 10.0]

        try:
            # full_output=True позволяет получить информацию о сходимости, но мы ее подавим
            sol = fsolve(self._residuals, x0, args=(P_res, interp), full_output=True)
            solution = sol[0]

            q1, q2, q3, P_man = solution

            # Фильтрация результатов: отрицательные дебиты невозможны
            qs = [max(0.0, q) for q in [q1, q2, q3]]
            P_man = max(self.P_line, P_man)

            return P_man, qs
        except Exception as e:
            # В случае полного краха возвращаем безопасные значения
            return self.P_line + 5.0, [0.0, 0.0, 0.0]

    def get_states(self, P_res, interp):
        P_man, qs = self.solve(P_res, interp)
        states = {}
        for i, q in enumerate(qs):
            well = self.wells[i]
            try:
                dP_well = well.pipe.calculate_dP(P_man, q, interp)
                P_bhp = P_man + dP_well
            except:
                P_bhp = P_man
            states[f'well_{i + 1}'] = {'q': q, 'P_bhp': P_bhp, 'P_man': P_man}

        q_total = sum(qs)
        try:
            dP_shlyf = self.shlyf.calculate_dP(P_man, q_total, interp)
        except:
            dP_shlyf = 0.0

        states['shlyf'] = {'P_in': P_man, 'P_out': self.P_line, 'dP': dP_shlyf}

        return states, P_man, q_total