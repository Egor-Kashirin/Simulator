import math
from src.state import NodeState
from src.fluid import Fluid


class Pipe:
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid,
                 vertical_depth: float = 0.0, name: str = ""):
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.H = vertical_depth
        self.name = name
        self.g = 9.81

    def _calc_lambda(self, Re: float, rel_rough: float) -> float:
        # Ламинарный режим (Re < 2300)
        if Re <= 2300.0:
            return 64.0 / Re if Re > 0 else 0.0

        # Турбулентный режим (Колбрук-Уайт)
        lam = 0.02  # Начальное приближение
        for _ in range(50):
            # Формула: 1/sqrt(lam) = -2 * log10( eps/3.7D + 2.51/(Re*sqrt(lam)) )
            term = rel_rough / 3.7 + 2.51 / (Re * math.sqrt(lam))
            lam_new = 1.0 / (-2.0 * math.log10(term)) ** 2
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def dp(self, P: float, q: float) -> NodeState:
        if q <= 0.0:
            return NodeState(name=self.name, P_in=P, P_out=P, dP=0.0,
                             q_std=0.0, q_res=0.0, v=0.0, rho=self.fluid.ro(P))

        # 1. Свойства газа при давлении P
        rho = self.fluid.ro(P)  # кг/м3
        mu_cP = self.fluid.mu(P)  # сП
        mu_Pa_s = mu_cP / 1000.0  # Па*с
        bg = self.fluid.bg(P)  # Объемный коэф. расширения

        # 2. Скорость потока (м/с)
        # v = (Q_stand * Bg) / (Area * 86400)
        area = math.pi * (self.D ** 2) / 4.0
        v = (q * bg) / (area * 86400.0)

        # 3. Число Рейнольдса
        Re = (rho * v * self.D) / mu_Pa_s if mu_Pa_s > 0 else 0

        # 4. Коэффициент трения
        rel_rough = self.roughness / self.D
        lam = self._calc_lambda(Re, rel_rough)

        # 5. Перепады давления (Па -> атм)
        # Формула Дарси-Вейсбаха: dP = lambda * (L/D) * (rho * v^2 / 2)
        dP_fric_Pa = lam * (self.L / self.D) * (rho * v ** 2 / 2.0)
        dP_hydro_Pa = rho * self.g * self.H
        dP_total_atm = (dP_fric_Pa + dP_hydro_Pa) / 101325.0

        return NodeState(
            name=self.name,
            P_in=P,
            P_out=P - dP_total_atm,
            dP=dP_total_atm,
            q_std=q,
            q_res=q * bg,
            v=v,
            rho=rho
        )

    def get_vlp(self, P_man: float, q_values: list) -> tuple[list, list]:
        """
        Построение кривой VLP: P_bhp(q) при фиксированном давлении на манифольде.
        Формула: P_bhp = P_man + ΔP_hydrostatic + ΔP_friction
        """
        qs = []
        pbhps = []

        for q in q_values:
            # Свойства газа при давлении на устье (P_man)
            rho = self.fluid.ro(P_man)
            bg = self.fluid.bg(P_man)
            mu_Pa_s = self.fluid.mu(P_man) / 1000.0

            # ⬇️ ГИДРОСТАТИКА СЧИТАЕТСЯ ВСЕГДА (даже при q=0)
            dP_hydro = (rho * self.g * self.H) / 101325.0

            if q <= 0.0:
                # При нулевом расходе: трение = 0, но гидростатика есть!
                P_bhp = P_man + dP_hydro
                qs.append(0.0)
                pbhps.append(P_bhp)
                continue

            # Расчет скорости и трения для q > 0
            area = math.pi * (self.D ** 2) / 4.0
            v = (q * bg) / (area * 86400.0)

            Re = (rho * v * self.D) / mu_Pa_s if mu_Pa_s > 0 else 0
            rel_rough = self.roughness / self.D
            lam = self._calc_lambda(Re, rel_rough)

            dP_fric = (lam * (self.L / self.D) * (rho * v ** 2 / 2.0)) / 101325.0

            # P_bhp = устье + гидростатика + трение
            P_bhp = P_man + dP_hydro + dP_fric

            qs.append(q)
            pbhps.append(P_bhp)

        return qs, pbhps