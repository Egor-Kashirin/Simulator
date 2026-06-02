class ReservoirProperties:
    """
    Контейнер для статических параметров пласта.
    """
    def __init__(self, P: float, T: float, phi: float, h: float, A: float):
        """
        Параметры
        ----------
        P : float
            Начальное пластовое давление [атм].
        T : float
            Температура пласта [К].
        phi : float
            Пористость (доля, 0..1).
        h : float
            Эффективная толщина пласта [м].
        A : float
            Площадь дренажа [м²].
        """
        self.P = P
        self.T = T
        self.phi = phi
        self.h = h
        self.A = A


class Reservoir:
    """
    Модель пласта. Хранит свойства и предоставляет методы для расчёта объёмов.
    """
    def __init__(self, props: ReservoirProperties):
        self.resprops = props

    def pore_volume(self) -> float:
        """Объём порового пространства, м³"""
        return self.resprops.A * self.resprops.h * self.resprops.phi

    def update_pressure(self, P_new: float):
        """Обновление пластового давления (для материального баланса)"""
        self.resprops.P = P_new