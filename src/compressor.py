class DCS:
    """
    Дожимная компрессорная станция (ДКС)
    """
    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        """
        Параметры
        ----------
        CR : float
            Степень сжатия (>= 1.0). CR = P_out / P_in.
        P_line : float
            Давление в магистральном газопроводе [атм].
        q_ext : float
            Сторонний приток газа в шлейф [ст.м³/сут].
        """
        self.CR = CR
        self.P_line = P_line
        self.q_ext = q_ext

    def P_in(self) -> float:
        """Давление на входе в ДКС [атм]"""
        if self.CR <= 1.0:
            return self.P_line  # Компрессор выключен
        return self.P_line / self.CR

    def q_total(self, q_wells: float) -> float:
        """Суммарный расход через ДКС [ст.м³/сут]"""
        return q_wells + self.q_ext