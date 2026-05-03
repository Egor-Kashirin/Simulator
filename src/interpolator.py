class LinearInterpolator:
    def __init__(self, xs: list[float], ys: list[float]):
        if len(xs) != len(ys):
            raise ValueError("xs и ys должны быть одинаковой длины")
        if any(xs[i] > xs[i + 1] for i in range(len(xs) - 1)):
            raise ValueError("xs должен быть отсортирован по возрастанию")
        self.xs = xs
        self.ys = ys

    def predict(self, xp: float) -> float:
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(f"Значение {xp} вне диапазона [{self.xs[0]}, {self.xs[-1]}]")

        for i in range(len(self.xs) - 1):
            if self.xs[i] <= xp <= self.xs[i + 1]:
                x0, x1 = self.xs[i], self.xs[i + 1]
                y0, y1 = self.ys[i], self.ys[i + 1]
                return y0 + (y1 - y0) * (xp - x0) / (x1 - x0)

        return self.ys[-1]