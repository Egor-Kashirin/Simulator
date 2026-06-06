class LinearInterpolator:
    def __init__(self, x, y):
        if len(x) != len(y):
            raise ValueError("Разная длина списков")
        self.x = list(x)
        self.y = list(y)

    def predict(self, xp):
        if xp == self.x[-1]:
            x0, x1 = self.x[-2], self.x[-1]
            y0, y1 = self.y[-2], self.y[-1]
            k = (y1 - y0) / (x1 - x0)
            return y1 + k * (xp - x1)

        for i in range(len(self.x) - 1):
            if self.x[i] <= xp < self.x[i + 1]:
                x0, x1 = self.x[i], self.x[i + 1]
                y0, y1 = self.y[i], self.y[i + 1]
                k = (y1 - y0) / (x1 - x0)
                return y0 + k * (xp - x0)

        if xp < self.x[0]:
            return self.y[0]
        return self.y[-1]