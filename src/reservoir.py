import math

class Reservoir:
    def __init__(self, k, h, re, rw, P_res, beta=0.00852702):
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.P_res = P_res
        self.beta = beta

    def get_productivity_index(self, mu_cP):
        return (self.beta * self.k * self.h) / (mu_cP * math.log(self.re / self.rw))

    def get_ipr_linear(self, P_bhp, C):
        if P_bhp >= self.P_res:
            return 0.0
        return max(0, C * (self.P_res - P_bhp))

    def get_ipr_square(self, P_bhp, C):
        if P_bhp >= self.P_res:
            return 0.0
        return max(0, C * (self.P_res ** 2 - P_bhp ** 2))