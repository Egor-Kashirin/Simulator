import numpy as np
import pandas as pd
from src.interpolator import LinearInterpolator


class Fluid:
    Pstd = 1
    Tstd = 293.15
    R = 0.00831451

    def __init__(self, rho_c: float, xa: float, xy: float):
        self.rho_c = rho_c
        self.xa = xa / 100
        self.xy = xy / 100
        self.M = 16.04

    def get_Z(self, P: float, T: float) -> float:
        P_mpa = P * 0.101325
        xa = self.xa
        xy = self.xy
        xz = 1 - xa - xy

        z_c = 1 - (0.0741 * self.rho_c - 0.006 - 0.063 * xa - 0.0575 * xy) ** 2
        M3 = (24.05525 * z_c * self.rho_c - 28.0135 * xa - 44.01 * xy) / xz
        H = 128.64 + 47.479 * M3

        B1 = (-0.425468 + 2.865e-3 * T - 4.62073e-6 * T ** 2 +
              (8.77118e-4 - 5.56281e-6 * T + 8.8151e-9 * T ** 2) * H +
              (-8.24747e-7 + 4.31436e-9 * T - 6.08319e-12 * T ** 2) * H ** 2)
        B2 = -0.1446 + 7.4091e-4 * T - 9.1195e-7 * T ** 2
        B23 = -0.339693 + 1.61176e-3 * T - 2.04429e-6 * T ** 2
        B3 = -0.86834 + 4.0376e-3 * T - 5.1657e-6 * T ** 2

        C1 = (-0.302488 + 1.95861e-3 * T - 3.16302e-6 * T ** 2 +
              (6.46422e-4 - 4.22876e-6 * T + 6.88157e-9 * T ** 2) * H +
              (-3.32805e-7 + 2.2316e-9 * T - 3.67713e-12 * T ** 2) * H ** 2)
        C2 = 7.8498e-3 - 3.9895e-5 * T + 6.1187e-8 * T ** 2
        C3 = 2.0513e-3 + 3.4888e-5 * T - 8.3703e-8 * T ** 2
        C223 = 5.52066e-3 - 1.68609e-5 * T + 1.57169e-8 * T ** 2
        C233 = 3.58783e-3 + 8.06674e-6 * T - 3.25798e-8 * T ** 2

        B_star = 0.72 + 1.875e-5 * (320 - T) ** 2
        C_star = 0.92 + 0.0013 * (T - 270)

        Bm = (xz ** 2 * B1 + xz * xa * B_star * (B1 + B2) -
              1.73 * xz * xy * (B1 * B3) ** 0.5 + xa ** 2 * B2 +
              2 * xa * xy * B23 + xy ** 2 * B3)
        Cm = (xz ** 3 * C1 +
              3 * xz ** 2 * xa * C_star * (C1 ** 2 * C2) ** (1 / 3) +
              2.76 * xz ** 2 * xy * (C1 ** 2 * C3) ** (1 / 3) +
              3 * xz * xa ** 2 * C_star * (C1 * C2 ** 2) ** (1 / 3) +
              6.6 * xz * xa * xy * (C1 * C2 * C3) ** (1 / 3) +
              2.76 * xz * xy ** 2 * (C1 * C3 ** 2) ** (1 / 3) +
              xa ** 3 * C2 +
              3 * xa ** 2 * xy * C223 +
              3 * xa * xy ** 2 * C233 +
              xy ** 3 * C3)

        b = 1000 * P_mpa / (2.7715 * T)
        C0 = b ** 2 * Cm
        B0 = b * Bm
        A1 = 1 + B0
        A0 = 1 + 1.5 * (B0 + C0)
        A2 = (A0 - (A0 ** 2 - A1 ** 3) ** 0.5) ** (1 / 3)
        z = (1 + A2 + A1 / A2) / 3
        return z

    def get_Bg(self, P: float, T: float) -> float:
        Z = self.get_Z(P, T)
        Bg = (self.Pstd * Z * T) / (P * self.Tstd)
        return Bg


def load_pvt_data(csv_path: str, rho_c: float, xa: float, xy: float, T: float):
    df = pd.read_csv(csv_path, sep=';')

    fluid = Fluid(rho_c=rho_c, xa=xa, xy=xy)

    col_pressure = [c for c in df.columns if 'pressure' in c.lower()][0]
    col_viscosity = [c for c in df.columns if 'viscosity' in c.lower() or 'mu' in c.lower()][0]

    if 'Z' not in df.columns and 'z' not in df.columns:
        df['Z'] = df[col_pressure].apply(lambda P: fluid.get_Z(P, T))
    if 'Bg' not in df.columns and 'bg' not in df.columns:
        df['Bg'] = df[col_pressure].apply(lambda P: fluid.get_Bg(P, T))

    col_Z = [c for c in df.columns if c.strip().lower() == 'z'][0]
    col_Bg = [c for c in df.columns if c.strip().lower() == 'bg'][0]

    x_data = df[col_pressure].values

    interpolators = {
        'visc': LinearInterpolator(x_data, df[col_viscosity].values),
        'Z': LinearInterpolator(x_data, df[col_Z].values),
        'Bg': LinearInterpolator(x_data, df[col_Bg].values)
    }

    return fluid, interpolators