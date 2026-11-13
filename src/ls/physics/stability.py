"""Infinite-slope Mohr-Coulomb factor of safety with suction stress (proposal eq.):

    FS = [c' + (W_t cos^2(beta) - sigma_s) tan(phi')] / [W_t sin(beta) cos(beta)]

sigma_s (kPa) is the net pore-water term on the failure plane:
    sigma_s = u_w - sigma_suction
      u_w            = gamma_w * h_w * cos^2(beta)          (positive pore pressure, perched layer)
      sigma_suction  = gamma_w * psi_f * Se * (1 - z_w/H)   (suction stress, Lu & Likos form;
                                                             decays as the wetting front reaches H)
"""
import numpy as np


def factor_of_safety(cohesion_kpa, phi_deg, slope_rad, bulk_density_gcm3, porosity,
                     theta_i, z_w, h_w, psi_f_m, soil_depth_m=2.0, gamma_w=9.81):
    """All per-cell inputs are [N] (broadcast to [N, T]); theta_i, z_w, h_w are [N, T]."""
    c = np.asarray(cohesion_kpa, dtype=np.float64)[:, None]
    tanphi = np.tan(np.radians(np.asarray(phi_deg, dtype=np.float64)))[:, None]
    b = np.asarray(slope_rad, dtype=np.float64)[:, None]
    n = np.asarray(porosity, dtype=np.float64)[:, None]
    psi = np.asarray(psi_f_m, dtype=np.float64)[:, None]
    gamma_d = np.asarray(bulk_density_gcm3, dtype=np.float64)[:, None] * 9.81    # kN/m3
    H = float(soil_depth_m)

    th = np.minimum(np.asarray(theta_i, dtype=np.float64), n)
    zw = np.clip(np.asarray(z_w, dtype=np.float64), 0.0, H)
    hw = np.clip(np.asarray(h_w, dtype=np.float64), 0.0, H)

    # mean water content of the column: saturated above the front, antecedent below
    theta_col = (zw * n + (H - zw) * th) / H
    theta_col = np.maximum(theta_col, (hw * n + (H - hw) * theta_col) / H)
    W = (gamma_d + gamma_w * theta_col) * H                                        # kPa
    cos2 = np.cos(b) ** 2
    se = th / n
    sigma_suction = gamma_w * psi * se * (1.0 - zw / H)
    u_w = gamma_w * hw * cos2
    sigma_s = u_w - sigma_suction

    num = c + (W * cos2 - sigma_s) * tanphi
    den = W * np.sin(b) * np.cos(b)
    return num / np.maximum(den, 1e-6)


def failure_probability(fs, k=10.0):
    """Modified sigmoid centred at FS = 1: p = 1 / (1 + exp(k (FS - 1)))."""
    z = np.clip(k * (np.asarray(fs) - 1.0), -60, 60)
    return 1.0 / (1.0 + np.exp(z))
