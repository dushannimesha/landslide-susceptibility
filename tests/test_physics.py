import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ls.physics import factor_of_safety, failure_probability, green_ampt_series, rawls_psi_f  # noqa: E402


def test_dry_cohesionless_matches_tanphi_over_tanbeta():
    beta = np.radians([30.0])
    phi = 35.0
    fs = factor_of_safety([0.0], [phi], beta, [1.6], [0.4], theta_i=np.zeros((1, 1)),
                          z_w=np.zeros((1, 1)), h_w=np.zeros((1, 1)), psi_f_m=[0.0], soil_depth_m=2.0)
    assert np.isclose(fs[0, 0], np.tan(np.radians(phi)) / np.tan(beta[0]), rtol=1e-6)


def test_hand_calculation_with_cohesion():
    # W = gamma_d * H = 1.6*9.81*2 = 31.392 kPa, beta 30, phi 30, c 5
    b = np.radians(30.0)
    W = 1.6 * 9.81 * 2
    expected = (5 + W * np.cos(b) ** 2 * np.tan(np.radians(30))) / (W * np.sin(b) * np.cos(b))
    fs = factor_of_safety([5.0], [30.0], [b], [1.6], [0.4], np.zeros((1, 1)), np.zeros((1, 1)),
                          np.zeros((1, 1)), [0.0], soil_depth_m=2.0)
    assert np.isclose(fs[0, 0], expected, rtol=1e-6)


def test_fs_decreases_with_rain():
    N, T = 1, 30
    rain = np.full((N, T), 60.0)
    theta = np.full((N, T), 0.30)
    slope = np.radians([35.0])
    psi = rawls_psi_f(np.array([40.0]), np.array([30.0]), np.array([0.45]))
    ga = green_ampt_series(rain, theta, np.array([5e-6]), np.array([0.45]), psi, slope,
                           soil_depth_m=2.0, tau_days=20)
    fs = factor_of_safety([6.0], [30.0], slope, [1.5], [0.45], theta, ga["z_w"], ga["h_w"], psi, 2.0)
    assert fs[0, -1] < fs[0, 0]
    assert np.all(np.diff(ga["z_w"][0]) >= -1e-9)          # front advances monotonically under constant rain
    assert ga["infil_mm"].max() <= 60.0 + 1e-6            # cannot infiltrate more than rain


def test_infiltration_capped_by_capacity():
    rain = np.full((1, 5), 500.0)
    ga = green_ampt_series(rain, np.full((1, 5), 0.2), np.array([1e-7]), np.array([0.4]),
                           np.array([0.1]), np.radians([20.0]), tau_days=5)
    assert ga["infil_mm"][0, -1] < 500.0


def test_sigmoid_centre():
    assert np.isclose(failure_probability(1.0, 10), 0.5)
    assert failure_probability(0.8, 10) > 0.5 > failure_probability(1.2, 10)
