from optical_system import zernike_polynomial as zp
import pytest
import numpy as np


def test_zernike_index():
    assert zp.ZernikeIndex((0, 0)).fringe == 1
    assert zp.ZernikeIndex((1, 1)).fringe == 2
    assert zp.ZernikeIndex((4, -2)).fringe == 13
    assert zp.ZernikeIndex(6).radius == 2
    assert zp.ZernikeIndex(6).azimuth == -2
    assert zp.ZernikeIndex(4).azimuth == 0
    assert zp.ZernikeIndex(4).radius == 2

    with pytest.raises(AssertionError):
        zp.ZernikeIndex((1, 2))


def test_zernikecoeeficient():
    zernike_coef = zp.ZernikeCoefficient(9, 'Focus', np.sqrt(3))
    assert zernike_coef == zp.ZernikeCoefficient(9, 'Focus', np.sqrt(3))
    assert zernike_coef.radius == 4


def test_calculate():
    zernike_cos = zp.ZernikeCoefficients(20)
    _rho = 0.5
    _theta = 0.3

    zernike_co = zernike_cos[1]
    assert zernike_co.calculate_zernike_polynomial(np.asarray(_rho),
                                                   np.asarray(_theta)) == 1

    zernike_co = zernike_cos[4]
    assert zernike_co.radial_polynomial(np.asarray(_rho)) == 2 * _rho ** 2 - 1

    zernike_co = zernike_cos[5]
    manual_calculation = _rho ** 2 * np.cos(2 * _theta) * zernike_co.norm
    assert zernike_co.calculate_zernike_polynomial(
        np.asarray(_rho), np.asarray(_theta)) == manual_calculation


def test_normalization():
    zernike_cos = zp.ZernikeCoefficients(6)
    zernike_co = zernike_cos[5]

    r = np.linspace(0, 1, 1000)
    theta = np.linspace(0, 2 * np.pi, 1000)
    [r2, theta2] = np.meshgrid(r, theta)

    zernike_values = zernike_co.calculate_zernike_polynomial(r2, theta2)
    l2_norm = np.sum(zernike_values ** 2 * r2) / r2.size * 2 * np.pi
    assert l2_norm == pytest.approx(np.pi, abs=1e-2)


def test_init():
    with pytest.raises(zp.ZernikeCoefficientsError):
        zp.ZernikeCoefficients(300)

    zernike_cos = zp.ZernikeCoefficients(3)
    assert zernike_cos[1] == zp.ZernikeCoefficient(1, "Piston", 1)
    assert zernike_cos[1].radius == 0
    assert zernike_cos[(0, 0)] == zp.ZernikeCoefficient(1, "Piston", 1)
