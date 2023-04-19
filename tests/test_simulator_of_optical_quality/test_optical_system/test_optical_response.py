import numpy as np
import pytest

from simulator_of_optical_quality.optical_system.optical_response import (
    Aperture,
    Optical_psf,
)


def test_initialization():
    ap = Aperture()
    assert isinstance(ap.center[0], int)
    assert isinstance(ap.center[1], int)
    assert ap[ap.center] == 0


def test_copy():
    ap = Aperture(np.zeros((11, 11)), radius=4, obstruction_radius=2)
    bp = ap.copy_with(radius=3)
    assert ap.radius == 4
    assert bp.radius == 3
    assert np.array_equal(bp.array, np.zeros((11, 11)))


def test_disk():
    ap = Aperture.disk(30)
    assert ap[ap.center] == 255
    assert ap[(0, 0)] == 0


def test_from_aperture():
    """
     According to Airy disk Theory, https://en.wikipedia.org/wiki/Airy_disk#Cameras
     we should hit the first minimum at 1.22 * λ * efl / diameter.
     For a perfect aperture without obstruction
    """
    radius = 175  # mm
    wave_length = 560e-9  # m
    focal_length = 2.3  # m
    unit = 1e-3

    ap = Aperture.disk(radius, unit)
    ap.add_padding(4 ** 2)
    psf = Optical_psf.from_aperture(ap, wave_length, focal_length)
    theoretical_first_minimum_ring_radius = 1.22 * wave_length * focal_length / (2 * ap.radius * ap.unit)
    first_0 = psf.get_first_zero()
    assert first_0 == pytest.approx(theoretical_first_minimum_ring_radius, abs=1e-7)


def test_mtf():
    """
    According to William B. Wetherell in its paper
    "THE USE OF OF IMAGE IMAGE QUALITY QUALITY CRITERIA CRITERIA THE USE DESIGNING A DIFFRACTION LIMITED IN DESIGNING
    A DIFFRACTION LIMITED LARGE SPACE TELESCOPE"
    MTF of a perfect round aperture will hit zero at a cut off frequency.
    """
    radius = 15  # mm
    wave_length = 800e-9  # m
    focal_length = 91e-3  # m
    pixel_pitch = 6.5e-6  # m

    ap = Aperture.disk(radius)
    ap.add_padding(4)
    psf = Optical_psf.from_aperture(ap, wave_length, focal_length)
    mtf_sampling, mtf = psf.mtf(unit=pixel_pitch)

    theoretical_cut_off_frequency = (2 * ap.radius * 1e-3) / (wave_length * focal_length) * pixel_pitch

    mtf_zero_index = (mtf < 1e-5).argmax()
    assert mtf_sampling[mtf_zero_index] == pytest.approx(theoretical_cut_off_frequency, abs=0.1)
