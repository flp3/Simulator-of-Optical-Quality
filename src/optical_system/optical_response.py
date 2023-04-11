import numpy as np
from base_class import BaseClass
import cv2


class ApertureError(Exception):
    pass


class Aperture:
    def __init__(self,
                 array=np.zeros((1, 1),dtype=np.complex128),
                 radius=0,
                 obstruction_radius=0,
                 sampling=None,
                 unit=1e-3):
        self.array = array
        self.radius = radius
        self.center = self.get_center
        self.obstruction_radius = obstruction_radius
        self.spatial_sampling = sampling
        self.unit = unit

    def __getitem__(self, item):
        return self.array[item]

    def __repr__(self):
        return 'Aperture(array=%s, radius=%s, center=%s, obstruction_radius=%s)' % (self.array,
                                                                                    self.radius,
                                                                                    self.center,
                                                                                    self.obstruction_radius)

    def copy_with(self, **kwargs):
        init_args = {
            'array': self.array,
            'radius': self.radius,
            'obstruction_radius': self.obstruction_radius,
            'sampling': self.spatial_sampling
            }
        init_args.update(kwargs)
        return Aperture(**init_args)

    def __copy__(self):
        return self.copy_with()

    @property
    def get_center(self):
        if self.array.shape[0] % 2 != 1 or self.array.shape[1] % 2 != 1:
            raise ApertureError('prefered odd size')
        return (self.array.shape[0] // 2, self.array.shape[1] // 2)

    @classmethod
    def disk(cls, radius, unit=1e-3):
        shape = (radius * 2 + 1, radius * 2 + 1)
        array = np.zeros(shape)
        spatial_sampling = np.linspace(-radius, radius, shape[0]) * unit
        output = cls(array, radius, sampling=spatial_sampling, unit=unit)
        image = cv2.circle(array, output.get_center, radius, 255, -1)
        output.array = image
        return output

    @classmethod
    def disk_with_obstruction(cls, radius, obstruction_radius):
        disk_aperture = cls.disk(radius)
        new_array = disk_aperture.array.copy()
        new_array = cv2.circle(new_array, disk_aperture.get_center, obstruction_radius, 0, -1)
        return disk_aperture.copy_with(array=new_array)

    @classmethod
    def disk_obstruction_spider(cls, radius, obstruction_radius):
        disk_obstruction = Aperture.disk_with_obstruction(radius, obstruction_radius)
        new_array = disk_obstruction.array.copy()
        new_array[disk_obstruction.center[0], :] = 0
        new_array[:, disk_obstruction.center[1]] = 0
        return disk_obstruction.copy_with(array=new_array)

    def add_phase(self, zernike_phase, factor=1):
        self.array = self.array.astype(np.complex128) * np.exp(-1j * 2 * np.pi * zernike_phase * factor)

    def add_padding(self, padding=16):
        npad_y = self.array.shape[0] * padding // 4
        npad_x = self.array.shape[1] * padding // 4
        self.array = np.pad(self.array, ((npad_y, npad_y), (npad_x, npad_x)), 'constant')
        self.center = self.get_center


class Optical_psf(BaseClass):
    def __init__(self, array=np.zeros(0), aperture=Aperture(), wave=450e-9, foc=2.1, sampling=None, unit=1e-3):
        super().__init__(array, sampling, unit)
        self.ap = aperture
        self.wave_length = wave
        self.focal_length = foc

    @classmethod
    def from_aperture(cls, aperture, wave_length=450e-9, focal_length=2.1):
        Fap = np.abs(np.fft.fftshift(np.fft.fftn(aperture.array))) ** 2
        spatial_step = np.mean(np.diff(aperture.spatial_sampling))
        psf_sampling = np.fft.fftfreq(Fap.shape[0], d=spatial_step) * wave_length * focal_length
        psf_sampling = np.fft.fftshift(psf_sampling)
        return Optical_psf(array=Fap, aperture=aperture, sampling=psf_sampling, unit=aperture.unit)

    def mtf(self, unit=4.5e-6):
        mtf = self.fourier_transform_2D(output_unit=unit).line_profile()
        return mtf.sampling, mtf.array
