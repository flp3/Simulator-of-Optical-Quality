import numpy as np


class BaseClassError(Exception):
    pass


class BaseClass:
    def __init__(self, array=np.zeros(0), sampling=None, unit=1e-3):
        self.array = array
        self.center = self.get_center
        self.sampling = sampling
        self.unit = unit

    @property
    def get_center(self):
        return np.unravel_index(np.argmax(self.array, axis=None), self.array.shape)

    def zoom(self, roi=50):
        return self.array[self.center[0] - roi:self.center[0] + roi, self.center[1] - roi:self.center[1] + roi]

    def line_profile(self, direction='horizontal'):
        if direction not in ['horizontal', 'vertical']:
            raise BaseClassError('direction: %s is not yet Implemented' % direction)
        elif direction is 'horizontal':
            return BaseClass(array=self.array[self.center[0], self.center[1]:],
                             sampling=self.sampling[self.center[1]:])
        elif direction is 'vertical':
            return BaseClass(array=self.array[self.center[0]:, self.center[1]],
                             sampling=self.sampling[self.center[0]:])

    def fourier_transform_2D(self, output_unit=1):
        mtf = np.abs(np.fft.fftshift(np.fft.fftn(self.array)))
        mtf /= np.amax(mtf)
        spatial_step = np.mean(np.diff(self.sampling))

        mtf_sampling = np.fft.fftfreq(mtf.shape[0], d=spatial_step) * output_unit
        mtf_sampling = np.fft.fftshift(mtf_sampling)

        return BaseClass(array=mtf, sampling=mtf_sampling, unit=output_unit)


