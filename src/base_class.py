import numpy as np
import typing
import numpy.typing as npt


class BaseClassError(Exception):
    pass


class BaseClass:
    def __init__(self,
                 array:  npt.NDArray[np.float64]=np.zeros(0),
                 sampling:npt.NDArray[np.float64]=np.zeros(0),
                 unit: float=1e-3):
        """
        BaseClass creates an object in the spatial or frequency domain,
        with a 2D array, a 2D sampled track and a unit in meters.

        Args:
            array (np.ndarray, optional): It represents the physical value of the considered object, 
            or in spatial or in the frequency domain. Defaults to np.zeros(0).
            sampling (np.ndarray, optional): It represents the sampling track, 
            taking the np.diff(sampling) will give the sampled frequency. Defaults to None.
            unit (float, optional): It is referenced with meters, then we can select mm by setting it to 1e-3. Defaults to 1e-3.
        """        
        self.array: npt.NDArray[np.float64] = array
        self.center: tuple[np.int64, ...] = self.get_center
        self.sampling: npt.NDArray[np.float64] = sampling
        self.unit: float = unit

    @property
    def get_center(self) -> tuple[np.int64, ...]:
        """
        We consider here the maximum of self.array as the center.
        """
        return np.unravel_index(np.argmax(self.array, axis=None), self.array.shape)

    def zoom(self, roi: int=50) -> npt.NDArray[np.float64]:
        return self.array[self.center[0] - roi:self.center[0] + roi, self.center[1] - roi:self.center[1] + roi]

    def line_profile(self, direction: str='horizontal') -> typing.Union["BaseClass", typing.NoReturn]:
        """
        Dealing with 2D Modulation Transfer Functions (MTF) could be cumbersome,
        assuming some symetric properties of the considered payload, we might prefer to study only a 1D MTF vector.

        Returns:
            a 1D BaseClass or BaseClassError
        """        
        if direction is 'horizontal':
            return BaseClass(array=self.array[self.center[0], self.center[1]:],
                             sampling=self.sampling[self.center[1]:])
        elif direction is 'vertical':
            return BaseClass(array=self.array[self.center[0]:, self.center[1]],
                             sampling=self.sampling[self.center[0]:])
        else:
            raise Exception('direction: %s is not yet Implemented' % direction)

    def fourier_transform_2D(self, output_unit: float=1) -> "BaseClass":
        """
        Considering the BaseClass object, its self.array attribute is a physical value.
        Here we calculate its fourier transfrom in 2D, and we keep track of the sample 2D array in this Fourier transfrom domain.

        Args:
            output_unit (float, optional): We may want to have the MTF in a pixel space. 
            Then, we could set output_unit={pixel_size} Defaults to 1.

        Returns:
            BaseClass: as a BaseClass object, it has a physical 2D value and a 2D sample track
        """        
        mtf = np.abs(np.fft.fftshift(np.fft.fftn(self.array)))
        mtf /= np.amax(mtf)
        spatial_step = np.mean(np.diff(self.sampling))

        mtf_sampling = np.fft.fftfreq(mtf.shape[0], d=spatial_step) * output_unit
        mtf_sampling = np.fft.fftshift(mtf_sampling)

        return BaseClass(array=mtf, sampling=mtf_sampling, unit=output_unit)
