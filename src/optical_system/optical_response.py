import numpy as np
from numpy.typing import NDArray
from base_class import BaseClass
import cv2
import typing


class ApertureError(Exception):
    pass


class Aperture:
    def __init__(self,
                 array: NDArray[np.complex128] = np.zeros((1, 1), dtype=np.complex128),
                 radius: int = 0,
                 obstruction_radius: int = 0,
                 sampling: NDArray[np.float64] = np.zeros(0),
                 unit: float = 1e-3):
        """
        Aperture class simulates different sort of optical aperture in
        the spatial domain: Lenses and Cassegrain telescopes.
        As this Class will be used to simulate optical acquisition system,
        we can add a phase which represents the manufacturing unperfections.

        Args:
            array (np.ndarray, optional): It represents the optical system
            aperture shape.
            Defaults to np.zeros((1, 1),dtype=np.complex128).
            radius (int, optional): So far we consider circular aperture shape,
            with a radius. Defaults to 0.
            obstruction_radius (int, optional): If we want to add an
            obstruction due to a secondary mirror in a Telescope. Defaults to 0
            sampling (np.ndarray, optional):  It represents the sampling track,
            taking the np.diff(sampling) will give the sampled frequency.
            Defaults to None.
            unit (float, optional):  It is referenced with meters,
            then we can select mm by setting it to 1e-3. Defaults to 1e-3.
        """

        self.array: NDArray[np.complex128] = array
        self.radius: int = radius
        self.center: tuple[int, int] = self.get_center
        self.obstruction_radius: int = obstruction_radius
        self.spatial_sampling: NDArray[np.float64] = sampling
        self.unit: float = unit

    def __getitem__(self, item: int) -> typing.Any:
        return self.array[item]

    def __repr__(self) -> str:
        return 'Aperture(array=%s,radius=%s,center=%s, obstruction_radius=%s)'\
            % (self.array, self.radius, self.center, self.obstruction_radius)

    def copy_with(self, **kwargs: typing.Any) -> "Aperture":
        init_args = {
            'array': self.array,
            'radius': self.radius,
            'obstruction_radius': self.obstruction_radius,
            'sampling': self.spatial_sampling
            }
        init_args.update(kwargs)
        return Aperture(**init_args)

    def __copy__(self) -> "Aperture":
        return self.copy_with()

    @property
    def get_center(self) -> tuple[int, int]:
        """
        We consider the self.array spatial center as the center here.
        Its principal use is for zerro padding the self.array
        """
        if self.array.shape[0] % 2 != 1 or self.array.shape[1] % 2 != 1:
            raise ApertureError('prefered odd size')

        return (self.array.shape[0] // 2, self.array.shape[1] // 2)

    @classmethod
    def disk(cls, radius: int, unit: float = 1e-3) -> "Aperture":
        """
        This method creates an Aperture Class where the optical system
        aperture is a simple hole letting the light in through a disk.
        Principal usage: to simulate a Lens.

        Args:
            radius (int): The aperture is circular, a radius is its unique
            parameter.
            unit (float, optional): depending if we want to consider mm,
            in this case we can set 1e-3, or in pixel space,
            in this case {pixel_size}. Defaults to 1e-3.

        Returns:
            Aperture
        """
        shape = (radius * 2 + 1, radius * 2 + 1)
        array = np.zeros(shape)
        spatial_sampling = np.linspace(-radius, radius, shape[0]) * unit
        output = cls(array.astype(np.complex128), radius, sampling=spatial_sampling, unit=unit)
        image = cv2.circle(array, output.get_center, radius, 255, -1)
        output.array = image.astype(np.complex128)
        return output

    @classmethod
    def disk_with_obstruction(cls,
                              radius: int,
                              obstruction_radius: int) -> "Aperture":
        """
        This method creates an Aperture Class where the optical system
        aperture is a donut.
        Args:
            radius (int): The aperture is a donut,
            a radius is its first parameter.
            obstruction_radius (int): The aperture is a donut,
            an obstruction_radius is its second parameter.

        Returns:
            Aperture
        """
        disk_aperture = cls.disk(radius)
        new_array = disk_aperture.array.copy()
        new_array = cv2.circle(new_array, disk_aperture.get_center, obstruction_radius, 0, -1)
        return disk_aperture.copy_with(array=new_array)

    @classmethod
    def disk_obstruction_spider(cls,
                                radius: int,
                                obstruction_radius: int) -> "Aperture":
        """
        This method creates an Aperture Class where the optical system
        aperture is a donut, mounted on a cross.
        Principal usage: to simulate a Cassegrain Telescope with a secondary
        mirror attached with 4 mechanical legs (called the spider).

        Args:
        Note: Here we assume the spider is 1 unit thick,
        so we did not add an extra parameter to modelize it.
            radius (int): The aperture is a donut, a radius is its first
            parameter.
            obstruction_radius (int): The aperture is a donut,
            an obstruction_radius is its second parameter.

        Returns:
            Aperture
        """
        disk_obstruction = cls.disk_with_obstruction(radius, obstruction_radius)
        new_array = disk_obstruction.array.copy()
        new_array[disk_obstruction.center[0], :] = 0
        new_array[:, disk_obstruction.center[1]] = 0
        return disk_obstruction.copy_with(array=new_array)

    def add_phase(self,
                  zernike_phase: NDArray[np.float64],
                  factor: float = 1) -> None:
        """
        To simulate manufacturing imperfections, we can add a phase to our Optical system aperture.

        Args:
            zernike_phase (NDArray[np.float64]): To modelize the optical
            system aberration, we consider the Zernike decomposition.
            factor (float, optional) Defaults to 1.
        """
        self.array = self.array.astype(np.complex128) * np.exp(-1j * 2 * np.pi * zernike_phase * factor)

    def add_padding(self, padding: int = 16) -> None:
        """
        Principal use: zerro padding is used to have more samples when
        applying a Fourier transfrom.
        Careful, it does not improve the resolution of the obtained Fourier
        transfromed signal.

        Args:
            padding (int, optional): by which factor do we want to extend the
            original signal with 0s. Defaults to 16.
        """
        npad_y = self.array.shape[0] * padding // 4
        npad_x = self.array.shape[1] * padding // 4
        self.array = np.pad(self.array, ((npad_y, npad_y), (npad_x, npad_x)), 'constant')
        self.center = self.get_center


class Optical_psf(BaseClass):
    def __init__(self,
                 array: NDArray[np.float64] = np.zeros(0, dtype=np.float64),
                 aperture: Aperture = Aperture(),
                 wave: float = 450e-9,
                 foc: float = 2.3,
                 sampling: NDArray[np.float64] = np.zeros(0),
                 unit: float = 1e-3):
        """
        Based on a Given optical system aperture, we can obtain the Point
        Spread Function, and therefore its Modulation Transfer Function.
        Args:
            array (NDArray[np.float64], optional): It represent the Point
            spread function in 2D of an optical system.
            Defaults to np.zeros(0, dtype=np.float64).
            aperture (Aperture, optional): It is the simulated optical system
            aperture. Defaults to Aperture().
            wave (float, optional): As light diffraction depends on th
            wavelength, this parameter could improve the accuracy of our
            simulations. Defaults to 450e-9.
            foc (float, optional): the focal length of the optical system,
            in meters. Defaults to 2.3.
            sampling (np.ndarray, optional):  It represents the sampling track,
            taking the np.diff(sampling) will give the sampled frequency.
            Defaults to None.
            unit (float, optional):  It is referenced with meters,
            then we can select mm by setting it to 1e-3. Defaults to 1e-3.
        """
        super().__init__(array, sampling, unit)
        self.ap: Aperture = aperture
        self.wave_length: float = wave
        self.focal_length: float = foc

    @classmethod
    def from_aperture(cls,
                      aperture: Aperture,
                      wave_length: float = 450e-9,
                      focal_length: float = 2.3) -> "Optical_psf":
        """
        From a simulated optical system aperture, if we specify the considered
        wavelength and the focal length of our system,
        this method gives the Point Spread function in Aperture.unit

        Args:
            aperture (Aperture)
            wave_length (float, optional): in meters. Defaults to 450e-9.
            focal_length (float, optional): in meters. Defaults to 2.3.

        Returns:
            Optical_psf
        """
        Fap = np.abs(np.fft.fftshift(np.fft.fftn(aperture.array))) ** 2
        spatial_step = np.mean(np.diff(aperture.spatial_sampling))
        psf_sampling = np.fft.fftfreq(Fap.shape[0], d=spatial_step) * wave_length * focal_length
        psf_sampling = np.fft.fftshift(psf_sampling)
        return Optical_psf(array=Fap.astype(np.float64), aperture=aperture, sampling=psf_sampling, unit=aperture.unit)

    def mtf(self, unit: float = 5.5e-6) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        From the Optical Point Spread Function, applying a Fourier Transform,
        we obtain the Modulation Transfer function of the optical system.

        Args:
            unit (float, optional): If we want to obtain the MTF in the
            pixel space. Defaults to 5.5e-6.

        Returns:
            tuple[NDArray[np.float64], NDArray[np.float64]]: the MTF sample
            track, and its values.
        """
        mtf = self.fourier_transform_2D(output_unit=unit).line_profile()
        return mtf.sampling, mtf.array
