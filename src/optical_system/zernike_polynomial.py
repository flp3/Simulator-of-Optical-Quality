import numpy as np
from math import factorial, copysign
import matplotlib.pyplot as plt
import typing
from numpy.typing import NDArray as nd


class ZernikeIndexError(Exception):
    pass


class ZernikeIndex():
    def __init__(self, index: typing.Union[tuple[int, float], int]):
        """
        Zernike polynomials are used to analytically represent the shape and
        local deformations of optical surfaces or wavefronts.
        Many instruments and commercial softwares had been using them such as
        Zemax. However, due to many practical and computational reasons,
        they have slightly different structures and orderings.
        Here we are considering the Fringe representation,
        from the university of Arizona.
        Principle usage: from a double index, i.e. radius and azimuth,
        obtain a unique Fringe index. and reciprocally.

        Args:
            index (tuple or int): if a tuple, then we use the radial and
            azimuthal representation,
            if a unique int, then we use Fringe representation.
        """
        if isinstance(index, tuple):
            assert index[0] >= 0
            assert abs(index[1]) <= index[0]
            self.radius: int = index[0]
            self.azimuth: float = index[1]
            self.fringe: float = self.get_fringe
        if isinstance(index, int):
            self.fringe = index
            self.radius, self.azimuth = self.get_radius_azimuth_index

    @property
    def get_fringe(self) -> float:
        """
        From radial and azimuthal representation to Fringe one
        """
        abs_az = np.abs(self.azimuth)
        sig_az = copysign(1, self.azimuth)
        abs_m = (self.radius + abs_az) / 2

        fringe_index: float = (1 + abs_m) ** 2 - 2 * abs_az + (1 - sig_az) / 2
        return fringe_index

    @property
    def get_radius_azimuth_index(self) -> tuple[int, float]:
        """
        From Fringe representation to radial and azimuthal one.
        """
        index_catalog = {}
        for n in range(10):
            for m in range(-n, n + 1):
                z = ZernikeIndex((n, m))
                f = z.get_fringe
                if f.is_integer():
                    index_catalog[f] = (n, m)

        if self.fringe not in index_catalog:
            raise ZernikeIndexError('Index: %s is not supported' % self.fringe)
        return index_catalog[self.fringe]


class ZernikeCoefficient(ZernikeIndex):
    def __init__(self,
                 index: typing.Union[tuple[int, float], int],
                 name: str,
                 norm: float):
        """
        The mathematical Zernike Polynomial were originally described by Frits
        Zernike in 1934. in:
        Zernike F, "Beugungstheorie des Schneidenver-Fahrens und Seiner
        Verbesserten Form der Phasenkontrastmethode",
        Physica, 1(1934)689-704; doi.org/10.1016/S0031-8914(34)80259-5.
        And were developed to describe the diffracted wavefront in
        phase contrast imaging, in:
        Zernike F, "The Diffraction Theory of Aberrations",
        in Optical Image Evaluation Circular 526,
        ( National Bureau of Standards , Washington, D. C), 1952.
        Moreover, Zernike won the 1953 Nobel Prize in Physics
        for developing Phase Contrast Microscopy.

        Their main application is for optical design and testing,
        describing complex shapes such as freeform surfaces
        and fabrication errors.
        Zernike polynomials have nice mathematical properties, indeed,
        they are orthogonal over the continuous unit circle.
        Meaning They can can represent arbitrarily complex continuous surfaces
        given enough terms.
        Furthermore, their first terms efficiently represent common errors
        (e.g. coma, spherical aberration) seen in optics.

        Args:
            index (tuple or int): if a tuple, then we use the radial and
            azimuthal representation,
            if a unique int, then we use Fringe representation.
            name (string): several Zernike terms (or coefficients) represent
            common optical aberration, such as coma, astigmatism, defocus...
            norm (float): We may want to normalize the Zernike
            coefficients, considering the L2 norm on the unit disk should be
            equal to pi.
            (integral over rho[0, 1] and theta[0, 2*pi] of (Z**2)*rho = pi)
        """
        super().__init__(index)
        self.name: str = name
        self.norm: float = norm
        self.max_rho = 1

    def __repr__(self) -> str:
        return 'ZernikeCoefficient(fringe_index= %s, radius= %s, azimuth= %s, ' \
               'name = %s, norm factor= %s' % (self.fringe,
                                               self.radius,
                                               self.azimuth,
                                               self.name,
                                               self.norm)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ZernikeCoefficient):
            return NotImplemented

        fring_bool = (self.fringe == other.fringe)
        name_bool = (self.name == other.name)
        norm_bool = (self.norm == other.norm)
        return fring_bool & name_bool & norm_bool

    def radial_polynomial(self, rho: nd[np.float64]) -> nd[np.float64]:
        """
        This method calculate the Zernike radial component following the
        ISO 24157 Zernikes standard, on a rectangular grid.

        Args:
            rho: 2d matrix. In order to evaluate the zernike radial
            component on a disk, we consider rho as a sqrt(x**2 + y**2) on a
            rectangular grid and then set everything to 0 above the max_rho,
            which represent the disk radius, usually 1.
        Returns:
            radial_polynomial(2D matrix): The Zernike radial component on a
            rectangular grid.
        """
        if (self.radius - abs(self.azimuth)) % 2 != 0:
            if self.fringe == 1:
                radial_polynomial = np.ones(rho.shape)
                radial_polynomial[rho > self.max_rho] = 0
                return radial_polynomial
            else:
                return np.zeros(rho.shape)

        param1 = (self.radius - abs(self.azimuth)) / 2.0
        param2 = (self.radius + abs(self.azimuth)) / 2.0
        radial_polynomial = np.zeros(rho.shape)

        for k in range(int(param1) + 1):
            rho_sq = rho ** (self.radius - 2 * k)
            numerator = (-1) ** k * factorial(self.radius - k)
            f_k = factorial(k)
            f_p1 = factorial(int(param1 - k))
            f_p2 = factorial(int(param2 - k))
            denominator: float = f_k * f_p1 * f_p2
            radial_polynomial += rho_sq * numerator / denominator

        radial_polynomial[rho > self.max_rho] = 0
        return radial_polynomial

    def calculate_zernike_polynomial(self,
                                     rho: nd[np.float64],
                                     theta: nd[np.float64]) -> nd[np.float64]:
        """
        This method calculate the Zernike coefficient following the ISO 24157
        Zernikes standard, on a rectangular grid.

        Args:
            rho (nd[np.float64]): 2d matrix. In order to evaluate the
            zernike coefficient on a disk,
            we consider rho as a sqrt(x**2 + y**2) on a rectangular quartesian
            grid and than set everything to 0 above the max_rho,
            which represent the disk radius, usually 1.
            theta (nd[np.float64]): 2d matrix. In order to evaluate the
            zernike coefficient on a disk,
            we consider theta as a arctan(y, x) on a rectangular quartesian
            grid

            Note: self.norm (float): We may want to normalize the
            Zernike coefficients, considering the L2 norm on the unit disk
            should be equal to pi.
            (integral over rho[0, 1] and theta[0, 2*pi] of (Z**2)*rho = pi)

        Returns:
            zernike_coefficient (nd[np.float64])
        """
        radial = self.radial_polynomial(rho)
        if self.azimuth >= 0:
            azimuthal = np.cos(self.azimuth * theta)
            zernike_coef: nd[np.float64] = radial * azimuthal * self.norm
        if self.azimuth < 0:
            azimuthal = np.sin(self.azimuth * theta)
            zernike_coef = radial * azimuthal * self.norm
        return zernike_coef

    def cartesian_matrix(self, N: int = 500) -> tuple[nd[np.float64], nd[np.float64]]:
        '''In order to evaluate the zernikes on a rectangular grid, we create
        the quartesian X, Y grid and return the Zernike coefficient'''
        x = np.linspace(-self.max_rho, self.max_rho, N)
        [x2, y2] = np.meshgrid(x, x)
        return x, self.eval_cartesian(x2, y2)

    def eval_cartesian(self, x: nd[np.float64], y: nd[np.float64]) -> nd[np.float64]:
        '''Evaluate the zernikes on predefined cartesian points'''
        rho = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return self.calculate_zernike_polynomial(rho, theta)


class ZernikeCoefficientsError(Exception):
    pass


class ZernikeCoefficients():
    def __init__(self, order: int):
        """Class that handles a set of Zernike Polynomials Coefficients
        using the Fringe representation.

        Args:
            order (int): From the self.max_fringe_order, consider the first
            {order} coefficients.
        """
        if order > self.max_fringe_order:
            raise ZernikeCoefficientsError(
                'order tasked is too high. Maximum supported order is'
                '%s, got %s' % (self.max_fringe_order, order))

        self.order: int = order
        self.coeffs: tuple[ZernikeCoefficient, ...] = self.supported[:order]

    def __getitem__(self,
                    item: typing.Union[tuple[int, float],
                                       int]) -> ZernikeCoefficient | None:
        output = None
        for coef in self.coeffs:
            if isinstance(item, tuple):
                if (coef.radius, coef.azimuth) == item:
                    output = coef
                    break
            else:
                if coef.fringe == item:
                    output = coef
                    break
        return output

    def illustrate_single(self, radius=2, azimuth=2):
        """
        Provide a nice visual too to understand how one Zernike polynomial represent optical aberrations.
        Example of usage:
        a = zp.ZernikeCoefficients(33)
        a.illustrate_single(5, 6)

        Args:
            radius (int, optional):. Defaults to 2.
            azimuth (int, optional):. Defaults to 2.
        """
        if self[(radius, azimuth)] is not None:
            plt.figure()
            _, img = self[(radius, azimuth)].cartesian_matrix()
            plt.title('Fringe: %s, %s' % (self[(radius, azimuth)].fringe, self[(radius, azimuth)].name))
            plt.imshow(img)
            plt.axis('off')
            plt.show()

    def illustrate_panel(self, max_radius=5, max_azimuth=6):
        """
        Provide a nice visual tool to understand how the different Zernike
        coefficients represent optical aberrations.
        Example of usage:
        a = zp.ZernikeCoefficients(33)
        a.illustrate_panel(5, 6)

        Args:
            max_radius (int, optional) Defaults to 5.
            max_azimuth (int, optional) Defaults to 6.
        """
        len_rows = max_radius + 1
        len_columns = 2 * max_azimuth + 1
        fig = plt.figure(figsize=(len_rows, len_columns))
        col = np.linspace(-max_azimuth, max_azimuth, len_columns)

        for radius in range(len_rows):
            for azimuth in range(-radius, radius + 1):
                if self[(radius, azimuth)] is not None:
                    i = np.where(col == azimuth)[0]
                    _, img = self[(radius, azimuth)].cartesian_matrix()
                    idx = int(radius * len_columns + i + 1)
                    ax = fig.add_subplot(len_rows, len_columns, idx)
                    ax.set_title('Fringe: %s, %s' % (self[(radius, azimuth)].fringe, self[(radius, azimuth)].name))
                    ax.imshow(img)
                    ax.axis('off')

        fig.text(0.5, 0.075, 'azimuth index', ha='center', fontsize=14)
        fig.text(0.05, 0.5, 'radial index', va='center', rotation='vertical',
                 fontsize=14)

        plt.show()

    @property
    def max_fringe_order(self):
        return len(self.supported)

    @property
    def supported(self):
        """
        We support the 33 first Fringe indexed Zernike coefficient, as we want
        to normalize them.
        For this we calculate the norm (float) factors, considering the L2
        norm on the unit disk should be equal to pi.
        (integral over rho[0, 1] and theta[0, 2*pi] of (Z**2)*rho = pi)

        Returns:
            tuple: a set of Zernike Coefficients ordered in Fringe index
            representation.
        """
        return (ZernikeCoefficient(1, "No aberration", 1),
                ZernikeCoefficient(2, "Tilt at 0°", 2),
                ZernikeCoefficient(3, "Tilt at 90°", 2),
                ZernikeCoefficient(4, "Defocus", np.sqrt(3)),
                ZernikeCoefficient(5, "Astigmatism at 0°", np.sqrt(6)),
                ZernikeCoefficient(6, "Astigmatism at 45°", np.sqrt(6)),
                ZernikeCoefficient(7, "Coma at 0°", np.sqrt(8)),
                ZernikeCoefficient(8, "Coma at 90°", np.sqrt(8)),
                ZernikeCoefficient(9, "3th spherical aberration", np.sqrt(5)),
                ZernikeCoefficient(10, "Trefoil at 0°", np.sqrt(8)),
                ZernikeCoefficient(11, "Trefoil at 90°", np.sqrt(8)),
                ZernikeCoefficient(12, "5th astigmatism at 0°", np.sqrt(10)),
                ZernikeCoefficient(13, "5th astigmatism at 45°", np.sqrt(10)),
                ZernikeCoefficient(14, "5th order coma at 0°", np.sqrt(12)),
                ZernikeCoefficient(15, "5th order coma at 90°", np.sqrt(12)),
                ZernikeCoefficient(16, "5th spherical aberration", np.sqrt(7)),
                ZernikeCoefficient(17, "Tetrafoil at 0°", np.sqrt(10)),
                ZernikeCoefficient(18, "Tetrafoil at 45°", np.sqrt(10)),
                ZernikeCoefficient(19, "7th order trefoil at 0°", np.sqrt(12)),
                ZernikeCoefficient(20, "7th trefoil at 90°", np.sqrt(12)),
                ZernikeCoefficient(21, "7th astigmatism at 0°", np.sqrt(14)),
                ZernikeCoefficient(22, "7th astigmatism at 45°", np.sqrt(14)),
                ZernikeCoefficient(23, "7th order coma at 0°", 4),
                ZernikeCoefficient(24, "7th order coma at 90°", 4),
                ZernikeCoefficient(25, "7th order spherical aberration", 3),
                ZernikeCoefficient(26, "Pentafoil at 0°", np.sqrt(12)),
                ZernikeCoefficient(27, "Pentafoil at 90°", np.sqrt(12)),
                ZernikeCoefficient(28, "9th tetrafoil at 0°", np.sqrt(14)),
                ZernikeCoefficient(29, "9th tetrafoil at 45°", np.sqrt(14)),
                ZernikeCoefficient(30, "9th order trefoil at 0°", 4),
                ZernikeCoefficient(31, "9th order trefoil at 90°", 4),
                ZernikeCoefficient(32, "9th astigmatism at 0°", np.sqrt(18)),
                ZernikeCoefficient(33, "9th astigmatism at 45°", np.sqrt(18))
                )
