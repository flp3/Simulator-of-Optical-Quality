import numpy as np
import math
import matplotlib.pyplot as plt


class ZernikeIndexError(Exception):
    pass


class ZernikeIndex():
    def __init__(self, index):
        if isinstance(index, tuple):
            assert index[0] >= 0
            assert abs(index[1]) <= index[0]
            self.radius = index[0]
            self.azimuth = index[1]
            self.fringe = self.get_fringe
        if isinstance(index, int):
            self.fringe = index
            self.radius, self.azimuth = self.get_radius_azimuth_index

    @property
    def get_fringe(self):
        return (1 + (self.radius + np.abs(self.azimuth)) / 2) ** 2 -\
               2 * np.abs(self.azimuth) + (1 - math.copysign(1, self.azimuth)) / 2

    @property
    def get_radius_azimuth_index(self):
        index_catalog = {}
        for n in range(10):
            for l in range(-n, n + 1):
                z = ZernikeIndex((n, l))
                f = z.get_fringe
                if f.is_integer():
                    index_catalog[f] = (n, l)

        if self.fringe not in index_catalog:
            raise ZernikeIndexError('this fringe index: %s is not supported yet' % self.fringe)
        return index_catalog[self.fringe]


class ZernikeCoefficient(ZernikeIndex):
    def __init__(self, index, name, normalization):
        super().__init__(index)
        self.name = name
        self.normalization = normalization
        self.max_rho = 1

    def __repr__(self):
        return 'ZernikeCoefficient(fringe_index= %s, radius= %s, azimuth= %s, ' \
               'name = %s, normalization factor= %s' % (self.fringe, self.radius, self.azimuth, self.name,
                                                        self.normalization)

    def __eq__(self, other):
        return (self.fringe == other.fringe) & (self.name == other.name) & (self.normalization == other.normalization)

    def calculate_radial_polynomial(self, rho):
        """
        :param rho: 2d matrix
        :return:
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
            numerator = (-1) ** k * math.factorial(self.radius - k)
            denominator: float = math.factorial(k) * math.factorial(int(param2 - k)) * math.factorial(int(param1 - k))
            radial_polynomial += rho_sq * numerator / denominator

        radial_polynomial[rho > self.max_rho] = 0
        return radial_polynomial

    def calculate_zernike_polynomial(self, rho, theta):
        radial_polynomial = self.calculate_radial_polynomial(rho)
        if self.azimuth >= 0:
            return radial_polynomial * np.cos(self.azimuth * theta) / self.normalization
        if self.azimuth < 0:
            return radial_polynomial * np.sin(self.azimuth * theta) / self.normalization

    def cartesian_matrix(self, N=100):
        '''Evaluate the zernikes on a rectangular grid, and return the X, Y grid and the values'''
        x = np.linspace(-self.max_rho, self.max_rho, N)
        [x2, y2] = np.meshgrid(x, x)
        return x, self.eval_cartesian(x2, y2)

    def eval_cartesian(self, x, y):
        '''Evaluate the zernikes on predefined cartesian points'''
        rho = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return self.calculate_zernike_polynomial(rho, theta)


class ZernikeCoefficientsError(Exception):
    pass


class ZernikeCoefficients():
    """
    Return a set of Zernike Polynomials Coefficient

    The particular set we are using here is the "Fringe" (or "University of Arizona"), and are
    NOT orthonormal. They are normalized to be one at radius=1.
    """

    def __init__(self, order):
        if order > self.max_fringe_order:
            raise ZernikeCoefficientsError('order tasked is too high.'
                                           'Maximum supported order is %s and got %s' % (self.max_fringe_order, order))
        self.order = order
        self.coefficients = self.get_supported_zernike[:self.order]

    def __getitem__(self, item):
        for coef in self.coefficients:
            if isinstance(item, tuple):
                if (coef.radius, coef.azimuth) == item:
                    return coef
            else:
                if coef.fringe == item:
                    return coef

    def illustrate(self, max_radius=5, max_azimuth=6):
        len_rows = max_radius + 1
        len_columns = 2 * max_azimuth + 1
        fig = plt.figure(figsize=(len_rows, len_columns))
        col = np.linspace(-max_azimuth, max_azimuth, len_columns)

        for radius in range(len_rows):
            for azimuth in range(-radius, radius + 1):
                if self[(radius, azimuth)] is None:
                    continue
                i = np.where(col == azimuth)[0]
                x_axis, img = self[(radius, azimuth)].cartesian_matrix()
                ax = fig.add_subplot(len_rows, len_columns, int(radius * len_columns + i + 1))
                ax.set_title('Fringe: %s' % self[(radius, azimuth)].fringe)
                ax.imshow(img)
                ax.axis('off')

        fig.text(0.5, 0.075, 'azimuth index', ha='center', fontsize=14)
        fig.text(0.05, 0.5, 'radial index', va='center', rotation='vertical', fontsize=14)

        plt.show()

    @property
    def max_fringe_order(self):
        return len(self.get_supported_zernike)

    @property
    def get_supported_zernike(self):
        return (ZernikeCoefficient(1, "Piston", 1),
                ZernikeCoefficient(2, "Tilt at 0°", 2),
                ZernikeCoefficient(3, "Tilt at 90°", 2),
                ZernikeCoefficient(4, "Focus", np.sqrt(3)),
                ZernikeCoefficient(5, "Astigmatism at 0°", np.sqrt(6)),
                ZernikeCoefficient(6, "Astigmatism at 45°", np.sqrt(6)),
                ZernikeCoefficient(7, "Coma at 0°", np.sqrt(8)),
                ZernikeCoefficient(8, "Coma at 90°", np.sqrt(8)),
                ZernikeCoefficient(9, "3th order spherical aberration", np.sqrt(5)),
                ZernikeCoefficient(10, "Trefoil at 0°", np.sqrt(8)),
                ZernikeCoefficient(11, "Trefoil at 90°", np.sqrt(8)),
                ZernikeCoefficient(12, "5th order astigmatism at 0°", np.sqrt(10)),
                ZernikeCoefficient(13, "5th order astigmatism at 45°", np.sqrt(10)),
                ZernikeCoefficient(14, "5th order coma at 0°", np.sqrt(12)),
                ZernikeCoefficient(15, "5th order coma at 90°", np.sqrt(12)),
                ZernikeCoefficient(16, "5th order spherical aberration", np.sqrt(7)),
                ZernikeCoefficient(17, "Tetrafoil at 0°", np.sqrt(10)),
                ZernikeCoefficient(18, "Tetrafoil at 45°", np.sqrt(10)),
                ZernikeCoefficient(19, "7th order trefoil at 0°", np.sqrt(12)),
                ZernikeCoefficient(20, "7th order trefoil at 90°", np.sqrt(12)),
                ZernikeCoefficient(21, "7th order astigmatism at 0°", np.sqrt(14)),
                ZernikeCoefficient(22, "7th order astigmatism at 45°", np.sqrt(14)),
                ZernikeCoefficient(23, "7th order coma at 0°", 4),
                ZernikeCoefficient(24, "7th order coma at 90°", 4),
                ZernikeCoefficient(25, "7th order spherical aberration", 3),
                ZernikeCoefficient(26, "Pentafoil at 0°", np.sqrt(12)),
                ZernikeCoefficient(27, "Pentafoil at 90°", np.sqrt(12)),
                ZernikeCoefficient(28, "9th order tetrafoil at 0°", np.sqrt(14)),
                ZernikeCoefficient(29, "9th order tetrafoil at 45°", np.sqrt(14)),
                ZernikeCoefficient(30, "9th order trefoil at 0°", 4),
                ZernikeCoefficient(31, "9th order trefoil at 90°", 4),
                ZernikeCoefficient(32, "9th order astigmatism at 0°", np.sqrt(18)),
                ZernikeCoefficient(33, "9th order astigmatism at 45°", np.sqrt(18))
                )
