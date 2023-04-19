import click
from optical_system import zernike_polynomial as zp
from optical_system.optical_response import Aperture, Optical_psf


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    pass


APERTURES = {
    'l': 'lense',
    't': 'Cassgrain Telescope',
}


@cli.command(short_help='visualize the instrument aperture')
@click.argument("instrument", type=click.Choice(list(APERTURES.keys())), default='l')
@click.option("--radius",
              type=click.INT,
              default=150,
              help='To define the aperture radius, assuming aperture is circular. We prefer to express it in mm.'
              )
@click.option("--obstruction",
              type=click.INT,
              default=50,
              help='For Cassgrain Telescope, To define the obstruction radius, we prefer to express it in mm.'
              )
@click.option("--unit",
              type=click.FLOAT,
              default=1e-3,
              help='To define the unit of the aperture radius, we prefer to express it in mm.'
              )
def visualize_apertures(instrument, radius, obstruction, unit):
    if instrument == 'l':
        ap = Aperture.disk(radius, unit)
    if instrument == 't':
        ap = Aperture.disk_obstruction_spider(radius, obstruction, unit)
    ap.add_padding(2)
    ap.illustrate_magnitude()


@cli.command(short_help='visualize the light diffraction through an aperture')
@click.argument("instrument", type=click.Choice(list(APERTURES.keys())), default='l')
@click.option("--radius",
              type=click.INT,
              default=150,
              help='To define the aperture radius, assuming aperture is circular. We prefer to express it in mm.'
              )
@click.option("--obstruction",
              type=click.INT,
              default=50,
              help='For Cassgrain Telescope, To define the obstruction radius, we prefer to express it in mm.'
              )
@click.option("--unit",
              type=click.FLOAT,
              default=1e-3,
              help='To define the unit of the aperture radius, we prefer to express it in m.'
              )
@click.option("--wavelength",
              '-wl',
              type=click.FLOAT,
              default=560e-9,
              help='To define the wavelength considered, we prefer to express it in m.'
              )
@click.option("--focal_length",
              '-fl',
              type=click.FLOAT,
              default=0.9,
              help='To define the instrument focal length, we prefer to express it in m.'
              )
@click.option("--pixel_size",
              '-px',
              type=click.FLOAT,
              default=5e-6,
              help='To define the instrument pixel size, we prefer to express it in m.'
              )
def visualize_diffraction(instrument, radius, obstruction, unit, wavelength, focal_length, pixel_size):
    if instrument == 'l':
        ap = Aperture.disk(radius, unit)
    if instrument == 't':
        ap = Aperture.disk_obstruction_spider(radius, obstruction, unit)
    ap.add_padding(4 ** 2)
    psf = Optical_psf.from_aperture(ap, wavelength, focal_length)
    mtf_sampling, mtf = psf.mtf(unit=pixel_size)
    mtfs = [(mtf_sampling, mtf, 'no aberrations')]
    psf.illustrate_psf_and_mtfs(mtfs)


@cli.command(short_help='visualize the optical aberrations')
@click.option("--panel",
              type=click.BOOL,
              default=True,
              help=''
              )
@click.option('-za',
              '--zernike_azimuth',
              default=0,
              type=click.INT,
              help='')
@click.option('-zr',
              '--zernike_radial',
              default=2,
              type=click.INT,
              help='')
def visualize_aberrations(panel, zernike_azimuth, zernike_radial):
    zernike_coeffs = zp.ZernikeCoefficients(33)
    if panel:
        zernike_coeffs.illustrate_panel(zernike_radial, zernike_azimuth)
    else:
        zernike_coeffs.illustrate_single(zernike_radial, zernike_azimuth)


@cli.command(short_help='visualize the light propagation through an acquisition system with optical aberrations')
@click.argument("instrument", type=click.Choice(list(APERTURES.keys())), default='l')
@click.option("--radius",
              type=click.INT,
              default=150,
              help='To define the aperture radius, assuming aperture is circular. We prefer to express it in mm.'
              )
@click.option("--obstruction",
              type=click.INT,
              default=50,
              help='For Cassgrain Telescope, To define the obstruction radius, we prefer to express it in mm.'
              )
@click.option("--unit",
              type=click.FLOAT,
              default=1e-3,
              help='To define the unit of the aperture radius, we prefer to express it in mm.'
              )
@click.option("--wavelength",
              '-wl',
              type=click.FLOAT,
              default=560e-9,
              help='To define the wavelength considered, we prefer to express it in m.'
              )
@click.option("--focal_length",
              '-fl',
              type=click.FLOAT,
              default=0.9,
              help='To define the instrument focal length, we prefer to express it in mm.'
              )
@click.option("--pixel_size",
              '-px',
              type=click.FLOAT,
              default=5e-6,
              help='To define the instrument pixel size, we prefer to express it in m.'
              )
@click.option('-za',
              '--zernike_azimuth',
              default=0,
              type=click.INT,
              help='')
@click.option('-zr',
              '--zernike_radial',
              default=2,
              type=click.INT,
              help='')
@click.option('--aberration_weight',
              '-aw',
              default=0.1,
              type=click.FLOAT,
              help='between o and 1, to give more or less impact of the aberration')
def visualize_system_with_aberrations(instrument,
                                      radius,
                                      obstruction,
                                      unit,
                                      wavelength,
                                      focal_length,
                                      pixel_size,
                                      zernike_azimuth,
                                      zernike_radial,
                                      aberration_weight):
    if instrument == 'l':
        ap = Aperture.disk(radius, unit)
    if instrument == 't':
        ap = Aperture.disk_obstruction_spider(radius, obstruction, unit)

    ap_no_aberration = ap.copy_with()
    ap_no_aberration.add_padding(4 ** 2)
    psf_no_aberration = Optical_psf.from_aperture(ap_no_aberration, wavelength, focal_length)
    no_abmtf_sampling, no_ab_mtf = psf_no_aberration.mtf(unit=pixel_size)

    zernike_coeffs = zp.ZernikeCoefficients(33)
    phase = zernike_coeffs[(zernike_radial, zernike_azimuth)].cartesian_matrix(ap.array.shape[0])[1]
    aberration_name = zernike_coeffs[(zernike_radial, zernike_azimuth)].name
    ap.add_phase(phase, aberration_weight)
    ap.add_padding(4 ** 2)
    psf = Optical_psf.from_aperture(ap, wavelength, focal_length)
    ab_mtf_sampling, ab_mtf = psf.mtf(unit=pixel_size)

    mtfs = [(no_abmtf_sampling, no_ab_mtf, 'no aberrations'),
            (ab_mtf_sampling, ab_mtf, aberration_name)
            ]

    extra_title = 'with ' + str(aberration_weight * 100) + '% ' + aberration_name + ' aberration'

    psf.illustrate_psf_and_mtfs(mtfs, sub_title=extra_title)


if __name__ == '__main__':
    cli()
