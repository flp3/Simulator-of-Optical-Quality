# Simulator of optical systems quality


This project aims to democratize physical and mathematical concepts that various industries uses to design optical systems, from smartphones to Eath Observation satellites.

With this python project, the user could:
* Simulate different optical systems, as telescopes or lenses, with key parameters, i.e. aperture radius, secondary mirror obstruction, focal length, pixel size, light wavelength.
* Generate optical quality metrics visualizations, i.e. Point spread Function, Modulation Transfer Function.
* Add and understand optical aberrations which are inherent to any acquisition system manufacturing.


## Acknowledgements
In this project we are using physical and mathematical concepts developped in the past 2 centuries:

From [Joseph von Fraunhofer](https://en.wikipedia.org/wiki/Joseph_von_Fraunhofer) who helped to modelize the propagation of the light, as a wave, through a partly blocked plane, with the Fraunhofer diffraction equation.
![Fraunhofer_diffraction_pattern_image](https://github.com/flp3/Simulator-of-Optical-Quality/tree/master/Screenshots/light propagation through an aperture.png?raw=true)

To [Frits Zernike](https://www.nobelprize.org/prizes/physics/1953/zernike/facts/), who was awarded with the Nobel Prize for Physics in 1952 for his work on modelizing optical aberrations through phase contrast polynomial.

The code is mainly based on those 3 scientific papers / books:
 - To simulate light diffraction this book helped me a lot, [Modeling the Imaging Chain of Digital Cameras](https://www.spiedigitallibrary.org/ebooks/TT/Modeling-the-Imaging-Chain-of-Digital-Cameras/eISBN-9780819483362/10.1117/3.868276?SSO=1)
 - To test that the algorithms modelize correctly the physical phenomena, this paper helped me, [The Use Of Image Quality Criteria In Designing A Diffraction Limited Large Space Telescope](https://spie.org/Publications/Proceedings/Paper/10.1117/12.953525?SSO=1)
 - To add optical aberrations with the simple and beautiful Zernike polynomial, this book helps me, [Zernike F, The Diffraction Theory of Aberrations, in Optical Image Evaluation Circular 526, ( National Bureau of Standards , Washington, D. C), 1952](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular526.pdf)

## Installation

As a recommendation, you could install this project with a virtual environment as [venv](https://docs.python.org/3/library/venv.html)

```bash
  git clone https://github.com/flp3/Simulator-of-Optical-Quality.git
  cd Simulator-of-Optical-Quality

  pip install virtualenv
  python -m venv venv
  venv\Scripts\activate

  pip install .
```
    
## Demo
For anyone interested into simulating and getting nice visualization without getting into the code, you could find the command line interface useful: [cli.py](https://github.com/flp3/Simulator-of-Optical-Quality/blob/master/src/cli.py)

Some nice commands are already deployed:
```bash
python .\src\cli.py --help
```
```bash
Options:
  -h, --help  Show this message and exit.

Commands:
  visualize-aberrations                    visualize the optical aberrations
  visualize-apertures                      visualize the instrument aperture
  visualize-diffraction                    visualize the light diffraction through an aperture
  visualize-system-with-aberrations        visualize the light propagation through an
                                             acquisition system with optical aberrations
```
For example we could visualize the light propagation through a lens of 15mm of radius and see the Point spread function, and its Modulation function:
```bash
python .\src\cli.py visualize-system-with-aberrations l --radius 15 --wavelength 800e-9 --focal_length 0.1 --pixel_size 1.5e-6
```
![light propagation example](https://github.com/flp3/Simulator-of-Optical-Quality/blob/master/Screenshots/light%20propagation%20with%20a%20lense,%2015mm%20radius,%20wavelength%20800nm,%20focal%20length%2010cm,%20pixel%20size%201.5um.png?raw=true)
