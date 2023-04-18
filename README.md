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




