from setuptools import setup, find_packages


setup(     
    name='Simulator of Optical Quality',
    author='Florian Pignol',
    author_email="florian.pignol2@gmal.com",
    version='0.1',
    package_dir={"": "src"},
    packages=find_packages(where="src"), 
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=[         
        'opencv-contrib-python',         
        'numpy',
        'matplotlib',
    ],
)