from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyThermoFlash'
VERSION = '0.1.1'
AUTHOR = 'Sina Gilassi'
EMAIL = '<sina.gilassi@gmail.com>'
DESCRIPTION = "PyThermoFlash is a Python package for performing thermodynamic flash calculations and property estimations for various fluid systems."
LONG_DESCRIPTION = "PyThermoFlash is a comprehensive computational tool designed for engineers, researchers, and scientists working with thermodynamic systems. The package focuses on flash calculations - a fundamental procedure in process simulation that determines the equilibrium state of a mixture at specified conditions. PyThermoFlash implements various thermodynamic models to accurately predict phase behavior and properties of complex fluid systems."

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    include_package_data=True,  # Make sure to include non-Python files
    # Add both config and data files
    package_data={'': ['config/*.yml', 'data/*.csv',
                       'templates/*.html', 'static/*']},
    # Add license file
    license='MIT',
    license_files=[],
    install_requires=['pandas', 'requests',
                      'urllib3', 'numpy', 'PyYAML', 'scipy', 'pycuc',
                      'PyThermoModels', 'PyThermoLinkDB'],
    extras_require={},
    keywords=['python', 'chemical engineering', 'thermodynamics',
              'flash-calculations', 'equilibrium', 'vle-calculations',
              'raoult-law', 'phase-equilibrium', 'dew-point-calculation',
              'bubble-point-calculation', 'vapor-liquid-equilibrium',],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.10",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.10',
)
