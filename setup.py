from setuptools import setup


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pydentate',
    version='0.1.1',
    description='pydentate is a biophysically realistic computational model of the dentate gyrus.',
    long_description=readme,
    author='Daniel Müller-Komorowska',
    author_email='danielmuellermsc@gmail.com',
    url='https://github.com/danielmk/pydentate',
    license=license,
    packages=['ouropy', 'pydentate'],
    install_requires=[
          'elephant',
          'numpy',
          'pandas',
          'matplotlib',
          'scipy',
          'neo',
          'neuron ; platform_system=="Linux"',
          'scikit-image'],
    package_data={'pydentate': ['*.txt', 'win64\\*', 'x86_64\\*']})
