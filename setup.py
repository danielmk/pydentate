from setuptools import setup


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pyDentate',
    version='0.1.1',
    description='pyDentate is a biophysically realistic computational model of the dentate gyrus.',
    long_description=readme,
    author='Daniel Müller-Komorowska',
    author_email='danielmuellermsc@gmail.com',
    url='https://github.com/danielmk/pyDentate',
    license=license,
    packages=['ouropy', 'pyDentate'],
    install_requires=[
          'elephant',
          'numpy',
          'pandas',
          'matplotlib',
          'scipy',
          'neo',
          'neuron ; platform_system=="Linux"',
          'scikit-image'
      ])
