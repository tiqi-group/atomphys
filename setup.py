from setuptools import find_packages, setup
setup(
    name='TIQIatomphys',
    packages=find_packages(),
    package_data = {'': ['data/*.json'],},
    version='0.1.0',
    description='',
    author='Matt Grau with Wojteks edits',
    license='MIT',
)
