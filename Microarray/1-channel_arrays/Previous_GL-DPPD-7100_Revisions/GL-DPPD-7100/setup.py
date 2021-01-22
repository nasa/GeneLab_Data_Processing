from setuptools import setup
import setuptools


setup(name='GeneLab-Microarray',
    version='0.6.4', 
    description='Standardized processing pipeline for microarray data on GeneLab.', 
    long_description='Standardized processing pipeline for microarray data on GeneLab.',
    url='https://github.com/jdrubin91/GeneLab-Microarray', 
    author='Jonathan Rubin & Daniel Mattox', 
    author_email='jonathan.d.rubin@nasa.gov', 
    license='GPL-3.0',
    packages=setuptools.find_packages(), 
    include_package_data=True,
    package_data={'': ['R_scripts/*','logos/*']},
    install_requires=['scipy','mpld3','jinja2','matplotlib','sklearn'],
    zip_safe=False,
    scripts=['bin/GeneLab-Microarray'])
