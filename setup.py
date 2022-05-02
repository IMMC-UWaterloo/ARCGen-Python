from setuptools import setup, find_packages

with open("README.md","r", encoding="utf-8") as fh:
      long_description = fh.read()

setup(name='arcgen-python',
      version='2022.1',
      description='Arc-length Response Corridor Generation',
      long_description= long_description,
      long_description_content_type="text/markdown",
      author='Devon Hartlen',
      author_email='hartlendc@gmail.com',
      url="https://github.com/IMMC-UWaterloo/ARCGen-Python",
      license='GPL GNU v3',
      packages=find_packages(include=['arcgen', 'arcgen.*']),
      keywords=['averaging', 'statistics', 'corridors', 'arc-length', 
            'response corridors'],
      install_requires=[
            'numpy>=1.22.0',
            'scipy>=1.8.0',
            'matplotlib>=3.5.1'],
      python_requires='~=3.8',
      project_urls={
            'Documentation': 'https://github.com/IMMC-UWaterloo/ARCGen-Python',
            'Source': 'https://github.com/IMMC-UWaterloo/ARCGen-Python',
            'tracker': 'https://github.com/IMMC-UWaterloo/ARCGen-Python/issues'},
      classifiers=[
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Topic :: Scientific/Engineering"],
      )