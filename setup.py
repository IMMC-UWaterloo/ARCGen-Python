from setuptools import setup, find_packages

setup(name='arcgen-python',
      version='2022.1',
      description='Arc-length Response Corridor Generation',
      author='Devon Hartlen',
      author_email='hartlendc@gmail.com',
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
      )