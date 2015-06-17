from setuptools import setup

from masai.version import get_version


setup(
    name='masai',
    version=get_version(),
    packages=['masai', 'masai.core', 'masai.plugins'],
    install_requires=['numpy', 'scipy', 'matplotlib', 
					  'traits', 'traitsui', 'apptools',
                      'pandas'],
    url='',
    license='',
    author='Christian Fernandez',
    author_email='christian.fernandez0@ensicaen.fr',
    description=''
)
