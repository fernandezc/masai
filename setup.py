from setuptools import setup 
from masai.version import get_git_version


setup(
    name='masai',
    version=get_git_version(),
    packages=['masai', 'masai.core', 'masai.plugins'],
    install_requires=['nmrglue', 'numpy', 'traits', 'matplotlib', 'addict', 'pandas'],
    url='',
    license='',
    author='Christian Fernandez',
    author_email='christian.fernandez.14440@gmail.com',
    description=''
)
