from setuptools import setup
from setuptools.command.install import install
from distutils.errors import DistutilsSetupError
import shutil

class CustomInstallCommand(install):
    """Customized setuptools install command - checks for SNID installation."""
    def run(self):
        snid_install = shutil.which("snid")
        if snid_install is not None:
            install.run(self)
        else:
            raise DistutilsSetupError("Cannot find SNID, aborting pySNID installation")

setup(name='pySNID',
      version='0.1.dev0',
      description='python code for running the SuperNova IDentifcation code (SNID) [Blondin and Tonry 2007]',
      url='https://github.com/benstahl92/pySNID',
      author='Benjamin Stahl',
      author_email='benjamin_stahl@berkeley.edu',
      license='MIT',
      packages=['pySNID'],
      install_requires=['numpy'],
      cmdclass = {
          'install' : CustomInstallCommand,
      },
    )
