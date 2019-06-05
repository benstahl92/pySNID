from setuptools import setup

setup(name='pySNDB',
      version='0.1.dev0',
      description='python code for running the SuperNova IDentifcation code (SNID) [Blondin and Tonry 2007]',
      url='https://github.com/benstahl92/pySNID',
      author='Benjamin Stahl',
      author_email='benjamin_stahl@berkeley.edu',
      license='MIT',
      packages=['pySNID'],
      install_requires=['numpy']
    )