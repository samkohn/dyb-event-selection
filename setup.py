from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))

setup(
        name='dyb-analysis',
        description='Daya Bay analysis',
        packages=find_packages(),
)
