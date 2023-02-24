#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(name="ParticlePostProcessor",
      version="0.1",
      description="Post-processing for particle tracking model",
      packages=find_packages(),
      entry_points={
          "console_scripts": ["ppp=src.main:main"]
      }
      )
