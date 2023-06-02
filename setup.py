#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn

# setup

from setuptools import setup, find_packages

setup(
    name='mahtpy',
    version='1.0.0',
    description='Draw fancy Manhattan plot with matplotlib!',
    author='lilinxuan',
    author_email='lilinxuan@genomics.cn',
    url='https://github.com/yourusername/mypackage',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib>=3.5.1'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)