# -*- coding: utf-8 -*-
"""Setup module."""
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def get_requires():
    """Read requirements.txt."""
    requirements = open("requirements.txt", "r").read()
    return list(filter(lambda x: x != "", requirements.split()))


def read_description():
    """Read README.md and CHANGELOG.md."""
    try:
        with open("README.md") as r:
            description = "\n"
            description += r.read()
        with open("CHANGELOG.md") as c:
            description += "\n"
            description += c.read()
        return description
    except Exception:
        return '''OPR is an open-source Python package designed to simplify and streamline primer design and analysis for biologists and bioinformaticians. OPR enables users to design, validate, and optimize primers with ease, catering to a wide range of applications such as PCR, qPCR, and sequencing. With a focus on user-friendliness and efficiency, OPR aims to bridge the gap between biological research and computational tools, making primer-related workflows faster and more reliable.'''


setup(
    name='opr',
    packages=[
        'opr', ],
    version='0.1',
    description='OPR: Optimized Primer',
    long_description=read_description(),
    long_description_content_type='text/markdown',
    author='OPR Development Team',
    author_email='opr@openscilab.com',
    url='https://github.com/openscilab/opr',
    download_url='https://github.com/openscilab/opr/tarball/v0.1',
    keywords="primer biology bioinformatics genome dna pcr",
    project_urls={
            'Source': 'https://github.com/openscilab/opr',
    },
    install_requires=get_requires(),
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='MIT',
)
