<div align="center">
    <img src="https://github.com/openscilab/opr/raw/main/otherfiles/logo.png" width="250">
    <h1>OPR: Optimized Primer</h1>
    <br/>
    <a href="https://codecov.io/gh/openscilab/opr"><img src="https://codecov.io/gh/openscilab/opr/branch/main/graph/badge.svg" alt="Codecov"></a>
    <a href="https://badge.fury.io/py/opr"><img src="https://badge.fury.io/py/opr.svg" alt="PyPI version"></a>
    <a href="https://www.python.org/"><img src="https://img.shields.io/badge/built%20with-Python3-green.svg" alt="built with Python3"></a>
    <a href="https://discord.gg/8mBspwXqcA"><img src="https://img.shields.io/discord/1064533716615049236.svg" alt="Discord Channel"></a>
</div>

----------


## Overview
<p align="justify">
<b>OPR</b> is an open-source Python package designed to simplify and streamline primer design and analysis for biologists and bioinformaticians. <b>OPR</b> enables users to design, validate, and optimize primers with ease, catering to a wide range of applications such as PCR, qPCR, and sequencing. With a focus on user-friendliness and efficiency, <b>OPR</b> aims to bridge the gap between biological research and computational tools, making primer-related workflows faster and more reliable.
</p>
<table>
    <tr>
        <td align="center">PyPI Counter</td>
        <td align="center">
            <a href="https://pepy.tech/projects/opr">
                <img src="https://static.pepy.tech/badge/opr">
            </a>
        </td>
    </tr>
    <tr>
        <td align="center">Github Stars</td>
        <td align="center">
            <a href="https://github.com/openscilab/opr">
                <img src="https://img.shields.io/github/stars/openscilab/opr.svg?style=social&label=Stars">
            </a>
        </td>
    </tr>
</table>
<table>
    <tr> 
        <td align="center">Branch</td>
        <td align="center">main</td>
        <td align="center">dev</td>
    </tr>
    <tr>
        <td align="center">CI</td>
        <td align="center">
            <img src="https://github.com/openscilab/opr/actions/workflows/test.yml/badge.svg?branch=main">
        </td>
        <td align="center">
            <img src="https://github.com/openscilab/opr/actions/workflows/test.yml/badge.svg?branch=dev">
            </td>
    </tr>
</table>


## Installation

### PyPI
- Check [Python Packaging User Guide](https://packaging.python.org/installing/)
- Run `pip install opr==0.1`
### Source code
- Download [Version 0.1](https://github.com/openscilab/opr/archive/v0.1.zip) or [Latest Source](https://github.com/openscilab/opr/archive/dev.zip)
- Run `pip install .`

## Usage

### Load
```pycon
>>> from opr import Primer
>>> primer1 = Primer(sequence="CTGGAGGACGGAAGAGGAAGTAA")
>>> primer1.sequence
'CTGGAGGACGGAAGAGGAAGTAA'
```

### Properties

#### Molecular weight
```pycon
>>> primer1.molecular_weight
7235.79
```
#### GC content
```pycon
>>> primer1.gc_content
0.5217391304347826
```
#### Melting temperature
```pycon
>>> primer1.melting_temperature()
57.056521739130446
```
### Operations

#### Reverse
```pycon
>>> primer1_reversed = primer1.reverse()
>>> primer1_reversed.sequence
'AATGAAGGAGAAGGCAGGAGGTC'
```
#### Complement
```pycon
>>> primer1_complemented = primer1.complement()
>>> primer1_complemented.sequence
'GACCTCCTGCCTTCTCCTTCATT'
```

## Issues & bug reports

Just fill an issue and describe it. We'll check it ASAP! or send an email to [opr@openscilab.com](mailto:opr@openscilab.com "opr@openscilab.com"). 

- Please complete the issue template
 
You can also join our discord server

<a href="https://discord.gg/8mBspwXqcA">
  <img src="https://img.shields.io/discord/1064533716615049236.svg?style=for-the-badge" alt="Discord Channel">
</a>

## References
					
<blockquote>1- <a href="http://biotools.nubic.northwestern.edu/OligoCalc.html">Oligo Calc: Oligonucleotide Properties Calculator</a></blockquote>

<blockquote>2- Marmur, Julius, and Paul Doty. "Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature." <i>Journal of molecular biology</i> 5.1 (1962): 109-118.</blockquote>

<blockquote>3- Wallace, R. Bruce, et al. "Hybridization of synthetic oligodeoxyribonucleotides to Φ X 174 DNA: the effect of single base pair mismatch." <i>Nucleic acids research</i> 6.11 (1979): 3543-3558.</blockquote>

<blockquote>4- Panjkovich, Alejandro, and Francisco Melo. "Comparison of different melting temperature calculation methods for short DNA sequences." <i>Bioinformatics 21.6</i> (2005): 711-722.</blockquote>


## Show your support


### Star this repo

Give a ⭐️ if this project helped you!

### Donate to our project
If you do like our project and we hope that you do, can you please support us? Our project is not and is never going to be working for profit. We need the money just so we can continue doing what we do ;-) .			

<a href="https://openscilab.com/#donation" target="_blank"><img src="https://github.com/openscilab/opr/raw/main/otherfiles/donation.png" height="90px" width="270px" alt="OPR Donation"></a>