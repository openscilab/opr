<div align="center">
    <img src="https://github.com/openscilab/opr/raw/main/otherfiles/logo.png" width="250">
    <h1>OPR: Optimized Primer</h1>
    <br/>
    <a href="https://codecov.io/gh/openscilab/opr"><img src="https://codecov.io/gh/openscilab/opr/branch/dev/graph/badge.svg" alt="Codecov"></a>
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
<table>
    <tr> 
        <td align="center">Code Quality</td>
        <td align="center"><a href="https://app.codacy.com/gh/openscilab/opr/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade"><img src="https://app.codacy.com/project/badge/Grade/0a819f6eb6ae483695ad4934eff42df9"></a></td>
        <td align="center"><a href="https://www.codefactor.io/repository/github/openscilab/opr"><img src="https://www.codefactor.io/repository/github/openscilab/opr/badge" alt="CodeFactor"></a></td>
    </tr>
</table>


## Installation

### PyPI
- Check [Python Packaging User Guide](https://packaging.python.org/installing/)
- Run `pip install opr==0.5`
### Source code
- Download [Version 0.5](https://github.com/openscilab/opr/archive/v0.5.zip) or [Latest Source](https://github.com/openscilab/opr/archive/dev.zip)
- Run `pip install .`

## Usage

### Load
```pycon
>>> from opr import Primer, MeltingTemperature
>>> primer1 = Primer(sequence="CTGGAGGACGGAAGAGGAAGTAA", salt=50)
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
#### GC clamp
```pycon
>>> primer1.gc_clamp
1
```
#### Single run length
```pycon
>>> primer1.single_runs
{'A': 2, 'T': 1, 'C': 1, 'G': 2}
```
#### Double run length
```pycon
>>> primer1.double_runs
{'TA': 1, 'TC': 0, 'TG': 1, 'AT': 0, 'AC': 1, 'AG': 2, 'CT': 1, 'CA': 0, 'CG': 1, 'GT': 1, 'GA': 1, 'GC': 0}
```
#### Repeats
```pycon
>>> primer1.repeats(sequence="GG", consecutive=False)
4
```
```pycon
>>> primer1.repeats(sequence="GG", consecutive=True)
1
```
#### Melting temperature
##### Basic
```pycon
>>> primer1.melting_temperature()
57.056521739130446
>>> primer1.melting_temperature(MeltingTemperature.BASIC)
57.056521739130446
```
##### Salt-adjusted
```pycon
>>> primer1.melting_temperature(MeltingTemperature.SALT_ADJUSTED)
64.64203250676053
```
##### Nearest neighbor
```pycon
>>> primer1.melting_temperature(MeltingTemperature.NEAREST_NEIGHBOR)
66.42693710590595
```
#### Thermodynamic Constants
##### Enthalpy change (ΔH)
```pycon
>>> primer1.delta_h
-174.99999999999997
```
##### Entropy change (ΔS)
```pycon
>>> primer1.delta_s
-0.44210000000000005
```
#### Extinction Coefficient at 260 nm (E260)
```pycon
>>> primer1.E260
248.40000000000006
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

#### To RNA
```pycon
>>> oprimer_rna = oprimer.to_rna()
>>> oprimer_rna
'CUGGAGGACGGAAGAGGAAGUAA'
```

#### To protein
```pycon
>>> oprimer_protein = oprimer.to_protein(frame=3)
>>> oprimer_protein
'GGRKRK*'
>>> oprimer_protein_aa3 = oprimer.to_protein(frame=3, multi_letter=True)
>>> oprimer_protein_aa3
'Gly-Gly-Arg-Lys-Arg-Lys-Stop'
```

ℹ️ When the `frame=3` it starts from the third nucleotide. The list of codons is then [`GGA`(G), `GGA`(G), `CGG`(R), `AAG`(K), `AGG`(R), `AAG`(K), `TAA`(STOP)]

ℹ️ Stop signal is marked as `*` in OPR

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