# MATLAB Code Examples of Simulation and Compuational Inference for Stochastic Biochemical Reaction Networks

This repository contains useful example MATLAB functions and scripts as an introduction to 
stochastic simulation and computational inference methods applied to stochastic models of biochemical reaction networks.

## Developer

David J. Warne (david.warne@qut.edu.au),
                School of Mathematical Sciences, 
                Science and Engineering Faculty, 
                Queensland Univeristy of Technology 
                
Google Scholar: (https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en)

## Citation Information

This code is provided as supplementary information to the paper,

David J Warne, Ruth E Baker, and Matthew J Simpson. Implementations of simulation 
and inference algorithms for stochastic reaction networks: from basic concepts 
to state-of-the-art Journal XXXX, X:pp--pp, 2019. 

## Licensing
This source code is licensed under the GNU General Public License Version 3.
Copyright (C) 2018 David J. Warne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Contents

This folder contains a number of instructive MATLAB implementations 
for Stochastic simulation and computational inference. Demonstration scripts
showing typical usage are also provided
```bash
The directory structure is as follows
|-- init.m                                  Adds all functions to the MATLAB Path
|-- Functions
    |-- BCRN_Definition                     Creation of BCRN data structures
        |-- MichaelisMenten.m
        |-- MonoMolecularChain.m
    |-- Simulation                          ESSA and ASSA methods and Monte Carlo
        |-- GillespieDirectMethod.m
        |-- ModifiedNextReactionMethod.m
        |-- TauLeapingMethod.m
        |-- CorTauLeapingMethod.m
        |-- CMEsolMonoMol.m
        |-- MonteCarlo.m
        |-- MonteCarloTauLeap.m
        |-- MonteCarloBiasCorrection.m
        |-- MultilevelMonteCarlo.m
    |-- Inference                           Bayesian and ABC samplers
        |-- GenerateObservations.m
        |-- likelihoodCMEsolMOnoMol.m
        |-- ABCRejectionSampler.m
        |-- ABCRejectionSamplerV2.m
        |-- ABCMCMCSampler.m
        |-- ABCSMCSampler.m
        |-- ABCMLMC.m
|-- Demostrations
    |-- Simulation                          Demonstration of simulation algorithms
        |-- DemoGillespie.m
        |-- DemoMNRM.m
        |-- DemoTauLeap.m
        |-- DemoRealisationsMonoMol.m
        |-- DemoRealisationsMichMent.m
        |-- DemoCMEMeanVar.m
        |-- DemoStationaryDist.m
        |-- DemoCorTauLeap.m
        |-- DemoMonteCarlo.m
    |-- Inference                           Demonstration of inference algorithms
        |-- DemoDirectBayesCME.m
        |-- DemoABCConvergence.m
        |-- DemoABCProcess.m
        |-- DemoABCMethodsMonoMol.m
        |-- DemoABCMethodsMichMent.m
```
## Usage

Follow these steps to run the demonstrations:

1. Start MATLAB
2. In MATLAB browse to the repository folder Warne2018
3. In the MATLAB command prompt, run 
   `>> init` 
   to set up the paths of all the example implementations.
4. Type in the name of any demo script in Demonstrations/Simulation 
   or Demonstrations/Inference. For example,
   `>> DemoGillespie` 
   generates Fig. 1A and 1B in the paper.

## List of examples

The following list of examples shows how to reproduce the figures in the main paper. For more computationlly intensive examples approximate run times are given for an Intel(R) Core(TM) i7-5600U CPU (2.6 GHz).

### Figure 1

For Fig. 1A and 1B
`>> DemoGillespie`

For Fig. 1C and 1D
`>> DemoMNRM`

For Fig. 1E and 1F
`>> DemoTauLeap`

For Fig. 1G
`>> DemoRealisationsMichMent`

For Fig. 1H
`>> DemoRealisationsMonoMol`

To obtain Fig. 1I (resp. Fig. 1J), `edit DemoRealisationsMichMent.m` 
(resp. `DemoRealisationsMonoMol.m`) and change 
`N` to `100` (line 17) and `alpha` to `0.05` (line 24), and re-run the script.

### Figure 2

For mean and variance plot in Fig. 2A 
`>> DemoCMEMeanVar`

To plot stationary distribution in F. 2B and the full CME solution
in Fig. 2C--2F run (Warning! this will take about 10 minutes)
`>> DemoStationaryDist`

### Figure 3

For Fig. 1A
`>> DemoCorTauLeap`

For Fig. 1B (Warning! takes about 24 hours)
`>> DemoMonteCarlo`

### Figure 4

For Fig. 4A--4C (Warning! takes about 2.5 hours)
`>> DemoABCConvergence`

### Figure 5

For Fig. 5A--5N (Warning! takes about 3.5 hours)
`>> DemoABCprocess`

### Figure 6 and tables 1 and 2

For Fig. 6A--6C and table 1 (Warning! takes about 5 hours) 
`>> DemoABCMethodsMonoMol`

For Fig. 6D--6F and table 2 (Warning! takes about 1.5 hours)
`>> DemoABCMethodsMichMent`
