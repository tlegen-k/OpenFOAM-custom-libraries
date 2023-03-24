# OpenFOAM-custom-libraries
This repository contains the computational modeling solver and associated libraries I have developed in my PhD journey.
All code was developed for OpenFOAM-v1812, nevertheless should be portable to more modern versions of software. The code was tested on macOS Monterey 12.1, CentOS v7, and generally should be easily portable to any Linux distribution (as long as OpenFOAM itseld was correctly installed). 
## interPlasmaFoam
Main solver `interPlasmaFoam` models turbulent gas flow, gas-liquid dynamic interface, multiphase reactive chemical species transport.
Volume-of-Fluid (VoF) approach is employed for capturing the multiphase flow and is coupled with the Chien's low-Reynolds k-$\epsilon$ turbulence model.

Multiphase species transport is modeled using 3 implemented approaches: 
- Haroun's model employing Harmonic averaging of diffusivities
- Marschall's model employing Arithmetic averaging of diffusivities
- The Continuous Species Transport model developed by Deising et al. Essentially, blends both approaches

## Custom libraries
`TurbulenceModels`

`transportModels`

`reactionMultiphase`

`multiphaseChemistryModel`