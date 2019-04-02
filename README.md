# BeamParameterEvo
Evolution of parameters during (HL-)LHC collision process.

## Description
Using the GSL library calculate the time evolution of beam & machine parameters during the collision process of (HL-)LHC. The emittances are evolved using the emittance model by F. Antoniou *et. al*, "LHC Luminosity Modeling for RUNII", Proceedings, 7th International Particle Accelerator Conference (IPAC 2016): Busan, Korea, May 8-13, 2016.

This work uses the statically compiled advantages and implements the [beamCal](http://www.github.com/nkarast/beamCal) code.

## Usage

The printout happens in the stdout which depending on the settings in `BeamParameters.hh` the evolution of the machine & beam parameters is followed.
$$\beta^*$$ and offset leveling are easily performed in terms of root finders, and the crossing angle (anti-)leveling is user based (usually Dynamic Aperture driven).

The code assumes the LHC ring, i.e. 4 Interaction Points (IPs). The two high luminosity experiments are always taken into account (IP1/5), while there is also the option of including the impact of IP2/8.

## Output

Based on the choices the user "ticked" in the `BeamParameters.hh` and the `.cpp` the code returs:
- time [h]
- $$\beta^*$$ in IP1/5 [m]
- Half-Crossing angle in IP1/5 [$\mu$rad]
- Integrated Luminosity in IP1/5 [$fm^{-1}$]
- Instantaneous Luminosity in IP1/5 [$10^{34}$ $Hz/cm^{2}$]
- Total Pileup in IP1/5
- r.m.s. luminous region in IP1/5 [cm]
- Peak Pileup in IP1/5 [evt/mm]
- Bunch Intensity [$10^{11}$ ppb]
- $\varepsilon_{n,x}$ [\mu$m]
- $\varepsilon_{n,y}$ [\mu$m]
- Bunch Length [m]
- $$\beta^*$$ in IP8 [m]
- Half-Crossing angle in IP8 [$\mu$rad]
- Parallel Separation in IP8 [$\sigma_{beam}$]
- Instantaneous Luminosity in IP8 [$10^{34}$ $Hz/cm^{2}$]
- Integrated Luminosity in IP8 [$fm^{-1}$]
- r.m.s. luminous region in IP8 [cm]
- Total Pileup in IP8
- Peak Pileup in IP8 [evt/mm]
- $$\beta^*$$ in IP2 [m]
- Half-Crossing angle in IP2 [$\mu$rad]
- Parallel Separation in IP2 [$\sigma_{beam}$]
- Instantaneous Luminosity in IP2 [$10^{34}$ $Hz/cm^{2}$]
- Integrated Luminosity in IP2 [$fm^{-1}$]
- r.m.s. luminous region in IP2 [cm]
- Total Pileup in IP2
- Peak Pileup in IP2 [evt/mm]


## Contact
N. Karastathis, nkarast .at. cern .dot. ch
