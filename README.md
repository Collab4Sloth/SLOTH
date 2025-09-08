SLOTH 
=====

[![Build](https://github.com/Collab4Sloth/SLOTH/actions/workflows/build-and-tests.yml/badge.svg?color=green)](https://github.com/Collab4Sloth/SLOTH/actions/workflows/build-and-tests.yml)
![C++17](https://img.shields.io/badge/C%2B%2B-17|20-cyan.svg)
[![License](https://img.shields.io/github/license/Collab4Sloth/SLOTH.svg)](https://github.com/Collab4Sloth/SLOTH/blob/main/LICENSE)

_____________

Phase-field methods represent a versatile and effective approach to modelling the spatiotemporal evolution of complex physical systems that exhibit significant heterogeneities, such as phase transitions, coalescence, and cracking. These methods have found extensive application in the field of materials science, including in the modelling of the behaviour of nuclear fuels. 

The phase-field approach has also been employed in recent studies conducted with the multiphysics computational tools of the PLEIADES platform. Broadly speaking, the `PLEIADES`/Fuel Performance Codes aim at providing a state-of-the-art multiphysics multiscale description of the fuel under irradiated conditions. The inclusion of advanced multiphysics and multiscale modelling such as phase-field, is a logical step forward.
Consequently, the CEA is developing **`SLOTH`**, a **multiphase-field multicomponent framework**, dedicated to the study of fuel behaviour at different scales of description, from nominal operating conditions to severe accident conditions. 

`PLEIADES/SLOTH` is developed at CEA under LPGL license (Version 3). It is based on the the `MFEM` Finite Element library, which already includes the main features to have a robust, flexible and **massively parallel** implementation of the solution algorithms.

# Some useful links

- Documentation
  - [Website](https://collab4sloth.github.io/Documentation/) 
  - [GitHub repository](https://github.com/Collab4Sloth/Documentation) 

- [Gallery](https://github.com/Collab4Sloth/Gallery)
- [Sloth project and Team members](https://collab4sloth.github.io/Documentation/About/index.html)
- [Installation guidelines](https://collab4sloth.github.io/Documentation/Started/Installation/linux.html)

# Community Guidelines 

For more details, see [CONTRIBUTING.md](CONTRIBUTING.md). Main guidelines are:

- For any bug, please create an issue and add the label "bug". We welcome all feedback to make `SLOTH` as robust as possible.
- If you would like to participate and add functionality to `SLOTH`, you can find instructions for coding style, tests and pull request process in [CONTRIBUTING.md](CONTRIBUTING.md).
- If you have any support-related / collaboration questions, please contact the team at clement.introini@cea.fr.


# Contributors

-   __Modelling & Development Team__
  
    ---    
    - [Clément Introïni](https://www.researchgate.net/profile/Clement-Introini) (Phase-Field, Computer Science, Material Science)
    - [Raphaël Prat](https://www.researchgate.net/profile/Raphael-Prat) (Computer Science, HPC)



-   __Students Team__

    ---
    - [Alessandro Scapini (PhD 2024-2027)]()
    - [Clément Plumecocq (PhD 2023-2026)]()
    - [Victor Navarre (Master 2025)]()
    - [Jules Czuckermand (Master 2025)]()     
    - [Mouad Bakhkakh (Master 2024)]()
    - [Etienne Delobre (Master 2023)]()


# Acknowledgment

`SLOTH` is part of the `PLEIADES` platform jointly developed by the French Atomic Energy Commission in collaboration with its industrial partners, mainly EDF and Framatome.

# License
The code is developed under [GNU LESSER GENERAL PUBLIC LICENSE Version 3](LICENSE)
