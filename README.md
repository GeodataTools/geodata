# Geodata

[![DOI](https://zenodo.org/badge/218690319.svg)](https://zenodo.org/badge/latestdoi/218690319)

**Geodata** is a Python library of geospatial data collection and
analytical tools. Through the creation of shared scripts and
documentation for analysis-ready physical variables, geodata streamlines
the collection and use of geospatial datasets for natural science,
engineering, and social science applications.

![png](docs/source/_static/images/geodata_workflow_chart.png)


# Motivation
The main motivation is the difficulty in working with high temporal and
spatial resolution datasets of physical variables from earth system
models and combining them with GIS datasets (land use, geographic
features, etc.). The primary analytical questions addressed here are
generating profiles of energy and environmental variables of interest
(solar PV, wind power, pollution distribution) subject to suitability and
weighting criteria. Additional applications are under development.

Working with these datasets has startup costs and computational barriers
due to diverse sources, formats, resolutions, and large memory
requirements. To solve this, geodata provides an all-in-one Python
interface for downloading, subsetting, and transforming large earth
systems datasets into relevant physical variables and flexibly
incorporating GIS datasets to mask these variables and generate
“analysis-ready” datasets for use in regression, plotting, and energy
models. Geodata builds off the structure of the [**atlite**](https://github.com/PyPSA/atlite) package.

# Installation
**Geodata** has been tested to run with python3 (>= 3.10). Read the [package setup instructions](https://geodata.readthedocs.io/en/latest/quick_start/packagesetup.html) to configure and install the package.

In short, you can install the latest version of the package by running the following commands:

```bash
pip install geodata-re
```

# Documentation
You can find detailed documentation for Geodata [here](https://geodata.readthedocs.io/en/latest/). This includes installation instructions, tutorials, and API documentation.

# Contributing
We welcome suggestions for feature enhancements and the identification of bugs. Please make an issue or contact the [authors](https://pwrlab.org/about.html) of geodata.


# License
Geodata is licensed under the GNU GENERAL PUBLIC LICENSE Version 3 (2007). This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](LICENSE) for more details.

# Support
The Geodata team would like to thank the Center for Global Transformation at UC San Diego for providing financial support to the project.
