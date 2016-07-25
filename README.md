# BioRuler

**BioRuler** is a library developed on top of the **ReGraph** library, which serves as a bridge between biological data formats describing protein-protein interactions and graph-based representation of knowledge.

Our main motivation while developing this library was to enable the user to automatically convert protein-protein interaction data in any standard format to its graph representation following a domain-specific meta-model of choice. We would also like to give the user the possibility to manipulate this graph representation, and finally, to automatically export it to the respective Kappa code, ready for simulations, static analysis, etc.

## Environment configs 

### Requirement

The following `Python 3` packages are required

```
cycler==0.10.0
decorator==4.0.10
JPype1==0.6.1
matplotlib==1.5.1
networkx==1.11
numpy==1.11.1
pyparsing==2.1.5
python-dateutil==2.5.3
pytz==2016.4
ReGraph==0.1
six==1.10.0
wheel==0.24.0
```

To avoid manual installation and to easily set up development environment you may consider following the instructions below:

### Create virtual environment

```
virtualenv venv -p path/to/your/python3
```

### Setup environment

To activate the virtual environment
```
source venv/bin/activate
```

To install required dependencies
```
pip install -r requirement.txt
```

## Installation

In order to install the **ReGraph** library you have to clone this repository using SSH
```
git clone git@github.com:eugeniashurko/ReGraph.git
```
or using HTTPS
```
https://github.com/eugeniashurko/ReGraph.git
```
Install the library with
```
python setup.py install
```

## Run BioPAXImporter tests

The tests are run on the **PID: the Pathway Interaction Database**

```
python test.py
```
