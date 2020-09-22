## Setup
Requires Python 3.7 or newer (definitely won't work with Python 3.4, 3.5 might).
```
# in repository root directory

virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

For the assignments dealing with protein surface, a molecular surface calculator library, used by BioPython, has to be installed. 
Go to [http://mgltools.scripps.edu/downloads#msms](http://mgltools.scripps.edu/downloads#msms) and select your binary. 
On linux, I had to extract the archive and for BioPython
 to find it, I linked the binary to `/usr/local/bin/msms` (should be in PATH under `msms`).
 
    sudo ln /path_to_extracted_archive/msms.x86_64Linux2.2.6.1 /usr/local/bin/msms1


If you don't have jupyter notebook already installed you can install it for your Python, or for this virtualenv:
    
    pip install notebook  # (run either in the venv or not)


## Overview

The repository includes a Python package at `/assignments`. It contains the library functions as well as source code for all of the
 assignments. 

Jupyter notebook `/assignment_demos.ipynb` contains solutions for all the assignments, shows that the required tasks are completed,
 comments on the solution and, if applicable, discusses the results.
 
The code in the notebook was compiled from the assignment solutions in the Python modules in `/assignments`. They are located in `if
 __name__ == '__main__':` blocks and also serve as example usage for the library functions in the module. 
 
All assignments have their own package/module, except assignment 'Conservation determination from multiple aligned sequences', that
 functionality is included in `/assignments/clustal_parser` ('Processing multiple sequence alignment').
 
Third party Python libraries used include BioPython (+ non-Python msms, molecular surface calculator), rmsd (Kabsch algorithm) and plotting
 libraries (matplotlib, seaborn).
 
## Run

Launch `jupyter notebook` and open `assignment_demos.ipynb` in the web interface

_or_ check the examples at the end of the sources and run the modules. They import from other packages of this library (relative imports
), so
 you have to run the modules with respect to the root package (root of the Python library). In the virtual environment, run `python -m
  assignments.path.to.module` e.g. `assignments.structure_related_properties.lib`.  Run this at the _repository_
root level.


The 'frontend' for this library is Python (possibly iPython notebooks). It's a good scripting language already and writing just a
 wrapper for shell commandline execution doesn't make sense for me.
