# Setup
- tested with Python 3.7

```
# in project's root directory
virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

if you dont have jupyter notebook already installed you can install it for your python, or for this virtualenv
pip install notebook # (either in the venv or not)

# Try me
Launch `jupyter notebook`

__or__ check the examples in if __main__  IDE and run the modules. They import from other packages of this library (relative imports), so
 you have to run the modules with respect to the root package (root of the python library) -- '-m flag for running directly the module
 ''`python -m assignments.path.to.module ` e.g. 'assignments.structure_related_properties.lib'  For this to work be at the _repository_
root  
# Code

for now, python modules have the testing code in `if __name__ == '__main__'` blocks, run only when the python file is run directly, not
 just imported from another python module
 
Fasta and PDB files are parsed using BioPython library. 

Clustal I parsed manually.
 
high-level description, location in directory structure
## 1 Edit Distance (Sequence alignment using edit distance)

## 2 Fasta (Processing FASTA files)

## 3 Hamming (Measuring sequence similarity using Hamming distance)

## 4 Clustal Parser (Processing multiple sequence alignment)

## 6 Pdb (Processing PDB files)
- I had to 'hack' the PDB library, so that it could work with my subclasses of Structure/Model/Chain/Residue
    - I patched (replaced) the classes used in the Bio.parser's StructureBuilder with my own subclasses
- The Bio.PDB library adopts Structure/Model/Chain/Residue/Atom scheme. Where structure can have multiple models (e.g. an NMR structure)

## 5 Conservation determination from multiple aligned sequences
## 7 Computing structure-related properties
## 8 combinining
## 9 superposition combinining


The 'frontend' for this library is the python (possible ipython notebooks). It's a good scripting language already and writing just a
 wrapper for shell commandline execution doesn't make sense for me.
 
 n papíře doma, co mám ještě udělat ze web. zadání. Udělám to znova

Follows a list of requirements and recommendations. please go through them before you submit your code:

    Provide a README (README.md in case of git repository) with a high-level description of the architecture of the library and with instructions about how to compile the program (if necessary) and how to run the individual tasks which you implemented.
   [done] Provide test data for all the implemented tasks and in README-.
   [document structure here, check documentation in code] The structure of the program should include separate directories with source code
   , with test data, documentation and so on. The directory structure has to be documented as well.
    The code has to be commented.
    
    [done and why should see higher conservation] If you decide to implement the "Combining structure and sequence" task, include also short
     analysis whether on your data you  indeed see higher conservation for active site residues.