# Setup
- tested with Python 3.7

```
# in project's root directory
virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

# Try me
You can run:
```
./tryme.sh | less
```

__or rather__ open this project in your favourite IDE and run all .py files (excl. `structure-related_properties/lib.py` and `conservation
.py`
, they
 are work in progress)
 
# Code

for now, python modules have the testing code in `if __name__ == '__main__'` blocks, run only when the python file is run directly, not
 just imported from another python module
 
Fasta and PDB files are parsed using BioPython library. 

Clustal I parsed manually.
 
## 1 Edit Distance (Sequence alignment using edit distance)

## 2 Fasta (Processing FASTA files)

## 3 Hamming (Measuring sequence similarity using Hamming distance)

## 4 Clustal Parser (Processing multiple sequence alignment)

## 6 Pdb (Processing PDB files)
- I had to 'hack' the PDB library, so that it could work with my subclasses of Structure/Model/Chain/Residue
    - I patched (replaced) the classes used in the Bio.parser's StructureBuilder with my own subclasses
- The Bio.PDB library adopts Structure/Model/Chain/Residue/Atom scheme. Where structure can have multiple models (e.g. an NMR structure)

## 5 todo (Conservation determination from multiple aligned sequences)
## 7 todo (Computing structure-related properties)


The 'frontend' for this library is the python (possible ipython notebooks). It's a good scripting language already and writing just a
 wrapper for shell commandline execution doesn't make sense for me.