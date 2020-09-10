#!/usr/bin/env sh

# ipython kernel install --user --name=projectname

for script in "edit_distance/distance.py"\
  "fasta/fasta.py"\
  "hamming/hamming.py"\
  "clustal_parser/clustal_parser.py" \
  "pdb/pdb.py"; do

  printf %${#script}s____'\n' | tr " " "_"
  echo "| $script |"
  printf %${#script}s===='\n' | tr " " "="

  grep -E "if __name__ == '__main__':" -A 1000 "$script"
  echo "****OUTPUT*************************************"
  "$script"
  echo ""
  echo ""
done