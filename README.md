
Import ConfAnalyser and run from another python script:
```python
from ConfAnalyser import ConfAnalyser

result = ConfAnalyser(paths_file=paths_to_pdbs_path, names_file=atom_names_path,
             molecule_type=ConfAnalyser.MoleculeType.Cyclohexane).result()
```

Run as a command-line application:
`python3 ConfAnalyser.py -l paths_to_pdbs -n atom_names --(cyclohexane|cyclopentane|benzene|oxane)`

