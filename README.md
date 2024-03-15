
Import ConfAnalyser and run from another python script:
```python
from ConfAnalyser import ConfAnalyser

result = ConfAnalyser(paths_file=paths_to_pdbs_path, names_file=atom_names_path,
             molecule_type=ConfAnalyser.MoleculeType.Cyclohexane).result()
```

Run as a command-line application:
`python3 ConfAnalyser.py -l paths_to_pdbs -n atom_names --(cyclohexane|cyclopentane|benzene|oxane)`

# Electron density coverage analysis
Determines electron density coverage of a cycle. By default, uses
trilinear interpolation to infer the electron density value, and
considers an atom to be covered by the electron density when the
corresponding intensity is MORE than the threshold for the isosurface     
(1.5 sigma).

## Prerequisites
Install gemmi version 0.4.5, e.g. using `pip`:
```
pip install gemmi==0.4.5
```

## Arguments description

### Positional arguments
  - `input_cycle_pdb`      Input PDB file with coordinates of a cycle, produced by `PatternQuery`
  - `input_density_ccp4`   Input electron density file for the corresponding protein structure from the PDB in CCP4 format      

### Optional arguments
  - `-s`                   Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle
  - `-d`                   Detailed mode - output a sequence of paried values, where the first member of each pair is a serial atom number, and the second member is "y" - if atom is covered, and "n" - if it is not.      
  - `-m`, `--more_or_equal`  Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface       
  - `-c`, `--closest_voxel`  Instead of trilinear interpolation, the intensity of the closest voxel is used


## Example

Running the the script in `simple` mode with `--more_or_equal` setup:
```
python electron_density_coverage_analysis.py -sm 3biu_0.pdb 3biu.ccp4
```
