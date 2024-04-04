# ConfAnalyser

Import ConfAnalyser and run from another python script:
```python
from ConfAnalyser import ConfAnalyser

result = ConfAnalyser(paths_file=paths_to_pdbs_path, names_file=atom_names_path,
             molecule_type=ConfAnalyser.MoleculeType.Cyclohexane).result()
```

Run as a command-line application:
`python3 ConfAnalyser.py -l paths_to_pdbs -n atom_names --(cyclohexane|cyclopentane|benzene)`

# Electron density coverage analysis scripts

## Installation

No installation is required.
After cloning the repository, change directory to `electron_density_coverage_analysis`.


## Requirements

Python version 3.10 or higher
`gemmi` package versions 0.6.5, can be installed using pip, e.g.:
```
pip install gemmi==0.6.5
```

## Running analysis for a single structure
Use `electron_density_coverage_analysis.py`, which is a program that determines electron density coverage of a cycle. By default, uses trilinear interpolation to infer the electron density value, and considers an atom to be covered by the electron density when the corresponding intensity is MORE than the threshold for the isosurface (1.5 sigma).

### Arguments description

#### Positional arguments (required)
  - `input_cycle_pdb`      Input PDB file with coordinates of a cycle, produced by `PatternQuery`
  - `input_density_ccp4`   Input electron density file for the corresponding protein structure from the PDB in CCP4 format      

#### Mandatory arguments
Choose one of the two modes:
  - `-s`                   Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle
  - `-d`                   Detailed mode - output a sequence of paried values, where the first member of each pair is a serial atom number, and the second member is "y" - if atom is covered, and "n" - if it is not.

#### Optional arguments
  - `-m`, `--more_or_equal`  Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface       
  - `-c`, `--closest_voxel`  Instead of trilinear interpolation, the intensity of the closest voxel is used


### Example of running analysis for a single structure

On Windows, we recommend to use `Anaconda PowerShell Prompt`

From `electron_density_coverage_analysis` directory run the script in e.g. `simple` mode with e.g. `--more_or_equal` setup:
```
python electron_density_coverage_analysis.py -sm tests/example_input/single_structure/3biu_0.pdb tests/example_input/single_structure/3biu.ccp4
```



## Running analysis for multiple structures
It is possible to run the analysis for multiple structures using `main.py` script.
By default, analysis uses trilinear interpolation to infer the electron density value, and considers an atom to be covered by the electron density when the corresponding intensity is MORE than the threshold for the isosurface (1.5 sigma).

Input should follow a specific format:
![Alt text](electron_density_coverage_analysis/image.png) 

The names of directories should be always as on the picture (`validation_data` as a main folder, then folders for three cycle types, folder `filtered_ligands` in each one, then folders for each ligand, then `patterns` folder in each of them), which contains `.pdb` files of cycles produced by `Pattern Query`. Filename format should always be as in the example: `{ligandID}_{pdbID}_{i}`, ligandID is ligand ID can be up to 3 characters long.

### Arguments description

#### Positional arguments (required)
  - `rootdir`      Root directory (e.g. "validation_data") with files containing coordinates of cycles, produced by `PatternQuery`. Directory hierarchy should follow a specific format (see above)

#### Optional arguments
  - `-m`, `--more_or_equal`  Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface       
  - `-c`, `--closest_voxel`  Instead of trilinear interpolation, the intensity of the closest voxel is used



### Example of running analysis for a multiple structures

On Windows, we recommend to use `Anaconda PowerShell Prompt`

Create `./ccp4` folder in `electron_density_coverage_analysis` folder:
```
mkdir ccp4
```

Download `.ccp4` files for each of the ligands in input folder and put them into `./ccp4` folder. E.g. for `tests/example_input/validation_data` input, you should download the following `.ccp4` files:
 - [4omc](https://www.ebi.ac.uk/pdbe/coordinates/files/4omc.ccp4)
 - [3v8d](https://www.ebi.ac.uk/pdbe/coordinates/files/3v8d.ccp4)
 - [8gdi](https://www.ebi.ac.uk/pdbe/coordinates/files/8gdi.ccp4)
 - [8gxh](https://www.ebi.ac.uk/pdbe/coordinates/files/8gxh.ccp4)
 - [8gxi](https://www.ebi.ac.uk/pdbe/coordinates/files/8gxi.ccp4)
 - [8gxp](https://www.ebi.ac.uk/pdbe/coordinates/files/8gxp.ccp4)


From `electron_density_coverage_analysis` directory run e.g.:

```
python main.py tests/example_input/validation_data -m
```
This will run the analysis in `--more_or_equal` mode.

Output of the script will be collected in `./output` folder including:
 - Three `.csv` files for each ligand type, e.g. `cyclopentane_params_m_analysis_output.csv`, in names of which parameters (`m` in this case) for analysis are specified.
 Each `.csv` file follow a specific format, e.g.:

  06L_8gxp_0,06L,2;5\
  0GV_3v8d_1,0GV,5;5\
  0GV_8gdi_0,0GV,5;5\
  0GV_3v8d_0,0GV,5;5

  where first column specifies ligand PDB file name without extension, second column specifies ligand, third column specifies coverage - first number is a number of covered atoms, then semicolon, then the total number of atoms in a cycle





