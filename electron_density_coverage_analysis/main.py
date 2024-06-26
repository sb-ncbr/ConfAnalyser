import csv
import logging
from multiprocessing import Pool, cpu_count
import os
import urllib.request as r
import argparse
from pathlib import Path
import shutil

# TODO: make it a module?
from electron_density_coverage_analysis import run_as_function

CCP4_DIR = Path('./ccp4')
EXE = Path('./electron_density_coverage_analysis.py')

OUTPUT_DIR = Path('./output')
# NO_CCP4_AVAILABLE_FILE = Path(OUTPUT_DIR / 'no_ccp4_pdb_ids.txt')
CPU_COUNT = cpu_count() / 2

def _create_output_folder():
    try:
        if OUTPUT_DIR.exists():
            shutil.rmtree(str(OUTPUT_DIR.resolve()))

        OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    
    return OUTPUT_DIR

def _download_single_file(data: tuple):
    link, path, pdb_id = data
    try:
        r.urlretrieve(link, path)
    except:
        return pdb_id

# def download_ccp4_multithread(args: argparse.Namespace):
#     try:
#         _create_output_folder()
        
#         input_dir = Path(args.rootdir)
#         pdbe_base_link = 'https://www.ebi.ac.uk/pdbe/coordinates/files/'
#         no_ccp4_available_file = str(NO_CCP4_AVAILABLE_FILE.resolve())

#         l = []

#         for dirname, dirs, files in os.walk(str(input_dir.resolve())):
#             for f in files:
#                 if f.endswith(".pdb"):
#                     pdb_id = f.split("_")[1]
#                     link = pdbe_base_link + pdb_id + ".ccp4"
#                     path = str(Path(CCP4_DIR / f'{pdb_id}.ccp4'))

#                     l.append((link, path, pdb_id))

#         # removing duplicates
#         l = set(l)
#         l = list(l)

#         pool = Pool(int(CPU_COUNT))
#         results = pool.map(_download_single_file, l)

#         results_without_none = [x for x in results if x is not None]

#         pool.close()
#         pool.join()

#         with open(no_ccp4_available_file, 'w') as w:
#             for i in results_without_none:
#                 w.write(i + '\n')
#                 # print(i)

#     except Exception as e:
#         logging.error(e, stack_info=True, exc_info=True)

def process_args(args: argparse.Namespace):
    try:
        a = argparse.Namespace()
        a.s = True
        a.d = False
        a.closest_voxel = False
        a.more_or_equal = False
        if args.closest_voxel:
            a.closest_voxel = True
        if args.more_or_equal:
            a.more_or_equal = True
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

    return a

def run_exe(ligand_filepath: Path, arguments: argparse.Namespace):
    try:
        PQ_pdb_name = ligand_filepath.name.split(".")[0]
        pdb_id = PQ_pdb_name.split('_')[1]
        residue_id = ligand_filepath.parent.parent.name
        ccp4_filepath = str(Path(CCP4_DIR / f'{pdb_id}.ccp4').resolve())
        
        arguments.input_cycle_pdb = str(ligand_filepath.resolve())
        arguments.input_density_ccp4 = ccp4_filepath
        output = run_as_function(arguments)
        result = (PQ_pdb_name, residue_id, output)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return result

def get_filepathes(rootdir: Path, ring_type: str):
    try:
        l = []

        all_ccp4_files = CCP4_DIR.glob('**/*')
        pdb_ids_for_which_ccp4_is_available = [x.stem for x in all_ccp4_files]
        for f in Path(rootdir / ring_type / 'filtered_ligands').rglob("*"):
            if f.is_file():
                # get path
                stem = f.stem
                pdb_id = stem.split('_')[1]
                if pdb_id in pdb_ids_for_which_ccp4_is_available:
                    l.append(f)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)
    return l

def run_analysis(args: argparse.Namespace):
    try:
        arguments = process_args(args)
        params = ''
        if arguments.closest_voxel:
            params = params + "c"
        if arguments.more_or_equal:
            params = params + "m"

        ring_types = ['cyclohexane', 'cyclopentane', 'benzene']
        for ring_type in ring_types:
            cvs_filename = ring_type + '_params_' + params + '_analysis_output.csv'
            with open(str(Path(OUTPUT_DIR / cvs_filename).resolve()), mode='w', newline='') as f:
                rows = []
                w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                filepathes = get_filepathes(Path(args.rootdir), ring_type)
                # print(filepathes)
                with Pool(int(CPU_COUNT)) as p:
                    modified_filepathes = [(f, arguments) for f in filepathes]
                    rows = p.starmap(run_exe, modified_filepathes)
                    w.writerows(rows)
    except Exception as e:
        logging.error(e, stack_info=True, exc_info=True)

def main():
    parser = argparse.ArgumentParser(description='ED coverage analysis')
    parser.add_argument('rootdir', type=str, help='Root directory (e.g. "validation_data")')
    
    parser.add_argument('-s', action='store_true', help='Simple mode - output is two numbers: first is the number of covered atoms, the second is the total number of atoms in a cycle')
    parser.add_argument('-m', '--more_or_equal', action='store_true', help='Atom is considered to be covered by the electron density when the corresponding intensity is MORE OR EQUAL to the threshold for the isosurface')
    parser.add_argument('-c', '--closest_voxel', action='store_true', help='Instead of trilinear interpolation, the intensity of the closest voxel is used')
    
    args = parser.parse_args()
    _create_output_folder()
    run_analysis(args)

if __name__ == '__main__':
    main()
