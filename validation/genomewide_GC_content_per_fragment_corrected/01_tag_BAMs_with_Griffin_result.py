#!/usr/bin/env python3

import shutil as sh
import subprocess as sp
import multiprocessing as mp

from pathlib import Path

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY AND TEMPORARY DIRECTORY HERE !!!!!!!!
MY_OUTPUT_DIRECTORY = Path('your_absolute_path_to_an_output_directory_with_huge_space')
MY_TEMPORARY_DATA_DIR = Path('your_absolute_path_to_a_temporary_directory_with_huge_space')
# TODO: SET YOUR GCparagon SOURCE DIRECTORY AND TEMPORARY DIRECTORY HERE !!!!!!!!
########################################################################################################################
MY_OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
MY_TEMPORARY_DATA_DIR.mkdir(parents=True, exist_ok=True)

SCRIPT_PARENT_PATH = Path(__file__).parent
SOURCE_CODE_ROOT_PATH = SCRIPT_PARENT_PATH.parent.parent.parent.parent
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

python_path = Path(sh.which('python'))
assert python_path.is_file()
gcparagon_script_path = SOURCE_CODE_ROOT_PATH / 'src/GCparagon/correct_GC_bias.py'
assert gcparagon_script_path.is_file()
twobit_reference_fasta = SOURCE_CODE_ROOT_PATH / 'src/GCparagon/2bit_reference/hg38.analysisSet.2bit'

# order of samples is important! must correspond between bam_path, weight_paths, and masks_paths !!!

input_path = SOURCE_CODE_ROOT_PATH / 'preset_computation'
bam_path = (input_path / 'B01.bam',
            input_path / 'H01.bam',
            input_path / 'C01.bam',
            input_path / 'P01.bam')  # TODO: you must DOWNLOAD these samples from EGA DB first !!!
weights_path = SOURCE_CODE_ROOT_PATH / 'validation/transformed_Griffin_bias_matrices'
weight_paths = (weights_path / 'B01_gc_weights.Griffin.txt.gz',
                weights_path / 'H01_gc_weights.Griffin.txt.gz',
                weights_path / 'C01_gc_weights.Griffin.txt.gz',
                weights_path / 'P01_gc_weights.Griffin.txt.gz')
masks_paths = (weights_path / 'B01_gc_bias_computation_mask.Griffin.txt.gz',
               weights_path / 'H01_gc_bias_computation_mask.Griffin.txt.gz',
               weights_path / 'C01_gc_bias_computation_mask.Griffin.txt.gz',
               weights_path / 'P01_gc_bias_computation_mask.Griffin.txt.gz')


def tag_bam_file(bam: Path, weights_matrix: Path, mask: Path):
    tagging_proc = sp.run([str(python_path), str(gcparagon_script_path),
                           '--bam', str(bam), '--two-bit-reference-genome', str(twobit_reference_fasta),
                           '--tag-only', '--correction-weights', str(weights_matrix), '--weights-mask', str(mask),
                           '--threads', '12', '--tag-name', 'GG', '--default-weight', '0.0',
                           '--out-dir', str(MY_OUTPUT_DIRECTORY), '--temporary-directory', str(MY_TEMPORARY_DATA_DIR)])
    tagging_proc.check_returncode()


if __name__ == '__main__':
    griffin_tagging_processes = []
    for bam_path, weight_matrix_path, mask_path in zip(bam_path, weight_paths, masks_paths):
        griffin_tagging_processes.append(mp.Process(
            target=tag_bam_file, kwargs={'bam': bam_path,
                                         'weights_matrix': weight_matrix_path,
                                         'mask': mask_path}))
    print("starting Griffin tagging processes...")
    _ = [tagging_proc.start() for tagging_proc in griffin_tagging_processes]
    _ = [tagging_proc.join() for tagging_proc in griffin_tagging_processes]
