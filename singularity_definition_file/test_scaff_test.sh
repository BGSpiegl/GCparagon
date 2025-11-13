# TEST gcparagon with scaffold checks:
# (DON'T PULL JUST YET!
#  Instead, test cmdline call using conda env from image build
#  Only if that worked, actually BUILD the image by running: 'sudo singularity build gcparagon_v0.6.15_feature-inputstdchrmcheck.sif gcparagon-testfeature-inputstdchrmcheckBranch.def'
#  DONT THIS -> after pulling sif from sylabs:
#  singularity pull library://bgspiegl/gcparagon/gcparagon_0.6.15:latest && singularity verify gcparagon_v0.6.15_latest.sif
#  run the following analyses:)
#
# WORKSTATION:
# 1) GCparagon - CONDA ENV CALLS (from within the singiularity building dir; after 'conda env create -f /mnt/NVMeScratch/PycharmProjects/GCparagon_public/conda_env/GCparagon_Ubuntu22_py3.10_for_sif_creation.yml --name GCparagon_test_inputstdchrmcheck')
# (should work!) -> WORKED.
python3 ../src/GCparagon/correct_GC_bias.py --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/B01.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_BAMscaffTest  --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24 --unclipped-min-aln-fraction 0.1
# (did not work but should work now!) -> WORKED.
python3 ../src/GCparagon/correct_GC_bias.py --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/Villanueva_DEBUG/CF24_0064.rh.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_BAMscaffTest-Villanueva --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24 --unclipped-min-aln-fraction 0.1
# (should also work now but the short fraction should be gone) -> WORKED.
python3 ../src/GCparagon/correct_GC_bias.py --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/Villanueva_DEBUG/CF24_0064.rh.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_BAMscaffTest-Villanueva --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24
#
# 2) GCparagon - SINGULARITY IMAGE CALLS (after building the image)
# (should work!) -> WORKED.
singularity run -B /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST -B /mnt/NVMeScratch/NVMe_data/tmp /mnt/NVMeScratch/PycharmProjects/GCparagon_public/singularity_definition_file/gcparagon_v0.6.15.sif --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/B01.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_BAMscaffTest  --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24 --unclipped-min-aln-fraction 0.1
# (did not work but should work now!)
singularity run -B /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST -B /mnt/NVMeScratch/NVMe_data/tmp /mnt/NVMeScratch/PycharmProjects/GCparagon_public/singularity_definition_file/gcparagon_v0.6.15.sif --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/Villanueva_DEBUG/CF24_0064.rh.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_singularity-BAMscaffTest-Villanueva  --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24 --unclipped-min-aln-fraction 0.1
# (should also work now but the short fraction should be gone)
singularity run -B /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST -B /mnt/NVMeScratch/NVMe_data/tmp /mnt/NVMeScratch/PycharmProjects/GCparagon_public/singularity_definition_file/gcparagon_v0.6.15.sif --bam /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/input_BAMs/Villanueva_DEBUG/CF24_0064.rh.bam --reference-genome-build hg38 --out-dir /media/benjamin/Analyses/GC_PARAGON_PUBLIC_TEST/GCparagon_output_v0.6.15_singularity-BAMscaffTest-Villanueva  --temporary-directory /mnt/NVMeScratch/NVMe_data/tmp --output-unaligned-reads --preset 1 --threads 24