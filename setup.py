#!/usr/bin/env python3
from setuptools import setup, find_namespace_packages



setup(
    name='GCparagon',
    version='0.6.15',
    description='commandline tool to correct GC bias in cfDNA WGS data',
    author='Benjamin Spiegl',
    author_email='benjamin.spiegl@medunigraz.at',
    url='https://github.com/BGSpiegl/GCparagon',
    packages=find_namespace_packages('src', exclude=('test*.py',)),
    package_dir={"": "src"},
    package_data={
        'GCparagon': ['2bit_reference/*.2bit',
                      '2bit_reference/*.chrom.sizes'],
        'GCparagon-AccessoryFiles': ['*_reference_GC_content_distribution.tsv',
                                     '*_minimalExclusionListOverlap_*_intervals_*OverlapLimited.FGCD.bed']},
    entry_points={
        'console_scripts': [
            'gcparagon=GCparagon.correct_GC_bias:main'  # terminal_command_name=python_script_name:main_method_name
        ]}
)

# simply install with "pip install ."
