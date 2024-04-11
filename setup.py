#!/usr/bin/env python3
from setuptools import setup, find_namespace_packages

setup(
    name='GCparagon',
    version='0.6.6',
    description='commandline tool to correct GC bias in cfDNA WGS data',
    author='Benjamin Spiegl',
    author_email='benjamin.spiegl@medunigraz.at',
    url='https://github.com/BGSpiegl/GCparagon',
    package_dir={"": "src"},
    packages=find_namespace_packages('src', exclude=('test*.py',)),
    entry_points={
        'console_scripts': [
            'gcparagon=GCparagon.correct_GC_bias:main'  # terminal_command_name=python_script_name:main_method_name
        ]
    }
)

# install with "python setup.py bdist_wheel && python -m pip install --force-reinstall dist/*.whl"
