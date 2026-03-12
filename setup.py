# resscan/setup.py

import setuptools
import re

def get_version():
    with open("resscan/resscan.py", "r") as f:
        for line in f:
            match = re.search(r'__version__ = "(.*?)"', line)
            if match:
                return match.group(1)
    raise RuntimeError("Version could not be found.")

setuptools.setup(
    name="resscan",
    version=get_version(),
    author="H. Soon Gweon",
    author_email="h.s.gweon@reading.ac.uk",
    description="A Bioinformatics Tool for AMR gene and variant detection from metagenomic data.",
    packages=setuptools.find_packages(),
    
    install_requires=[
        'pandas',
        'numpy',
    ],

    include_package_data=True,
    package_data={
        'resscan': ['databases/resscan_DB_SCG/*']
    },

    scripts=[
        'scripts/resscan_DB_USCG_download.py',
        'scripts/resscan_batch',
        'scripts/resscan_aggregate',
        'scripts/generate_merged_samplesheet'
    ],

    entry_points={
        'console_scripts': [
            'resscan = resscan.resscan:main',
            'resscan_build_db = resscan.resscan_build_db:main',
            'resscan_curate_metadata = resscan.resscan_curate_metadata:main',
        ],
    },
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
)