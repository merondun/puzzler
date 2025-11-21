from setuptools import setup, find_packages

scripts = [
    "bin/puzzler",
    "HapHiC/utils/filter_bam",  # Rust binary on PATH
]

setup(
    name="puzzler",
    version="1.9",
    packages=find_packages(include=["HapHiC", "HapHiC.*"]),
    package_data={
        "HapHiC": [
            "scripts/*.py",
            "utils/*.py",
            "utils/juicer_tools.1.9.9_jcuda.0.8.jar",
            "utils/agp_to_fasta",
            "utils/juicer",
        ]
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "haphic = HapHiC.haphic:main",
        ]
    },
    python_requires="==3.12.4",
)