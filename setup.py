from setuptools import setup, find_packages

setup(
    name="puzzler",
    version="1.9",
    description="Puzzler pipeline with HapHiC integration",
    packages=find_packages(include=["HapHiC", "HapHiC.*"]),
    # CLI / binaries that should end up in $PREFIX/bin
    scripts=[
        "bin/puzzler",
        "HapHiC/scripts/haphic",
        "HapHiC/utils/filter_bam",
        "HapHiC/utils/agp_to_fasta",
    ],
    # non-Python data that should be installed into site-packages/HapHiC/…
    package_data={
        "HapHiC": [
            "scripts/*.py",
            "utils/*",
        ]
    },
    include_package_data=True,
    python_requires=">=3.12,<3.13",
)
