from setuptools import setup, find_packages

scripts=[
    "bin/puzzler",          # standalone pipeline
    "HapHiC/scripts/haphic", # HapHiC wrapper
    "HapHiC/scripts/mock_agp_file.py" 
]

setup(
    name="puzzler",
    version="1.9",
    scripts=scripts,
    description="Puzzler pipeline with HapHiC integration",
    packages=find_packages(include=["HapHiC", "HapHiC.*"]),
    # non-Python data that should be installed into site-packages/HapHiC/
    package_data={
        "HapHiC": [
            "scripts/*.py",
            "utils/*",
        ]
    },
    include_package_data=True,
    python_requires=">=3.12,<3.13",
)
