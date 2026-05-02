from setuptools import setup, find_packages
from pathlib import Path

long_description = (Path(__file__).parent / "README.md").read_text(encoding="utf-8")

setup(
    name             = "gtf-merge",
    version          = "0.1.0",
    author           = "GTF Merge Contributors",
    description      = "High-performance GTF merger for long-read RNA-seq data",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url              = "https://github.com/yourusername/gtf-merge",
    packages         = find_packages(),
    python_requires  = ">=3.8",
    install_requires = [
        "intervaltree>=3.1.0",
    ],
    extras_require = {
        "dev": [
            "pytest>=7.0",
            "pytest-cov",
        ],
    },
    entry_points = {
        "console_scripts": [
            "gtf-merge = gtf_merge.cli:main",
        ],
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    keywords = ["gtf", "rna-seq", "long-read", "pacbio", "iso-seq", "transcriptome"],
)
