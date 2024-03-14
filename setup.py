from setuptools import *

setup(
    name="DupCaller",
    version="0.0.4",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/AlexandrovLab/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=[
        "DupCaller/DupCallerCall.py",
        "DupCaller/DupCallerTrim.py",
        "DupCaller/DupCallerSummarize.py",
    ],
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython>=1.78",
        "pysam>=0.19.0",
        "numpy>=1.21.5",
        "matplotlib>=3.8.0",
        "scipy>=1.11.3",
    ],
)
