from setuptools import setup

setup(
    name="DupCaller",
    version="0.0.3",
    description="A variant caller for barcoded DNA sequencing",
    url="https://github.com/AlexandrovLab/DupCaller",
    author="Yuhe Cheng",
    author_email="yuc211@ucsd.edu",
    scripts=["DupCaller/DupCallerCall.py", "DupCaller/DupCallerTrim.py","DupCaller/DupCallerLearn.py"],
    install_requires=["biopython>=1.78", "pysam>=0.19.0", "numpy>=1.21.5"],
)
