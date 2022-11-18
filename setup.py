import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="delfind",
    version="0.9.10",
    author="Damien Marsic",
    author_email="damien.marsic@aliyun.com",
    description="Find and quantify large deletions in populations of circular genomes from Illumina paired-end data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/damienmarsic/delfind",
    package_dir={'': 'delfind'},
    packages=setuptools.find_packages(),
    py_modules=["delfind"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.6',
)
