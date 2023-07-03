# Analysing mechanisms of DI RNA de novo generation and predicting the competitiveness of new DI RNA candidates

This repository includes the code for a master thesis about Defective Interfering Particles written at the Institute for Computational Systems Biology at the University of Hamburg.


## Installation and setup

To setup the repository for usage do the following steps:


1. Move into the desired folder and clone the repository using

```
git clone https://github.com/LohmannJens/MA_DIPs.git
```

2. Move into root folder and run the setup script

```
cd MA_DIPs
bash setup.sh
```

3. Download datasets [here](https://drive.google.com/drive/folders/1QAxqjZMCb7OJyK3GxBGogCJKs3rX65Ya) and unzip them

4. Move downloaded datasets into data folder

The folder structure should look like this afterwards:

<pre>
├── src
|    └── [...]
├── data
|    └── [...]
├── results
├── .gitignore
├── README.md
├── setup.sh
└── env.yml
</pre>


### python packages and conda environment

To run the scripts python version 3.9.12 including the following packages:

| Library         | Version | Website                               |
|-----------------|---------|---------------------------------------|
| Biopython       | 1.78    | www.biopython.org                     |
| Matplotlib      | 3.5.1   | www.matplotlib.org                    |
| Matplotlib-venn | 0.11.7  | www.pypi.org/project/matplotlib-venn/ |
| NumPy           | 1.21.5  | www.numpy.org                         |
| pandas          | 1.4.2   | www.pandas.pydata.org                 |
| scikit-learn    | 1.0.2   | www.scikit-learn.org/stable/          |
| SciPy           | 1.7.3   | www.scipy.org                         |
| ViennaRNA       | 2.5.0a5 | www.pypi.org/project/ViennaRNA/       |

In the root of the repository a conda environment called env.yml is given. It includes all necessary libraries and python packages. It can be installed by the following typing the following comand:

```bash
conda env create -n <new env name> --file env.yml
```

## Help
For questions and suggestions open a new ticket in this repository
