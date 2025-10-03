# notameStats - Non-targeted metabolomics feature prioritization

The notameStats package is the notame package that comes in handy when 
identifying interesting features in non-targeted metabolomics datasets. 
Provides streamlined univariate and multivariate functionality. Check out the 
[website](https://hanhineva-lab.github.io/notame/) for more on notame.

## Installation and getting started

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("notameStats")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("notameStats")
```

If you find any bugs or other things to fix, please submit an issue on GitHub! 
All contributions to notame are always welcome!