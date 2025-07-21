## LSGI: Local Spatial Gradient Inference

### Motivation
Cellular anatomy and signaling vary across niches, which can induce gradated gene expressions in subpopulations of cells. We introduce this framework, the Local Spatial Gradient Inference  (LSGI), to systematically identify spatial locations with prominent, interpretable spatial transcriptomic gradients (STGs) from spatial transcriptomic (ST) data.

### Install

```
# There is a dependency package 'anticlust' and please use older versions of the package (e.g., 0.6.1).

# Development version (recommended):

# install.packages("remotes")

#install from github
remotes::install_github("https://github.com/qingnanl/LSGI")
```

### Tutorial

For a tutorial, please refer to [Using LSGI to analyze 10x Visium data](http://htmlpreview.github.io/?https://github.com/qingnanl/LSGI/blob/master/vignette/LSGI-Tutorial.html)

We added a tutorial for single-gene type analysis [Using LSGI to analyze seqFISH data](http://htmlpreview.github.io/?https://github.com/qingnanl/LSGI/blob/master/vignette/LSGI-demo_seqFISH.html)

For a better explanation of this method, please refer to our manuscript: https://www.biorxiv.org/content/10.1101/2024.03.19.585725v1

### Related links

Here we have some tumor ST data processed with LSGI: https://zenodo.org/records/10626940

Please also refer to the manual on this Zenodo repo to explore the data (some visualization functions are tailored for highlighting tumor regions).
