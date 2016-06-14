# phylomice
Extension of R [mice](https://github.com/stefvanbuuren/mice) package adding new imputation methods that incorporate phylogenetic information.

You should be familiar with the mice package to use this one.

## Installation
You can install the package using [devtools](https://github.com/hadley/devtools):
```
devtools::install_github("pdrhlik/phylomice")
```

## Available methods
This is a list of methods that are available in the package. Each of these methods requires two additional arguments to the mice function: `psi` and `psiinv`.

* `psi` is a covariance matrix created from the phylogenetic tree,
* `psiinv` is the inverse of `psi`.

Both of these arguments need to be computed before running mice because of the computational cost.

### Creating `psi` and `psiinv`
Let's assume that you have your tree loaded in a `tree` variable. You can then create `psi` and `psiinv` using the `precomputePsi` helper function. You need to have the [ape](https://github.com/cran/ape) package installed.
```
prec <- precomputePsi(tree)

str(prec)
List of 2
 $ psi   : ...
 $ psiinv: ...
```

### phnorm
Imputes univariate continuous missing data using the generalized least square approach. It is based on the **mice.impute.norm** method in **mice** package.
```
library(phylomice)
imp <- mice(data, method = 'phnorm', psi = psi, psiinv = psiinv)
```

If you have precomputed `psi` and `psiinv` using `precomputePsi`, your call would look like this:
```
imp <- mice(data, method = 'phnorm', psi = prec$psi, psiinv = prec$psiinv)
```
### phpmm
This method is currently under development.

## Authors

* **Patrik Drhlik** - *Initial work* - [pdrhlik](https://github.com/pdrhlik)
* **Simon P. Blomberg** - *Theoretical background*
 
If you are interested in improving current or creating new imputation methods that use phylogenies, don't hesitate to contribute.

## License
This project is licensed under the GPL-3 License.

## Acknowledgments
This project wouldn't be possible without the following parties:
* The University of Queensland, Australia
* Technical University of Liberec, Czech Republic
* NESSIE Erasmus Mundus project
