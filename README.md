# runibic
Unibic Biclustering algorithm for R


# Compilation

After each modification in .cpp file, the package needs to be recompiled:
```R
> library(Rcpp)
> Rcpp::compileAttributes()
```

After recompiling the package, in order to build the package, type the following command
```
$ export PKG_CXXFLAGS='`Rscript -e "Rcpp:::CxxFlags()"` -std=c++11 -fopenmp'
$ make

```

# Usage

```
library(Rcpp)
a=replicate(10, rnorm(20))
a
sourceCpp("src/runibi.cpp")
b=unisort(a)
b
calculateLCS(b[1,],b[2,])
```
