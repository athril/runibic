# runibic
Unibic Biclustering algorithm for R


# Development

In order to be able to compile the package you need to prepare flags for the compiler:
```
$ export PKG_CXXFLAGS='`Rscript -e "Rcpp:::CxxFlags()"` -std=c++11 -fopenmp'
```

After each modification in .cpp file, the package needs to be recompiled:
```
$ make
$ R
> library(Rcpp)
> Rcpp::compileAttributes()
```

# Compiling

In order to compile runibic, type:
```
make
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
