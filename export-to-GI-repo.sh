### Export to LIVE package Github repository
### Does NOT export the DESCRIPTION, NAMESPACE and LICENSE files

### R files:
cp GI/R/*.R ../GI/R

### C++ source files:
cp GI/src/* ../GI/src

### Documentation:
cp GI/man/* ../GI/man
cp vignettes/*.Rmd ../GI/vignettes

### Data added to the package:
cp data/*.RData ../GI/data
cp data/data.R ../GI/R
