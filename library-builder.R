library(devtools)
library(roxygen2)

pkg.name <- 'GI'

desc <- list('Title' = 'Generation intervals',
             'Description'='Correctly handles generation interval distributions for epidemic models.',
             'License' = 'MIT',
             'Version' = '0.1',
             'Author' = 'David Champredon',
             "Maintainer" ="'David Champredon' <david.champredon@gmail.com>")

# Basic folder structure and files setup
devtools::create(path = pkg.name,
                 description = desc)

# Add source files:
rfiles <- c('GI.R', 'RESuDe_FCT.R', 'SEmInR_det.R', 'ct.R')
for(i in seq_along(rfiles))
    system(paste0('cp ',rfiles[i],' ',pkg.name,'/R/'),
           intern = TRUE)

# Documentation
devtools::document(pkg = pkg.name, clean = TRUE)

# ('Depends' is stronger than 'Imports')

# Import other packages:
pkg.imports <- c('deSolve','bbmle')
for(i in 1:length(pkg.imports))
    devtools::use_package(package = pkg.imports[i],
                          pkg = pkg.name,
                          type = 'Imports')
# Dependencies on other packages:
pkg.dep <- c('snow', 'snowfall')  # <-- need this bc of unsolved issue with 'snow': https://stackoverflow.com/questions/36284643/importing-snowfall-into-custom-r-package
for(i in seq_along(pkg.dep))
    devtools::use_package(package = pkg.dep[i],
                          pkg = pkg.name,
                          type = 'Depends')
