# MachineLearningGWAS: Machine Learning for Genome-Wide Association Studies

## Description

The `MachineLearningGWAS` R package implements machine learning tools for genome-wide association studies.

## Installation

- Install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/).

- Install the `caret` R package with all its dependencies.

```R
install.packages("caret", dependencies = TRUE)
```

- Install [Bioconductor](https://bioconductor.org/).

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

- Install the development version of `MachineLearningGWAS` from Github.

```R
# install.packages("devtools")
devtools::install_github("olivmrtn/MachineLearningGWAS")
```

- Load the package and you're good to go.

```R
library(MachineLearningGWAS)
```

## Usage

Read the vignette [here](https://cdn.rawgit.com/olivmrtn/MachineLearningGWAS/master/inst/doc/MachineLearningGWAS.html).

## Contact

* Olivier M. F. Martin, Pharm.D.
[olivmrtn@gmail.com](mailto:olivmrtn@gmail.com)

## License

MachineLearningGWAS 0.1.0
Copyright (C) 2016

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The GNU General Public License is available at http://www.gnu.org/licenses/

The source code can be found at https://github.com/olivmrtn/MachineLearningGWAS
