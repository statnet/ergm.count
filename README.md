# `ergm.count`: Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count Edges

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ergm.count?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/ergm.count)](https://cran.r-project.org/package=ergm.count)
[![Coverage status](https://codecov.io/gh/statnet/ergm.count/branch/master/graph/badge.svg)](https://codecov.io/github/statnet/ergm.count?branch=master)
[![R build status](https://github.com/statnet/ergm.count/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/ergm.count/actions)
[![R-universe](https://statnet.r-universe.dev/ergm.count/badges/version)](https://statnet.r-universe.dev/ergm.count)

A set of extensions for the 'ergm' package to fit weighted networks whose edge weights are counts. See Krivitsky (2012) <doi:10.1214/12-EJS696> and Krivitsky, Hunter, Morris, and Klumb (2023) <doi:10.18637/jss.v105.i06>.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/ergm.count`
* A private repository `statnet/ergm.count-private`

The intention is that all developments in `statnet/ergm.count-private` will eventually make their way into `statnet/ergm.count` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.org/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

[R-Universe](https://r-universe.dev) builds a set of binaries after every commit to the main branch of the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. To obtain the binaries from r-universe, navigate to the package page at https://statnet.r-universe.dev/ergm.count .
You may also want to install the corresponding latest binaries for packages on which `ergm.count` depends, in particular [`network`](https://github.com/statnet/network), [`statnet.common`](https://github.com/statnet/statnet.common), and [`ergm`](https://github.com/statnet/ergm).
