# Burgers 1D example

I am working on reproducing the results from
[this](http://rspa.royalsocietypublishing.org/content/470/2171/20140080) paper
with [Proteus code](https://github.com/erdc-cm/proteus). One objective is to learn
Proteus which is DG-FEM code, another is to understand the paper better.

## Structure

The repository is split onto several different directories

  * gpc - contains example generalized polynomial chaos setting for the
    Burgers 1D problegpc - contains example generalized polynomial chaos
    setting for the Burgers 1D problem.
  * graphs - sample code that I used to visualise some of the results
  * high_dimension - example setup for an exponential correlation function.
    The point is here to prove sensitivity of high dimensional IC on the
    forcing term.
  * notebooks - iPython notebooks which serve as documentation
  * riemann - Riemann problem examples

