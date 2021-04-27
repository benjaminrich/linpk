# Version 1.1.1

## Test environments

* Local:
  - NixOS (Linux), R 4.0.4 (x86_64-pc-linux-gnu)
  - Windows 10, R 4.0.5 (x86_64-w64-mingw32/x64 (64-bit))
* travis-ci:
  - Ubuntu Linux 16.04.6 LTS (xenial) (release and devel)
* win-builder:
  - Windows Server 2008 (release (4.0.4) and devel (unstable))
* R-hub builder (https://builder.r-hub.io)
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
    - This had the following problem, which is unrelated to my package:
      ```
      Error: Bioconductor version '3.13' requires R version '4.1'; R version is too new; see
      https://bioconductor.org/install
      ```

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.
