# Version 1.1.2

## Test environments

* Local:
  - Linux (NixOS), R 4.2.2 (x86_64-pc-linux-gnu)
  - Windows 10, R 4.2.2 (x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
  - Windows Server 2022 (release (4.2.2) and devel (unstable))
* R-hub builder (https://builder.r-hub.io)
  - Windows Server 2022 R-devel, 64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
    - Gave following note, which is unrelated to my package:
      ```
      * checking HTML version of manual ... NOTE
      Skipping checking HTML validation: no command 'tidy' found
      Skipping checking math rendering: package 'V8' unavailable
      ```

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.
