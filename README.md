# flooddrake
Firedrake DG shallow water model with wetting and drying for flood simulation

## Linux terminal installation
As a prerequisite, it is required to have Firedrake installed alongside it's dependencies. Instructions on how to this can be found at: http://firedrakeproject.org/download.html. Then, it is required to have the Firedrake virtualenv activated; instructions for this can be found by again following the aforementioned link.

To install, type the following commands into the terminal:

1. `git clone https://github.com/firedrakeproject/flooddrake`
2. `pip install -e ./flooddrake`

## Generating the documentation
To generate the documentation, `flooddrake_doc.pdf`, type the following commands into the terminal inside the `flooddrake` repository:

1. `cd ./docs`
2. `make docs`

## Contact
For any enquiries, please contact: a.gregory14@imperial.ac.uk
