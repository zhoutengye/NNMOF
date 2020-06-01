

# NNMOF

Nueral Network boosted Moment of Fluid


## Note:

(1) Execute with command

    mpirun -np (X) (executable) (input)

-   (X) should be consistant with px\*py in input file

-   (executable) is generated by Makefile

-   (input) is the name of the input namelist file, without suffix .namelis
    for example, for file test3.namelist, the (input) should be test3.
    If no (input) is given, it will search for the default input file input.namelist.

(2) Only HDF5 is supported for input and output field files.

(3) Only support one ghost cell layer.
    For MPI exchange, always 0 and n(i+1) as halo nodes.
    For centered variable, 0 and n(i+1) are boundary nodes.
    For face variable, 0 and n(i) are boundary nodes.

(4) Now boundary conditoin only support fixed value Dilichilet and Nuemann.
The `Field` type contains the BC information. `lohi` indicates the index of the boundary in each 
direction, `bound_type` is the type of boundary conditoin, 1 means Dilichlet, 2 means Nuemann;
`bound_value` is the value of the boundary conditoin. The shape of `lohi`, `bound_value` and 
`bound_type` are `(3,2)`. 
When calling


## GTD


### Basic VOF<code>[0/4]</code>

-   [ ] Basic MPI from CaNS
-   [ ] VOF-PLIC, APPLIC, THINC, THINCSW
-   [ ] Advection of centroid
-   [ ] MOF


### MOF <code>[0/3]</code>

-   [ ] Fkeas
-   [ ] MOF forward data generation
-   [ ] Parameters optimization


### Numerical tests<code>[0/1]</code>

-   [ ] Analytic slope
-   [ ] Resoncnstruction
-   [ ] Zalesak
-   [ ] Single vortex


### Write paper<code>[0/3]</code>

-   [ ] Abstract and introduction
-   [ ] Describe method
-   [ ] Numerical tests

