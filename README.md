This is a sparse collection of programs that implement a case-by-case solution to vary of different problems, using libraries available in the QcmPlab Organization. 

The programs are divided in different directories, named after the method used: **Dmft_ED** for Exact Diagonalization based DMFT method, **SlaveSpins** for the Slave Spins method, **Lattice_ED** for programs implementing lattice Exact Diagonalization, etc. 
 
### Build instructions

- Install in your system the required the QcmPlab libraries that you need (`SciFor`,`DMFT-tools`, etc.)
- Edit to your needs the upper part of the `makefile` in the suitable working directory where your code reside.  
This amounts to set variables for:
    - The driver name _without_ `.f90` extension, e.g. `EXE=ed_bhz_2d`    
    - The fortran compiler [`FC`] of choice: `ifort`, `gnu` or `mpif90` alias are supported.   
    - The platform of choice; either `gnu` or `intel`, to define compilation options.   
    - The target directory for the compiled executable [`DIREXE`]; remember to add it to the `PATH`[^*]  
- Open a terminal therein and run `$ make`: the executable will be built and placed in `DIREXE`
- Open in your working directory a terminal and run `<driver-name>`. 
    - A default input file will be dumped therein as `used.input<MODELNAME>.<extension>`
- Edit to your needs the input file with your text editor of choice and save it as `input<MODELNAME>.<extension>`
- Then you can run again the executable and start crunching numbers!

[^*]: Append ```export PATH=<DIREXE>:$PATH``` to your ```$HOME/.bashrc``` file
