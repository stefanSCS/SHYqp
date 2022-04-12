

#  Bezier5YS and SHYqp

This code implements the Bezier5YS and SHYqp models of the yielding and flow properties of metallic alloys. Applications range over the entire spectrum of lattice structures: e.g., magnesium, aluminum, iron. The only symmetry assumed is that of orthotropic material symmetry (featured in general by metal sheet).
The underlying theory and algorithms are presented in:

- Preprint: *Bezier5YS and SHYqp: A general framework for generating data and for modeling symmetric and asymmetric orthotropic yield surfaces* ()

- .


## USAGE

The code is written in python. It requires the following packages: [numpy](https://numpy.org/install/), [matplotlib](https://matplotlib.org/stable/users/installing.html), [cvxopt](https://cvxopt.org/install/index.html) and  [quadprog](https://pypi.org/project/quadprog/). Testing has been done on Windows only (but it should work on Linux also), so the instructions are specific to this OS (but they pertain to Linux also, modulo any reference to '.exe'). The simplest work-flow is:

- Install [Python](https://www.python.org/) (version >=3.7)
- Make sure python.exe is on the search path (i.e., it can be run from a command prompt); If not, add it to the 'PATH' environment variable
- Install the required packages via 'pip'; In a command prompt run:
    - pip install numpy, matplotlib, cvxopt, quadprog

Note, however, that this mixing of specific libraries (such as cvxopt and quadprog) onto the global python installation is not the recommended approach. The best work-flow is via a python virtual environment where libraries and specific applications contexts are kept in isolation from the main installation. Thus after installing Python, create a directory where work on material modeling is to take place. From now on I'll assume this folder is named 'MATM'. Open PowerShell, navigate to 'MATM' and execute:

-  ` ` `python -m venv OPTIM ` ` `
    - (the above command will create a subdirectory 'OPTIM' - the virtual environment - with a bare bones pseudo-installation of Python)
- navigate to 'OPTIM' and execute:
    - ` ` `.\Scripts\Activate.ps1` ` `
        - (the above activates the virtual environment 'OPTIM')
    - ` ` `pip install numpy` ` `
    - ` ` `pip install matplotlib` ` `
    - ` ` `pip install cvxopt` ` `
    - ` ` `pip install quadprog` ` `
        - (the last four commands install the required packages)
    - a check of the list of installed packages can be done with the command:
        - ` ` `pip-list` ` `                   
    - exit a virtual environment by executing at the prompt:
        - ` ` `deactivate` ` `

To work with SHYqp, copy the files 'SHYqpV1.py' and 'SHYqp_main.py' into the folder 'OPTIM'. The first file is the actual collection of numerical functions, while the second is the driver containing the basic calls required to calculate a Bezier5YS+SHYqp model. To calculate a model, simply execute (within the activated 'OPTIM' environment):

` ` `python SHYqp_main.py` ` `


## INPUT

The script requires the input file 'mat000File.txt'. Here are recorded the material data and other meta-parameters. This file must be at the same location as the script (i.e., in folder 'OPTIM'). Instructions for creating the input file are detailed in the demo file 'mat001File_Instructions.txt'. Sample input files for all the examples illustrated in the cited article are provided as 'mat\*.txt' files, e.g., 'matAZ31_Lou2007.txt', 'matDP980_Li2020.txt', 'matAA2090T3.txt', etc. Simply copy the content of any of these files and paste it over the content of 'mat000File.txt' to obtain the corresponding Bezier5YS and SHYqp models.


## OUTPUT

Within the 'OPTIM' working directory create a subfolder named 'FIGS'. This is the location where all output is saved (Note: The script aborts execution if it does not detect the 'FIGS' directory). If all goes well and a solution is found, the script will generate the following output:

- a report file '\*_Err_and_Coeff.txt' containing convexity and performance measures, predictions vs actual data, and the material parameters (coefficients) of the SHYqp model;
- plots of Bezier5YS model (directional properties and yield surface);
- plots of SHYqp model (directional properties, yield surface, biaxial sections).



## License

This code is released under the MIT license.
