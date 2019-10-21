
# CRAUT

`CRAUT` is a tool used to generate and rank linkers for proteins using using a rod and angle model.

Our software is designed as a simplified version of the code provided by the 2014 Heidelberg iGEM team. [Here is the their wiki page detailing their version of the software](http://2014.igem.org/Team:Heidelberg/Software/Linker_Software), and [here is their code](https://github.com/igemsoftware/Heidelberg_2014/tree/master/CRAUT).

The linker is constructed using predetermined alpha helix and angle sequences as building blocks. Our version of the software only considers using 2 alpha helices and 1 angle segment. Linkers are ranked using a linear regression model made by the Heidelberg team, described [here](http://2014.igem.org/Team:Heidelberg/Modeling/Linker_Modeling).

## Build Instructions
The software was developed using C++11, and uses the linear algebra libary [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). Version 3.3.7 of this library was used in development. If you wish to compile our code, then you will need to instruct the compiler on where to find the appropriate headers for Eigen. [Here is documentation on compiling Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html). The software was developed using Microsoft Visual Studio, but should still work for other compilers.

All other code is included in this repo.

## Usage
The program will output the linkers sorted from best to worst (smallest weight to largest weight), as well as information on their weighting components. To run the program, a pdb file must be passed as the first argument, and there is also an option of providing an integer for how many linkers are to be printed in the second argument. If no integer is given for the second argument or a negative integer is given, then all linkers will be printed. Passing in invalid input will result in undefined behaviour. The program will output to stdout.

Here is an example of running the software from the command line with the protein file "protein.pdb", which will print the 10 best linkers. The output is redirected to the newly created file "linkers.txt".
```
craut.exe protein.pdb 10 > linkers.txt
```

The output will consist of one line for each linker. Each line will contain the linker sequence, the overall weight, and the individual weight components. Here is an example of one such line:
```
AEAAAKEAAAKEAAAKAKTAAEAAAKEAAAK    weight=3.595392   length_weight=2.704839   angle_weight=5.645027   dist_weight=0.359070
```

If different formatting is desired, it should be easily modifiable in the code.

## Software Design
The main goal of the code was to be simple and efficient, and be easily extendable to include other weighting schemes. To tackle the problem of efficiency, C++ was used, with the math library Eigen for vector math. Most convoluted code translated from the Heidelberg team's code (in attempts to match weighting results) was removed to make the computation simple and easy to understand.

Most of the general functionality is provided in an abstract base class called `Craut`. To use it, another class must be made to extend it. All this derived class needs to do is call the base class's constructor, and implement a weighting function for any given linker. This allows for multiple different weighting schemes to be easily added and used. This was done with the intention of having two schemes, where one closely follows the computations done by the Heidelberg team, and another is a very simple version. The former was removed due to difficulty in understanding and replicating their algorithm.

A simple weighting function that uses their linear regression model was implemented in the class `CrautNew`. The different weight contributions will produce different results, but they are computed with the same idea in mind (i.e. still computing the total length, distance from protein, etc.). However, no binding site contribution is computed because of time constraints, as well as low priority since the binding site is not in the way of linkers for our protein.

To actually generate and rank linkers with this class, little external code in needed. In `main.cpp`, the alpha helix and angle segment building blocks are defined and passed to a newly created instance of `CrautNew`, along with the protein file and the protein subunit to consider. The subunit passed in is `A`, and different subunits of a protein may be observed in column 22 (1-indexed) of a line describing an atom in a `.pdb` file.

Due to time constraints, there is little error handling in the code. Using invalid `.pdb` files will likely result in crashes.

If it is found that having a runtime of a few seconds is too long, then optimizing the method for finding the distance from the linker to the protein will help. This can be done by implementing an axis-aligned bounding box tree for the atoms in the protein, and using that to efficiently find the closest atom to a given point. This should result in a dramatic speedup, since the minimum distance calculation was observed to be the main bottleneck in terms of runtime.
