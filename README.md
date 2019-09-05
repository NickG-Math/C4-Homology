# C4-Homology
This project computes the RO(C4) homology of a point as a Green Functor, and beyond...
 
Read the wiki for the full documentation/whitepaper. This is enough just to get you started

## Requirements
To run the code you need a sufficiently recent version of MATLAB.
 I have only tested it with v. R2019a but you can try your luck if you have previous versions of MATLAB or freeware like Octave. 
 If you have a University account, you might be able to get an academic MATLAB license for free through your institution.
 At some point the project might be ported to C++ for wider accessibility and speed.

## How do I compute the C4 homology of a point right now?
Open the test algorithms folder and any of the three test scripts.  



## The repository is organized into 4 folders:

### General: 
This contains general purpose functions such us the Smith Normal Form, the homology of a complex and the box product of equivariant chains. It also contains transfer and restriction functions that work with C_{2^n} that can be easily modified to work with  C_{p^n} for p a prime number. Eventually we expect them to get general enough to include all cyclic groups, but for now it's just C_{2^n}.

### C4 Specific: 
This contains the algorithms that (together with the general purpose ones) compute the C_4 homology of a point. They can serve as a template for the C_{2^n} case but more changes would be needed.

### C4 Data: 
When using the functions in General or C4 specific, there are a few calculations that are constantly repeated. This folder contains functions that can perform these calculations once and save the results in .m files. These results are then loaded into a single variable ```Data``` of type struct. This is input for many of the functions in the other two folders. ```Data``` also includes variables that are fixed, such as a Mackey functor list that includes the names of the Mackey functors corresponding to their Lewis diagrams. 

### Test Algorithms:
This contains scripts that test the output of the C4 specific algorithms against the results in our paper. They can be used to either test those results (if one trusts the algorithms), or the algorithms (if one trusts the math). Ultimately, the best use of it is regression testing: If a change is made in the core algorithms, it's best to run these tests to make sure it's not a regression. 
