# C4-Homology
This project computes the  RO(C<sub>4</sub>) homology of a point as a Green Functor, and beyond...
 
Read the wiki for the full documentation/whitepaper. This is enough just to get you started

## Requirements
To run the code you need a sufficiently recent version of MATLAB.
 I have only tested it with v. R2019a but you can try your luck if you have previous versions of MATLAB or freeware like Octave. 
 If you have a University account, you might be able to get an academic MATLAB license for free through your institution.
 At some point the project might be ported to C++ for wider accessibility and speed.

## How do I compute the C4 homology of a point right now?
Type the following in the command line:

```
write_Data(8,7,1);
load_Data;
test_Pure_Homology(8,7,1,Data);
```

This will print the homology of <img src="http://latex.codecogs.com/svg.latex?S^{n\sigma+m\lambda}" border="0"/> for n between 0 and 8 and m between 0 and 7 in the form 
```
The k homology of the (n,m) sphere is MackeyFunctorSymbol
```
where MackeyFunctorSymbol is our notation of the corresponding Mackey functor.
Since the MATLAB display output does not support Latex, MackeyFunctorSymbol is the Latex *code* of the symbol from our paper. For example, "overline Z/2" stands for <img src="http://latex.codecogs.com/svg.latex?\overline{\langle \mathbb{Z}/2\rangle }" border="0"/>

The non obvious differences are as follows:


## The repository is organized into 4 folders:

### General: 
This contains general purpose functions such us the Smith Normal Form, the homology of a complex and the box product of equivariant chains. It also contains transfer and restriction functions that work with C_{2^n} that can be easily modified to work with  C_{p^n} for p a prime number. Eventually we expect them to get general enough to include all cyclic groups, but for now it's just C_{2^n}.

### C4 Specific: 
This contains the algorithms that (together with the general purpose ones) compute the C_4 homology of a point. They can serve as a template for the C_{2^n} case but more changes would be needed.

### C4 Data: 
When using the functions in General or C4 specific, there are a few calculations that are constantly repeated. This folder contains functions that can perform these calculations once and save the results in .m files. These results are then loaded into a single variable ```Data``` of type struct. This is input for many of the functions in the other two folders. ```Data``` also includes variables that are fixed, such as a Mackey functor list that includes the names of the Mackey functors corresponding to their Lewis diagrams. 

### Test C4:
This contains functions that test the output of the C4 specific algorithms against the results in our paper. They can be used to either test those results (if one trusts the algorithms), or the algorithms (if one trusts the math). Ultimately, the best use of it is regression testing: If a change is made in the core algorithms, it's best to run these tests to make sure it's not a regression. 
