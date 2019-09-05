# C4-Homology
This project computes the  RO(C<sub>4</sub>) homology of a point as a Green Functor.
 
Read the wiki for the full documentation if you wish to contribute or check the source code. What follows is an FAQ that should be enough to verify the results of our paper in your preferred range.

## What are the requirements to run the program?
To run the code you need a version of MATLAB.
 I have only tested it with v. R2019a but you can try your luck if you have previous versions of MATLAB or freeware like Octave. 
 If you have a University account, you might be able to get an academic MATLAB license for free through your institution.
 At some point the project might be ported to C++ for wider accessibility and speed.

## How do I compute the C4 homology of a point with your program?
Download the repository and add it to your MATLAB path. Then type the following in the command line:

```
write_Data(1,1,1);
load_Data;
test_Pure_Homology(8,7,0,Data);
```

This will compute the homology of 
<img src="http://latex.codecogs.com/svg.latex?S^{n\sigma+m\lambda}" border="0"/> for n=0,...,8 and m=0,...,7, check the answer against the tables in our paper and print the answer in the form 
```
The k homology of the (n,m) sphere is MackeyFunctorSymbol
```
where MackeyFunctorSymbol is our notation of the corresponding Mackey functor.

Since the MATLAB display output does not support Latex, MackeyFunctorSymbol is the Latex *code* of the symbol from our paper. For example, "overline Z/2" stands for <img src="http://latex.codecogs.com/svg.latex?\overline{\langle \mathbb{Z}/2\rangle }" border="0"/>

If you want to check a different range n=0,...,rangeN, m=0,...,rangeM run
```
test_Pure_Homology(rangeN,rangeM,0,Data);
```

If you want to check the homology of <img src="http://latex.codecogs.com/svg.latex?S^{-n\sigma-m\lambda}" border="0"/> run
```
test_Pure_Cohomology(rangeN,rangeM,0,Data);
```
To check <img src="http://latex.codecogs.com/svg.latex?S^{m\lambda-n\sigma}" border="0"/> run
```
test_Sigma_Minus_Lambda(rangeN,rangeM,0,Data);
```
Finally to check <img src="http://latex.codecogs.com/svg.latex?S^{n\sigma-m\lambda}" border="0"/> run
```
test_Lambda_Minus_Sigma(rangeN,rangeM,0,Data);
```

## But what about the multiplicative structure?

To be added once it's user friendly enough


## The program runs too slowly for large ranges

Significant improvements can be had by using precomputed Data so as to avoid repeat calculations. Once you have decided on the n,m ranges, run the following commands:

```
write_Data(rangeN,rangeM,1);
load_Data;
```

After that you can run the test functions with the third input being 1 eg

```
test_Pure_Homology(rangeN,rangeM,1,Data);
```

If you are still not thrilled about runtime speed, you can use the parallel processing package (see the Performance section in the sidebar of the wiki).


## For more details see the wiki
