# C4-Homology
This project computes the  RO(C<sub>4</sub>) homology of a point as a Green Functor.
 
Read the [Wiki](https://github.com/NickG-Math/C4-Homology/wiki) for the full documentation. What follows is an **FAQ** that should be enough to verify the results of our paper in your preferred finite range.

## What do I need to run this code?
All you need is a version of MATLAB.
 I have only tested it with v. R2019a but you can try your luck if you have previous versions of MATLAB or freeware like Octave. 
 If you are associated with a University you might be able to get a free academic MATLAB license through your institution.
 At some point the project might be ported to C++ for wider availability and speed, but for now the code is in MATLAB

## How do I compute the C4 homology of a point?
Download the repository and add it to your MATLAB path. Then type the following in the command line:

```
write_Data(1,1,1,1);
load_Data;
test_Pure_Homology(8,7,0,Data);
```

This will perform three operations:

- Compute the homology of 
<img src="http://latex.codecogs.com/svg.latex?S^{n\sigma+m\lambda}" border="0"/> for n=0,...,8 and m=0,...,7

- check the answer against the tables in our paper, giving out an error if there is a mismatch and

- print the answer in the form 
```
The k homology of the (n,m) sphere is MackeyFunctorSymbol
```
where ```MackeyFunctorSymbol``` is our notation of the corresponding Mackey functor.

Since the MATLAB display output does not support Latex, ```MackeyFunctorSymbol``` is the Latex *code* of the symbol from our paper. For example, "overline Z/2" stands for <SPAN STYLE="text-decoration:overline"><Z/2></SPAN>

If you want to verify a different range, say n=0,...,rangeN, m=0,...,rangeM, run
```
test_Pure_Homology(rangeN,rangeM,0,Data);
```

If you want to check the homology of <img src="http://latex.codecogs.com/svg.latex?S^{-n\sigma-m\lambda}" border="0"/> (in the usual range n=0,...,rangeN and m=0,...,rangeM), run
```
test_Pure_Cohomology(rangeN,rangeM,0,Data);
```
To check <img src="http://latex.codecogs.com/svg.latex?S^{n\sigma-m\lambda}" border="0"/> run
```
test_Sigma_Minus_Lambda(rangeN,rangeM,0,Data);
```
Finally to check <img src="http://latex.codecogs.com/svg.latex?S^{m\lambda-n\sigma}" border="0"/> run
```
test_Lambda_Minus_Sigma(rangeN,rangeM,0,Data);
```

## But what about the multiplicative structure?

For the multiplicative generators of <img src="http://latex.codecogs.com/svg.latex?S^{n\sigma+m\lambda}" border="0"/> run
```
test_Pure_Homology_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
```
The range variables specify the sets that the exponents of the Euler and orientation classes ```asigma, u2sigma,usigma,alambda,ulambda``` are allowed to range in. The relationship with the previous range variables is ```rangeN=rangeN1+2*rangeN2``` and ```rangeM=rangeM1+rangeM2```.

This command does three things:

- Computes the product of these Euler and orientation classes and checks if it's a generator

- checks if the generator is also in our tables, giving out an error if there is a mismatch and

- prints the answer in a variable form depending on the generator, that can look like

```
Top Generator verified: asigma^2*alambda^3
```

or

```
Mid Generator verified: usigma^4*bar(ulambda^2)
```

We don't check bottom level generators since the nonequivariant homology of spheres has obvious multiplicative structure.

To get the other equivariant spheres, run any of the following

```
test_Pure_Cohomology_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
test_Sigma_Minus_Lambda_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
test_Lambda_Minus_Sigma_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
```


## Can I make the program go faster?

Significant speed improvements can be achieved by using precomputed Data so as to avoid repeat calculations. Once you have decided on the n,m ranges   ```rangeN,rangeM```, run the following commands:

```
write_Data(NumberN,NumberM,1,1);
load_Data;
```
where ```NumberN, NumberM``` should be at least  ```rangeN, rangeM``` respectively.

After that you can run the test functions with the third input being 1 eg
```
test_Pure_Homology(rangeN,rangeM,1,Data);
```

Warning: If you run
```
test_Sigma_Minus_Lambda_Mult(rangeN1,rangeN2,rangeM1,rangeM2,1,Data);
test_Lambda_Minus_Sigma_Mult(rangeN1,rangeN2,rangeM1,rangeM2,1,Data);
```
you will get an error 
```Index in position 3 exceeds array bounds (must not exceed 1).```
This is because you didn't precompute enough Data. So you must first run
```
write_Data(NumberN,NumberM,1,LargeNumber);
load_Data;
```
for a large enough ```LargeNumber``` depending on the ranges you are using.

The larger the ```LargeNumber``` you select, the wider the ranges you can check but the more memory ```Data``` consumes. Eg for ```LargeNumber=50```, ```Data``` consumes 1.1GB of RAM and your available range includes at least ```rangeN1=rangeN2=rangeM1=rangeM2=10``` but errors out with ```rangeN1=rangeN2=rangeM1=rangeM2=20```

## Even faster ?

- Try using the parallel processing package ([Wiki](https://github.com/NickG-Math/C4-Homology/wiki)).

- If you are using an AMD CPU, however recent, you should know that MATLAB uses the Intel MKL for matrix computations, which is optimized for Intel CPUs. So you should change the default BLAS (Basic Linear Algebra Subprograms) MATLAB uses to an open source one like OpenBLAS (I am not sure if this is possible with MATLAB; I believe it is with Octave but I have not tested it).



## For more details please consult the wiki.
