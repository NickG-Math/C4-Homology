# C4-Homology
This project computes the  RO(C<sub>4</sub>) homology of a point as a Green Functor and more. This has been supplanted by the C++ library <a href="https://github.com/NickG-Math/Mackey">Mackey</a>.
 
Read the [Wiki](https://github.com/NickG-Math/C4-Homology/wiki) for the full documentation. What follows is an **FAQ**; if you want to verify the results of our paper, in your preferred finite range, look no further than what's below.

## What do I need to run this code?
All you need is a version of MATLAB (the base installation).
 I have only tested it with v. R2019a but it should probably work with previous versions of MATLAB or freeware like Octave that can run MATLAB code.

## How do I compute the RO(C<sub>4</sub>) homology of a point?
Download the repository and add it to your MATLAB path. Then type the following in the command line:

```
write_Data(1,1);
Data=load_Data;
test_Pure_Homology(8,7,0,Data);
```

The last command will perform three operations:

- Compute the homology of 
<img src="http://latex.codecogs.com/svg.latex?S^{n\sigma+m\lambda}" border="0"/> for n=0,...,8 and m=0,...,7

- Check the answer against the tables in our paper, sending out an error message if there is a mismatch

- Print the answer in the form 
```
The k homology of the (n,m) sphere is MackeyFunctorSymbol
```
where ```MackeyFunctorSymbol``` is our notation of the corresponding Mackey functor.

Since the MATLAB display output does not support Latex, ```MackeyFunctorSymbol``` is more or less the Latex *code* of the symbol from our paper. For example, ```overline Z/2``` stands for <img src="http://latex.codecogs.com/svg.latex?\overline{\langle \mathbb{Z}/2\rangle }" border="0"/>. 

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
The input variables ```rangeN1, rangeN2, rangeM1, rangeM2``` specify the sets that the exponents of the Euler and orientation classes a<sub>&sigma;</sub>, u<sub>2&sigma;</sub>, a<sub>&lambda;</sub>, u<sub>&lambda;</sub> are allowed to range in, respectively.

The command above does three things:

- Computes the product of these Euler and orientation classes and checks if it's a generator

- Checks the answer against our tables sending out an error message if not

- Prints the answer in a form that looks like

```
Top Generator verified: asigma^2*alambda^3
```

or like

```
Mid Generator verified: usigma^4*bar(ulambda^2)
```

There's no reason to check the bottom level generators as that just involves the *nonequivariant* homology of spheres.

To check the rest of the tables run any of the following

```
test_Pure_Cohomology_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
test_Sigma_Minus_Lambda_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
test_Lambda_Minus_Sigma_Mult(rangeN1,rangeN2,rangeM1,rangeM2,0,Data);
```


## Can I make the program run faster?

Significant speed improvements can be achieved by precomputing data so as to avoid repeat calculations. Just run the test functions with the third input being 1 eg
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
This is because these two functions involve triple box products and extra data must be precomputed. So you must first run
```
write_Data(1,LargeNumber);
Data=load_Data;
```
for a large enough ```LargeNumber``` depending on the ranges you are using.

The larger the ```LargeNumber``` you select the more precomputed data you have available and the wider the ranges you can check. ```LargeNumber=100``` should be more than enough for any reasonable ranges (increasing it much further will make ```write_Data(1,LargeNumber)``` a lot slower).

## How about even faster?

Try [parallel processing](https://github.com/NickG-Math/C4-Homology/wiki/Parallel-Processing).


## For more details please consult the [Wiki](https://github.com/NickG-Math/C4-Homology/wiki)!
