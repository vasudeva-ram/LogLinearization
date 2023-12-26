# Overview
The **LogLinearizer** package generates symbolic log-linearized representations of dynamic nonlinear equations using a first-order Taylor Series approximation.
The tool is designed with macroeconomic DSGE models in mind, and supports the DYNARE notation.
However, while DYNARE _solves_ the model, this package focuses on providing the symbolic linearized equations, which are often helpful in driving intuition.

# How to use the LogLinearization package
There are two ways to use this package. 
The first method assumes that you have a fully built out DYNARE mod file.
This means, for each variable and each parameter, the file provides a TeX representation, steady state values and so on.
The second method is based on providing the necessary dictionaries to the package separately, and is more convenient when you only have a couple of equations to log-linearize.
Below, I provide examples and use cases for each method.

## Option 1: Using a .mod file
For this method, I assume that you have a fully written out .mod file available with the DYNARE code fully written out.
For best results, this .mod file should have the following features:
1. Defined TeX representations for each variable and parameter defined in the file, like so
``` 
var C   ${C}$
    Pi  ${\Pi}$
    A   ${A}$
    ...
    
parameters  alppha  ${\alpha}$
            betta   ${\beta}$
            ...
                
```
This allows the resulting log-linearized equations to look like LaTeX output.
In the absence of these TeX representations, the package uses the name of the variable/parameter directly, but this ends up representing percentage deviations of inflation from its steady state as $\hat{Pi}$ instead of $\hat{\Pi}$, which is worse looking.

2. A defined `steady_state_model` block
Providing steady state values for some of the variables helps the simplification process enormously and makes the end result look cleaner.
Even if you don't have the steady state fully derived, if you have some specific steady state values, you should provide those.
For example, if you are linearizing around a zero-inflation steady state, you gross inflation rate would be equal to 1. You should provide this value.
**IMPORTANT**: This package _log_-linearizes the equations in your model.
This means that the steady state values of your variables cannot be equal to zero (since $\log{0}$ is not defined).
If you have variables taking on a steady-state value of zero, redefine the variable to be equal to 1.
For example, rather than use a net inflation rate $\pi$ of zero, use a gross inflation rate $\Pi = 1 + \pi$ with a steady-state value of 1.

### Generating the log-linearized equations
For this tutorial, I will use a [(]simple NK model](https://github.com/JohannesPfeifer/DSGE_mod/blob/master/Gali_2008/Gali_2008_chapter_2.mod) with capital from Johannes Pfeifer's excellent collection of DSGE models in DYNARE.
Start by importing this file from wherever you have saved it.


```python
import LogLinearization as ll
```

    IPython console for SymPy 1.12 (Python 3.11.4-64-bit) (ground types: python)
    
    These commands were executed:
    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)
    >>> init_printing()
    
    Documentation can be found at https://docs.sympy.org/1.12/
    


Next, call the `createModEconomy()` function on this file.
This returns two objects: (1) a `ModEconomy` object which contains variables, parameters and steady-state values and their dynamic representations, and (2) a list of `Equation` objects that represent the equations of the model.


```python
nk, eqns = ll.createModEconomy('SimpleNK.mod')
```

`nk` here is the `ModEconomy` object, while eqns is a `list` containing the `Equation` objects.
To generate the log-linearized versions of these equations, I use the `loglin()` method of the Equation objects:


```python
for eqn in eqns[:6]: #only first 7 equations; more on this later
    display(eqn.loglin())
```


$\displaystyle \hat{{\frac{W}{P}}}_t = \hat{{C}}_t {\sigma} + \hat{{N}}_t {\phi}$



$\displaystyle - \hat{{R}}^n_t = \hat{C_t} {\sigma} - \hat{C_{t+1}} {\sigma} - \hat{{\Pi}}_{t+1}$



$\displaystyle \hat{{A}}_t + \hat{{N}}_t \left(1 - {\alpha}\right) = \hat{{C}}_t$



$\displaystyle \hat{{\frac{W}{P}}}_t = \hat{{A}}_t - \hat{{N}}_t {\alpha}$



$\displaystyle \hat{R^r_t} = \hat{R}^n_t - \hat{{\Pi}}_{t+1}$



$\displaystyle \hat{R^n_t} = \frac{\bar{{\Pi}}^{{\phi_{\pi}}} \hat{\Pi_t} {\phi_{\pi}} + \bar{{\varepsilon_m}} \hat{\varepsilon^m_t} {\beta}}{\bar{{\Pi}}^{{\phi_{\pi}}} + \bar{{\varepsilon_m}} {\beta}}$


## Option 2: Building a ModEconomy object

Sometimes, you may only have a couple of equations involving a few variables and parameters to log-linearize, in which case building an entire .mod file in DYNARE is a bit too much work.
For such cases, you can build out a `ModEconomy` object from its component pieces quite easily.
- Step 1: Define which symbols in your equation are variables, which are parameters and provide their respective representations.
- Step 2: Provide steady-state values for any variables which have simple steady-state values
You can do this in the form of dictionaries.

Let's say we want to log-linearize only the Euler equation from the previous .mod file.
``` 
euler = '1/R=betta*(C(+1)/C)^(-siggma)/Pi(+1)'
```
The only variables are `C`, `Pi` and `R`. 
The parameters are `betta`, and `siggma`.
Let's say we are linearizing around a zero-inflation steady state, such that `Pi = 1` in the steady state.

Let's give them the following representations in a dictionary:


```python
euler = '1/R=betta*(C(+1)/C)^(-siggma)/Pi(+1)'
vardict = {'C':'C', 'R':'R', 'Pi':r'\Pi'}
paramdict ={'siggma':r'\sigma', 'betta':r'\beta'}
ssdict = {'Pi':1};
```

Remember that when using a TeX representation for a variable, you should prefix the string with an `r` to indicate that it should be interpreted as a raw string.

Now we create a `ModEconomy` object and an `Equation` object separately (something the `createModEconomy()` function did automatically for us from the .mod file).


```python
econ = ll.ModEconomy(vardict, paramdict, ssdict)
eqn = ll.Equation(euler, econ)
display(eqn.loglin())
```


$\displaystyle - \hat{R_t} = \hat{C_t} \sigma - \hat{C_{t+1}} \sigma - \hat{\Pi}_{t+1}$


# Some caveats on what the package does not handle
Notice that I only log-linearized the first 7 equations earlier when using the mod file. 
This is for reasons I elaborate on here.
 
1. Equations that are linear in logs. 
Often in DSGE models, "shocks" are assumed to follow a process that is linear in logs. 
For example, in the DSGE file used above, equation 8 reads:
``` log(A)=rho*log(A(-1))+eps_A; ``` 
The package currently does not handle these.

2. Variables with a zero steady state value.
Equation 9 in the mod file reads uses a variable called `m_growth_ann` which is defined to have a 
steady state value of 0.
Since log-linearization requires taking logs, applying this prcocedure to this equation will yield a `NaN` or incorrect result.

3. "Simplification" is not a well-defined mathematical operation!!
The most useful simplified representation of an equation depends on what you are looking to solve for.
So some of the simplified log-linearized versions returned by this package may not immediately be in the optimal form for you.



```python
# Some 
```
