# MasterFunctions

This is a Wolfram Language package used in most of my scattering amplitude computations. 
Many of the functions related to graph generation, graph isomorphisms, 
pure YM Feynman rules, momentum conservation, relabeling, as well as
many of the structure definitions are thanks to Tristan Dennen. 
Much of this code builds off of his great foundation. Many thanks 
to Tristan Dennen, Scott Davies and Zvi Bern for code and/or the 
physical insights therein.

A file with examples is also provided. These examples include basic calculations as well as the 
verification of results presented in [arXiv: 1309.7416](http://arxiv.org/abs/1309.7416) 
and [arXiv: 1510.03448](http://arxiv.org/abs/1510.03448). These verifications demonstrate the 
use of many of analytical and numerical tools in the package.

Some functions are cleaner and more elegant than others. That being said, most
functions have been verified to reproduce known results in a variety of
different ways. I am sure there are bugs as well as scenarios in which
some functions will fail. Most functions were developed with specific
applications in mind. Please feel free to contact me with any questions!

Notes:
- The `aux` subdirectory stores previously generated results for later retrieval.
- Be sure to set the variable `MASTERFUNCTIONSDIRECTORY` in `Public Initializations` section of `MasterFunctions.m` to the appropriate directory.
- The function `MonsterSolve2` has been removed since it is a proprietary function. This should only affect color information generation if the result hasn't already been saved. Contact me for more information.
