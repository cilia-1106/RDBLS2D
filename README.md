# RDBLS2D
A reaction diffusion-based B-spline level set (RDBLS) method for structural topology optimization

This repository aims to introduce RDBLS2D, a 99-line MATLAB code that integrates the advantages of the B-spline based level set function and the reaction–diffusion update scheme to achieve efficient and effective structural topology optimization. It uses 66 and 33 lines for the main program and the B-spline level set representation method updated by the reaction–diffusion scheme, respectively. Also, repeated B-spline's knots on the boundary of the design domain were used in the code, which can naturally guarantee the connection between neighboring cells in the design of functionally graded materials.

Here are the descriptions of each file:
"RDBLS2D.m" is the standard code with B-spline simple knots.
"RDBLS2D_repeated.m" is the standard code with B-spline repeated knots, and "NN.mat"** is a mandatory file to run this program.
All the other files are the codes for the examples in Ref. [1] with default parameters.

Please cite the following article if you use this code in your publications:

[1] C. Wang, Y.M. Xie, X. Huang, X. Zhang, S. Zhou, MATLAB codes of the parametrized level set method for structural topology optimization using B-spline’s simple or repeated knots, Structural and Multidisciplinary Optimization, 67 (2024).
[2] C. Wang, Y.M. Xie, X. Lin, S. Zhou, A reaction diffusion-based B-spline level set (RDBLS) method for structural topology optimization, Computer Methods in Applied Mechanics and Engineering, 398 (2022) 1–36.

------------------------------------------------------------------------------------------------------------
**The code used to generate "NN.mat" is being sorted out and will be uploaded soon.
