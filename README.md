# afd_1d
numerical analysis of finite difference schemes:

How to use afd_1d:
------------------
1) Go to 'afd_1d' directory and type 'make clean' to remove previous object and executable files.
2) Then run the command 'make' to build the executable 'afd_1d'.
3) Edit 'input.in' to select the model problem and available finite-difference schemes. Finite-
difference schemes need to added in the code if not available and 'make' again.
4) Now run the code using command './afd_1d input.in'.

Code generates output file/files in vtk file format that can be visualized in VisIt or paraview.
Output file contains properties of the numerical method evaluated for range of 'kh' and 'Nc'
(and/or 'Pe').
