# IonExplorer2
# Overview
IonExplorer2 is a tool for automatic calculations of ion migration barriers and trajectories in crystalline ion conductors.
# How to run
The minimal required input information is the Crystallographic Information File (CIF) https://www.iucr.org/resources/cif and the type of migrating ion. An example
```
python IonExplorer2.py –i Li6PS5Cl.cif -m Li
```
To see the complete list of available options
```
python IonExplorer2.py –h
```
To make IonExplorer2.py executable:

1. Add a shebang line with a path to the python interpreter to the top of the file
```
#!/home/.../bin/python
```
2. Mark the file as executable
```
chmod +x IonExplorer2.py
```
3. Add a path to the directory with the program to your PATH variable. Execute the command below or add it to .bashrc or .bash_profile in your home directory.
```
export PATH=/path/to/IonExplorer2:$PATH
```
# Dependencies
To run the program, you need to have:
-	Python3 https://www.python.org/downloads/release/python-385/
-	Scipy https://docs.scipy.org/doc/scipy-1.7.1/reference/
-   Matplotlib https://matplotlib.org/
-	PyCifRW  https://pypi.org/project/PyCifRW/4.4.1/
-	Critic2 https://github.com/aoterodelaroza/critic2
# How to cite
1. Golov, A., Carrasco, J. Enhancing first-principles simulations of complex solid-state ion conductors using topological analysis of procrystal electron density. npj Comput Mater 8, 187 (2022). https://doi.org/10.1038/s41524-022-00877-6
2. Zolotarev, P.N., Golov, A.A., Nekrasova, N.A., Eremin R.A. Topological analysis of procrystal electron densities as a tool for computational modeling of solid electrolytes: A case study of known and promising potassium conductors. AIP Conf. Proc. 2163, 020007 (2019). https://doi.org/10.1063/1.5130086
