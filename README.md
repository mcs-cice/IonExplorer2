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

add a shebang line with a path to the python interpreter to the top of the file
```
#!/home/.../bin/python
```
mark the file as executable
```
chmod +x IonExplorer2.py
```
add a path to the directory with the program to your PATH variable. Execute the command below or add it to .bashrc or .bash_profile in your home directory.
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
