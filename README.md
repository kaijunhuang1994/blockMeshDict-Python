# blockMeshDict-Python
Using Python to generate a blockMesh file



## NACA-foil

To run a single case, execute

```
cd NACA-foil
pyhton blockmeshdict.py {foil} {alpha_deg}
```

For example, to simulate the flow around a NACA 0012 at 8 degrees angle of attack, run

```
pyhton blockmeshdict.py 0012 8
```



## NACA-foil-body-fitted

This case is modified based on NACA-foil 

To run a single case, execute

```
cd NACA-foil-body-fitted
pyhton blockmeshdict.py {foil}
```

For example, to simulate the flow around a NACA 0012, run

```
pyhton blockmeshdict.py 0012
```

