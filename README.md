# Heat Transfer of 18650 Li-ion Batteries
**Version 0.7.6** <br>
This program can be used for calculating heat transfer of li-ion batteries. This program uses this equation:<br>
![equation figure](https://github.com/mguluerler/Heat-Transfer-of-18650-Li-ion-Batteries/blob/master/li-ion%20batteries/figures/usedformula.png)<br>
### class solvingParams:
Before calculation, required and desired parameters are entered here.
### class BattProp:
Battery properties are entered here.
### class dataSavers:
For showing some datas in the end of calculation, some lists are opened here, changing this class may not be good for stability of program.
### class nodeFormulas: 
In this class, equations defined for each node types.<br>
### class k_find, class h_find:
In this class, convection and conduction coefficients are calculated.
### class R_find:
**def Area_finder:** *Calculates the area to be used in heat transfer between the selected point and calculated point.*<br>
**def dt_dr, def dt_dz:** *Calculates other values for calculating heat transfer.*
### class materials:
Calculates and saves materials properties for calculations.
#### def volume_node:
Calculates volume of the node to be calculated.
#### def e_generation:
Calculates heat generation for every time step.
### class graphs:
Shows calculated datas.
![graphs](https://github.com/mguluerler/Heat-Transfer-of-18650-Li-ion-Batteries/blob/master/li-ion%20batteries/figures/graphs.png)
