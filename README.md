<h1> Heat Transfer of 18650 Li-ion Batteries </h1>
<b>Version 0.7.6</b><br>
This program can be used for calculating heat transfer of li-ion batteries. This program uses this equation:<br>
<p align="center"><img src="https://github.com/mguluerler/Heat-Transfer-of-18650-Li-ion-Batteries/blob/master/figures/usedformula.png"></p>
<h3><i>class solvingParams: </i></h3>
Before calculation, required and desired parameters are entered here.
<h3><i>class BattProp:</i></h3>
Battery properties are entered here.
<h3><i>class dataSavers:</i></h3>
For showing some datas in the end of calculation, some lists are opened here, changing this class may not be good for stability of program.
<h3><i>class nodeFormulas: </i></h3>
In this class, equations defined for each node types.<br>
<h3><i>class k_find, class h_find:</i></h3>
In this class, convection and conduction coefficients are calculated.
<h3><i>class R_find:</i></h3>
<b>def Area_finder:</b> <i>Calculates the area to be used in heat transfer between the selected point and calculated point.</i><br>
<b>def dt_dr, def dt_dz:</b> <i>Calculates other values for calculating heat transfer.</i>
<h3><i>class materials:</i></h3>
Calculates and saves materials properties for calculations.
<h4><i>def volume_node:</i></h4>
Calculates volume of the node to be calculated.
<h4><i>def e_generation:</i></h4>
Calculates heat generation for every time step.
<h3><i>class graphs:</i></h3>
Shows calculated datas.
<p align="center"><img src="https://github.com/mguluerler/Heat-Transfer-of-18650-Li-ion-Batteries/blob/master/figures/graphs.png"></p><br>
