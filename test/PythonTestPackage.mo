within ;
package PythonTestPackage
  "Package used to generate .mat files for testing the python simman module"

connector HeatPort
Modelica.SIunits.Temperature T "potential";
flow Modelica.SIunits.HeatFlowRate Q_flow "flow variable";

/* Modelica.SIunits.pressure p;
  flow Modelica.SIunits.MassFlowrate m_flow;
*/
  annotation (Icon);
end HeatPort;

model Resistor
HeatPort heatPort_a;
HeatPort heatPort_b;
parameter Modelica.SIunits.ThermalResistance R=1;

equation
heatPort_a.Q_flow=(heatPort_a.T-heatPort_b.T)/R;
heatPort_a.Q_flow+heatPort_b.Q_flow=0;
end Resistor;

model Capacity
HeatPort heatPort;
parameter Modelica.SIunits.HeatCapacity C=800;
Modelica.SIunits.Temperature T;

equation
T=heatPort.T;
C*der(T)=heatPort.Q_flow;

end Capacity;

model LinkedCapacities
Capacity c1(T(start = 400),C=600);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=3);

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities;

model Array "Test model to check array handling in python"

parameter Integer n=4;
parameter Real[n] cvalues={600,1000,400,600};
parameter Real[n] Tstarts={500,350,600,200};
Capacity[n] c(C=cvalues, T(start = Tstarts));
Resistor[n] r;

/*Resistor[m,n,o] r; "3D matrix"*/

equation
for i in 1:n loop
connect(c[i].heatPort,r[i].heatPort_a);

connect(r[1].heatPort_b,r[i].heatPort_b);
end for;

end Array;
  annotation (uses(Modelica(version="3.1")));
model LinkedCapacities_A
Capacity c1(T(start = 400),C=800);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=3);
parameter String info_string = "Ref with changed c1.C";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_A;

model LinkedCapacities_B
Capacity c1(T(start = 400),C=1000);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=3);
parameter String info_string = "Ref with changed c1.C";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_B;

model LinkedCapacities_C
Capacity c1(T(start = 400),C=600);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=5.5);
parameter String info_string = "Ref with changed r.R";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_C;

model LinkedCapacities_D
Capacity c1(T(start = 400),C=800);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=5.5);
parameter String info_string = "as A but changed r.R";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_D;

model LinkedCapacities_E
Capacity c1(T(start = 400),C=800);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=8.15);
parameter String info_string = "as A but changed r.R";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_E;

model LinkedCapacities_D_TooShort
Capacity c1(T(start = 400),C=800);
Capacity c2(T(start = 350),C=1000);
Resistor r(R=5.5);
parameter String info_string = "as A but changed r.R";

equation
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_D_TooShort;

model Array_Big
    "Test model to check array handling in python - big array and big .mat file"

parameter Integer n=800;

Capacity[n] c;
Resistor[n] r;

/*Resistor[m,n,o] r; "3D matrix"*/

equation
for i in 1:n loop
connect(c[i].heatPort,r[i].heatPort_a);

connect(r[1].heatPort_b,r[i].heatPort_b);
end for;

end Array_Big;

model LinkedCapacities_F
Capacity c1(T(start = 400),C=800);
Capacity c2(T(start = 350),C=1000);
Resistor r(R = 3);
parameter String info_string = "as A with additional parameter";
parameter Real newParameter = 5 "a new parameter added here";
Real newVariable(start = 25);

equation
newVariable = 25 - time;
connect(c1.heatPort,r.heatPort_a);
connect(c2.heatPort,r.heatPort_b);

/*
connect-statement :
  - potentials are set equal,
  - Kirchoff to "flow" variables
*/

    annotation (experiment(StopTime=10000, Interval=200), experimentSetupOutput);
end LinkedCapacities_F;
end PythonTestPackage;
