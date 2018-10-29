model ExampleForESTEC
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T(displayUnit = "K") = 0.05) annotation(
    Placement(visible = true, transformation(origin = {-20, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Sources.Force force1 annotation(
    Placement(visible = true, transformation(origin = {-52, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.SpringDamper springDamper1(c = 1, d = 2, useHeatPort = true)  annotation(
    Placement(visible = true, transformation(origin = {22, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Mass mass1(m = 1)  annotation(
    Placement(visible = true, transformation(origin = {-14, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Translational.Components.Fixed fixed1 annotation(
    Placement(visible = true, transformation(origin = {66, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Pulse pulse1(amplitude = 0.01,nperiod = 1, period = 1, width = 100)  annotation(
    Placement(visible = true, transformation(origin = {-86, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.HeatingResistor heatingResistor1(R_ref = 100, useHeatPort = true)  annotation(
    Placement(visible = true, transformation(origin = {-72, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(visible = true, transformation(origin = {-42, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage1(V = 0.01, freqHz = 1)  annotation(
    Placement(visible = true, transformation(origin = {-72, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ModelicaRect_8_quad modelicaRect_8_quad1(Tstart = 0.05)  annotation(
    Placement(visible = true, transformation(origin = {10, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(fixedTemperature1.port, modelicaRect_8_quad1.IF_2) annotation(
    Line(points = {{-10, -42}, {10, -42}, {10, -10}, {10, -10}}, color = {191, 0, 0}));
  connect(heatingResistor1.heatPort, modelicaRect_8_quad1.IF_1) annotation(
    Line(points = {{-72, -18}, {-14, -18}, {-14, -2}, {10, -2}, {10, -2}}, color = {191, 0, 0}));
  connect(springDamper1.heatPort, modelicaRect_8_quad1.IF_1) annotation(
    Line(points = {{12, 54}, {10, 54}, {10, -2}, {10, -2}}, color = {191, 0, 0}));
  connect(sineVoltage1.n, ground1.p) annotation(
    Line(points = {{-62, 24}, {-42, 24}, {-42, 20}}, color = {0, 0, 255}));
  connect(sineVoltage1.p, heatingResistor1.p) annotation(
    Line(points = {{-82, 24}, {-82, -8}}, color = {0, 0, 255}));
  connect(sineVoltage1.n, heatingResistor1.n) annotation(
    Line(points = {{-62, 24}, {-62, -8}}, color = {0, 0, 255}));
  connect(springDamper1.flange_b, fixed1.flange) annotation(
    Line(points = {{32, 64}, {66, 64}, {66, 66}, {66, 66}}, color = {0, 127, 0}));
  connect(mass1.flange_b, springDamper1.flange_a) annotation(
    Line(points = {{-4, 62}, {12, 62}, {12, 64}, {12, 64}}, color = {0, 127, 0}));
  connect(force1.flange, mass1.flange_a) annotation(
    Line(points = {{-42, 62}, {-26, 62}, {-26, 62}, {-24, 62}, {-24, 62}}, color = {0, 127, 0}));
  connect(pulse1.y, force1.f) annotation(
    Line(points = {{-74, 62}, {-66, 62}, {-66, 62}, {-64, 62}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.2")));
end ExampleForESTEC;