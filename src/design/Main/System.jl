connect(fluid, sink.port)
connect(sink.port, valve.port_a)
connect(valve.port_b, vol.port)
connect(vol.flange, mass.flange)
connect(valve.area, ramp.output)