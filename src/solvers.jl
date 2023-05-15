using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible: HydraulicFluid, Valve, FixedPressure, DynamicVolume
using ModelingToolkitStandardLibrary.Blocks: Ramp
using ModelingToolkitStandardLibrary.Mechanical.Translational: Mass

using ModelingToolkitDesigner

@parameters t
D = Differential(t)

pars = @parameters begin
    g = -9.807
    r = 0.1
    m = 10
end

vars = @variables begin 
    x(t)=0.0
    dx(t)=0
    ddx(t)=0
    y(t) = -0.1
    dy(t) = 0
    ddy(t) = g
    θ(t) = -π/2
    T(t) = 0
end



eqs = [
    D(x) ~ dx
    D(dx) ~ ddx
    D(y) ~ dy
    D(dy) ~ ddy

    r*sin(θ) ~ y
    r*cos(θ) ~ x

    m*ddy ~ -T*sin(θ) + m*g
    m*ddx ~ -T*cos(θ)
    
]

@named system = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(system)
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 1.0))

sol = solve(prob, Rodas4())
sol = solve(prob, Rodas4(); initializealg = NoInit(), dt=1e-4, adaptive=false);



@component function System(; name)
    pars = []

    systems = @named begin
        fluid = HydraulicFluid()
        sink = FixedPressure(; p = 10e5)
        vol = DynamicVolume(5;p_int = 100e5, area=0.001, x_int=0.05, x_max=0.1, x_damp=0.02, x_min=0.01, direction=+1 )
        valve = Valve(; p_a_int = 10e5, p_b_int = 100e5, area_int = 0, Cd = 1e6)
        ramp = Ramp(; height = 0.001, duration = 0.001, offset = 0, start_time = 0.001,
                      smooth = true)
        mass = Mass(; m=1000, g=-9.807, s_0=0.1)
    end

    eqs = [connect(fluid, sink.port)
           connect(sink.port, valve.port_a)
           connect(valve.port_b, vol.port)
           connect(vol.flange, mass.flange)
           connect(valve.area, ramp.output)]

    ODESystem(eqs, t, [], pars; name, systems)
end

@named system = System()
s = complete(system)
sys = structural_simplify(system)
prob = ODEProblem(sys, [], (0, 1.0))
solve(prob, Rodas4())


NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt = 1e-4,
            initializealg = NoInit())

using CairoMakie
CairoMakie.activate!()
lines(sol[s.mass.s])

path = joinpath(@__DIR__, "design") # folder where visualization info is saved and retrieved
design = ODESystemDesign(system, path);
ModelingToolkitDesigner.view(design)