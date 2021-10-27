-- heatshield.lua
-- Simple job-specification file for e4prep -- for use with Eilmer4

-- We can set individual attributes of the global data object.
config.title = "Heatshield"
print(config.title)
config.dimensions = 2
config.axisymmetric = true

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel("tpair.lua")
T_inf = 75; -- Temperature 
P_inf = 290; -- Pressure 
rho_inf = 3e-4; -- Density 
mass_fraction = {N2=0.767, O2=0.233}; -- Mass fraction of each molecule

Q = GasState:new{gm}  -- Gas State 
Q.T = T_inf;
Q.p = P_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
gm:updateThermoFromRHOT(Q)

print("GasModel set to thermally perfect air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p = P_inf, T = T_inf, velx=0}  -- initial state
inflow = FlowState:new{ p = P_inf , T= T_inf, velx=1006.65} -- inflow parameters

-- Demo: Verify Mach number of inflow and compute dynamic pressure.
print("pressure=", inflow.p)
print("T=", inflow.T, "altitude=", 32000, "sound speed= ", inflow.a, "Pressure=", inflow.p)
print("u=", inflow.u, "rho=", Q.rho)


-- Geometry  
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x =0.05, y=0.0}
c = Vector3:new{x=0.045,y=0.03}
d = Vector3:new{x=0.11117, y=0.14}
g = Vector3:new{x=0.0, y = 0.03}
h = Vector3:new{x=0.11117, y=0.2}
e = Vector3:new{x=0.11507, y=0.128}
f = Vector3:new{x=0.12096, y=0.14}
i = Vector3:new{x=0.13507, y = 0.2}
j = Vector3:new{x=0.10507, y=0.1}
k = Vector3:new{x= 0.22, y = 0.09}
l = Vector3:new{x=0.22, y =0.15}

-- Creating Lines
ab = Line:new{p0=a, p1=b} 
bc = Bezier:new{points={b, c, d}}
ag = Arc3:new{p0=a, pmid=g, p1=h}
hd = Line:new{p0=h, p1=d}
def = Arc3:new{p0=d, pmid=e, p1=f}
fi = Line:new{p0=i, p1=f}
ih = Arc3:new{p0=h, pmid=j, p1=i}
kl = Line:new{p0=l, p1=k}
ik = Line:new{p0=i, p1=l}
fl = Line:new{p0=f, p1=k}

-- Creating Blocks
quad0 = makePatch{north=hd, east=bc, south=ab, west=ag}
quad1 = makePatch{north=fi, east=def, south=hd, west=ih}
quad2 = makePatch{north=kl, east= fl, south=fi, west=ik}

-- Mesh the patches, with particular discretisation.
cf_wall_normal = RobertsFunction:new{end0=false, end1=true, beta=1.08}
n0=30; n1=90; n2 = 15; n3 = 40 
grid0 = StructuredGrid:new{psurface=quad0, niv=n0+1, njv=n1+1, cfList = {north=cf_wall_normal, south=wall_normal}}
grid1 = StructuredGrid:new{psurface=quad1, niv=n0+1, njv=n2+1, cfList = {north=cf_wall_normal, south=cf_wall_normal}}
grid2 = StructuredGrid:new{psurface=quad2, niv=n0+1, njv=n3+1, cfList = {north=cf_wall_normal, south=cf_wall_normal}}

-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grid0, initialState=initial} -- Initial conditions
blk1 = FluidBlock:new{grid=grid1, initialState=initial} -- Initial conditions
blk2 = FluidBlock:new{grid=grid2, initialState=initial} -- Initial conditions 

-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList['west'] = InFlowBC_Supersonic:new{flowState=inflow}  -- Inflow from the west edge
blk1.bcList['west']= InFlowBC_Supersonic:new{flowState=inflow}  -- Inflow from west edge
blk0.bcList['east'] = WallBC_NoSlip_Adiabatic:new{}  -- Wall condition adiabatic with no slip
blk1.bcList['east'] = WallBC_NoSlip_Adiabatic:new{}  -- Wall Condition Adiabatic with no slip
blk2.bcList['east'] = WallBC_NoSlip_Adiabatic:new{}  -- Wall Condition Adiabatic with no slip
blk1.bcList['west']=OutFlowBC_Simple:new{}     -- OutFlow from west edge
blk2.bcList['north']= OutFlowBC_Simple:new{}   -- OutFlow from north edge
blk2.bcList['west']= OutFlowBC_Simple:new{}    -- OutFlow from west edge



-- add history point stagnation point
setHistoryPoint{x=0.02, y=0.00}
-- add history point at edge 
setHistoryPoint{x=0.08117, y=0.14}

-- Do a little more setting of global data.
config.fixed_time_step=false -- Time steps not fixed 
config.max_time = 1.2e-3  -- seconds
config.max_step = 1000000 -- Max Steps 
config.dt_init = 1.0e-8   -- Time Steps
config.cfl_value = 0.5  
config.dt_plot = 2.5e-5 -- plot ever 2.5e-5 seconds
config.dt_history = 1.0e-6   -- History Plot
config.extrema_clipping = false
config.viscous = true   -- Viscosity 
config.flux_calculator="adaptive"

--dofile("sketch-domain.lua")
