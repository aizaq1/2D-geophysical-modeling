import UWGeodynamics as GEO
from UWGeodynamics import visualisation as vis
import numpy as np



u = GEO.UnitRegistry

half_rate = 25.5 * u.millimeter / u.year
# model_length = 270e3 * u.meter
model_length  = 1. * u.kilometer
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1603.15 * u.degK
bodyforce = 3370 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT



Model = GEO.Model(elementRes=(180, 60), 
                  minCoord=(0. * u.kilometer, -75. * u.kilometer), 
                  maxCoord=(270. * u.kilometer, 15. * u.kilometer), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))



Model.outputDir="outputs_38_traction_free_boundary_neg5_by_0.1_weakZoneReRun"

Model.diffusivity = 9e-7 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)



import GenerateInitialStructures_for_Rey_et_al_2d as GIS

DATA = GIS.MakeStructure(threeD=False,twoD_lat=38.0,Model=Model,Moho='moho_data.txt',Topo='topo_data.txt',min_lat=27,max_lat=44,min_long=-118,max_long=-109,resolution_long=0.03333333333,resolution_lat=0.1,elsl=15,dcp=-75,mat=0,sat=5.,structure_resolution=10)



air = DATA.air
sair = DATA.sair
topo = DATA.topo
uc = DATA.uc
mc = DATA.mc
lc = DATA.lc
ml = DATA.ml
lc_cp = DATA.lc_cp
mc_cp = DATA.mc_cp
surfElevation = DATA.surfElevation



# weak = GEO.shapes.Box(3, -30, 0, 15)
weak = GEO.shapes.Box(3, -15, 0, 15)



weak_layer = Model.add_material(name = "weak", shape = weak)



from UWGeodynamics import visualisation as vis
Fig = vis.Figure(figsize=(1200,800))
Fig.Points(Model.swarm, Model.materialField, discrete=True)

Fig.show()



# import GenerateInitialStructures_for_Rey_et_al_2d as GIS

# DATA = GIS.MakeStructure(threeD=False,twoD_lat=39.0,Model=Model,Moho='moho_data.txt',Topo='topo_data.txt',min_lat=27,max_lat=44,min_long=-118,max_long=-109,resolution_long=0.03333333333,resolution_lat=0.1,elsl=15,dcp=-75,mat=0,sat=5.,structure_resolution=5)



# air = DATA.air
# sair = DATA.sair
# topo = DATA.topo
# uc = DATA.uc
# mc = DATA.mc
# lc = DATA.lc
# ml = DATA.ml
# lc_cp = DATA.lc_cp
# mc_cp = DATA.mc_cp
# surfElevation = DATA.surfElevation


# # Lat 39 N



# import visulization as vis
# Fig = vis.Figure(figsize=(1200,800))
# Fig.Points(Model.swarm, Model.materialField, discrete=True)

# Fig.show()



# import GenerateInitialStructures_for_Rey_et_al_2d as GIS

# DATA = GIS.MakeStructure(threeD=False,twoD_lat=36.0,Model=Model,Moho='moho_data.txt',Topo='topo_data.txt',min_lat=27,max_lat=44,min_long=-118,max_long=-109,resolution_long=0.03333333333,resolution_lat=0.1,elsl=15,dcp=-75,mat=0,sat=5.,structure_resolution=5)



# air = DATA.air
# sair = DATA.sair
# topo = DATA.topo
# uc = DATA.uc
# mc = DATA.mc
# lc = DATA.lc
# ml = DATA.ml
# lc_cp = DATA.lc_cp
# mc_cp = DATA.mc_cp
# surfElevation = DATA.surfElevation


# # Lat 36 N



# #import visulization as vis
# Fig = vis.Figure(figsize=(1200,800))
# Fig.Points(Model.swarm, Model.materialField, discrete=True)

# Fig.show()



air.density = 1. * u.kilogram / u.metre**3
air.diffusivity = 1e-6 * u.metre**2 / u.second
air.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)

sair.density = 1. * u.kilogram / u.metre**3
sair.diffusivity = 1e-6 * u.metre**2 / u.second
sair.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)



topo.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3

uc.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3

mc.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3

lc.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3

lc_cp.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3
mc_cp.radiogenicHeatProd = 7.67e-7 * u.watt / u.meter**3

topo.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

uc.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

mc.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

lc.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

lc_cp.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

mc_cp.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

ml.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)



weak_layer.density = GEO.LinearDensity(reference_density=2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)



rh = GEO.ViscousCreepRegistry()



Model.minViscosity = 1e+18 * u.pascal * u.second
Model.maxViscosity = 5e+23 * u.pascal * u.second

air.viscosity = 5e+18 * u.pascal * u.second
sair.viscosity = 5e+18 * u.pascal * u.second

topo.viscosity = 1.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)  

uc.viscosity = 1.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)  

mc.viscosity = 1.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)  

lc.viscosity = 1.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)  


ml.viscosity = 1.5 * GEO.ViscousCreep(name='Dry Olivine (Brace and Kohlstedt, 1995)',
                                 preExponentialFactor=7e+4/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=5.2e5 * u.joule/u.mole,
                                 f=1.0)  

mc_cp.viscosity = 20.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)

lc_cp.viscosity = 20.0 * GEO.ViscousCreep(name='Wet Quartz Dislocation (Goetze et al, 1978)',
                                 preExponentialFactor=5e-6/u.megapascal ** 3. /u.second,
                                 stressExponent=3.,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=1.9e5 * u.joule/u.mole,
                                 f=1.0)



weak_layer.viscosity = 1e+18 * u.pascal * u.second #minimum



# topo.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)

# uc.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)

# mc.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)

# lc.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)

# ml.plasticity = GEO.DruckerPrager(name="Upper Mantle",
#                                            cohesion=15. * u.megapascal,
#                                            cohesionAfterSoftening=3. * u.megapascal,
#                                            frictionCoefficient=0.44,
#                                            frictionAfterSoftening=0.088,
#                                            epsilon1=0.0, epsilon2=0.5)


# lc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)

# mc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
#                                                 cohesion=15. * u.megapascal,
#                                                 cohesionAfterSoftening=3. * u.megapascal,
#                                                 frictionCoefficient=0.44,
#                                                 frictionAfterSoftening=0.088,
#                                                 epsilon1=0.0, epsilon2=0.5)



topo.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

uc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

mc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

lc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

ml.plasticity = GEO.DruckerPrager(name="Upper Mantle",
                                           cohesion=5 * u.megapascal,
                                           cohesionAfterSoftening=1 * u.megapascal,
                                           frictionCoefficient=0.1,
                                           frictionAfterSoftening=0.01,
                                           epsilon1=0.0, epsilon2=0.5)


lc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

mc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5 * u.megapascal,
                                                cohesionAfterSoftening=1 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)



weak_layer.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=1 * u.megapascal,
                                                cohesionAfterSoftening=0.1 * u.megapascal,
                                                frictionCoefficient=0.01,
                                                frictionAfterSoftening=0.001,
                                                epsilon1=0.0, epsilon2=0.5)



Model.set_temperatureBCs(top=293.15 * u.degK, materials=[(air, 293.15*u.degK), (sair, 293.15*u.degK)], bottom=1380. * u.degK)



Model.set_heatFlowBCs(bottom=(-200 * u.microwatt / u.metre**2, ml))



# Model.set_velocityBCs(left=[(-2.0) * u.millimeter / u.year, 0.0 * u.millimeter / u.year],
Model.set_velocityBCs(left=[(0.0) * u.millimeter / u.year, None],
#                       right=[0.0 * u.millimeter / u.year, 0.0 * u.millimeter / u.year],
                      right=[(0.0) * u.millimeter / u.year, None],
                     top=[None, 0.0 * u.millimeter / u.year])
#                        right=[0.0 * u.millimeter / u.year, 0.0 * u.millimeter / u.year])
#                       right=[(25.5/1) * u.millimeter / u.year, None],
#                       bottom=GEO.LecodeIsostasy(reference_mat=ml, average=True))



Model.set_stressBCs(bottom=[5 * u.megapascal, 0.1 * u.megapascal], top = [None, None])
# Model.set_stressBCs(bottom=[-5 * u.megapascal, 0.1 * u.megapascal], top = [None, None])



sediment = DATA.sediment
sediment.density = 2300. * u.kilogram / u.metre**3
sediment.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3
sediment.diffusivity = 1.0e-6 * u.metre**2 / u.second  
sediment.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
sediment.viscosity = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
pl = GEO.PlasticityRegistry()
sediment.plasticity = pl.Huismans_et_al_2011_Crust
sediment.plasticity.epsilon1 = 0.01
sediment.plasticity.epsilon2 = 1.0



Model.init_model()



FigVisc = vis.Figure(figsize=(1200,400))
FigVisc.Points(Model.swarm, GEO.dimensionalise(Model.viscosityField, u.pascal*u.second), fn_size=3.0, logScale=True)
FigVisc.show()



import matplotlib.pyplot as plt

if GEO.nProcs == 1:

    distances, viscosities = GEO.extract_profile(Model.projViscosityField, 
                                                 line = [(200.* u.kilometer, 0.),
                                                         (200.* u.kilometer, Model.bottom)])

    Fig, ax1 = plt.subplots(1,1,figsize=(10,7))
    ax1.plot(GEO.dimensionalise(viscosities, u.pascal*u.second), GEO.dimensionalise(distances, u.kilometer))
    ax1.set_xscale("log")
    ax1.set_xlabel("Log Viscosity (Pa.s)")
    ax1.set_ylabel("Depth in kms")
    ax1.set_ylim(75, 0)
    ax1.set_xlim(1e20, 1e24)
    ax1.set_title("Viscosity Profile")



Model.solver.set_inner_method("mumps")
Model.solver.set_penalty(1e6)



Model.run_for(30000000.* u.years, checkpoint_interval=100000. * u.years,dt=100000.0*u.years)
















# 
