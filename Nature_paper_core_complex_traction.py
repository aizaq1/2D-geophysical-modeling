import UWGeodynamics as GEO
import numpy as np
import math
import underworld.function as fn
import os 



GEO.rcParams["swarm.particles.per.cell.2D"]= 20



min_lat = 27
max_lat = 44
min_long = -124
max_long = -103
spacing_long = 0.05
spacing_lat = 0.1

scale_coef = 1

#  define spacing for making your mesh in x and y.
# n = 1
n = 2

#  define spacing for making your mesh in z direction in km.
# nz = 0.5
nz = 1
# spacing for making 2D structures like topo, crust, etc.
n_2 = 1
# spacing for making badlands model. It should be similar to n_2 for the first step of simulations, during first and sedcond steps.
n_3 = 1

dfs = int(n)
dfs_2 = int(n_2)
dfs_3 = int(n_3)

latitude_input = 38.0

# provide information for your 2D model.
# maximum_elevation_of_2D_model_from_sealevel km
elsl = 15
# depth_of_compensation_for_2D_model km
dcp = -75.
# mantle_asthenosphere_thickness_for_2D_model km
mat = 0.
# upper crust bottom depth km
ucd = -10.0
# middle crust bottom depth km
mcd = -25.0
# thickness of sticky air
sat = 5.0

# time interval for horizontal velocity (years)
t_v = 500000
# time interval for solving the model for (years)
t_s = 500000

# velocity magnitude threshold for fix points to be less than
v_limit = +0.11

# if choosing boundary velocity from continental litho then put the crustal thickness limit km. If it doesn't matter then put 0.0
c_limit = 21

# if adding temperature to boudnary condition temperature at the bottom of mesh then add a value. If not then put 0.0
heat_add = mat

# air and sticky air density control to multiply density of 1 by that factor.
factor = 1000

print("longitude range is: {}, {}".format(min_long, max_long))
print("latitude range is: {}, {}".format(min_lat, max_lat))
print("spacing long x and y axes is: {}, {} degree".format(spacing_long, spacing_lat))
print("you are making a 2D profile for latitude: {}".format(latitude_input))



ll = round(latitude_input)

if ll > latitude_input:
    lat_decimal_profile = 1.0 - (ll - latitude_input)
    lat_int_profile = latitude_input - lat_decimal_profile
elif ll == latitude_input:
    lat_int_profile = latitude_input
    lat_decimal_profile = 0.0
else:
    lat_decimal_profile = latitude_input - ll
    lat_int_profile = latitude_input - lat_decimal_profile
    
la = lat_int_profile + lat_decimal_profile
inc_t = t_v/t_s



spacing_x = 1/spacing_long
spacing_y = 1/spacing_lat
resx = int(((max_long-min_long) * spacing_x)/n)
resy =  int((((elsl + (-1*dcp) + mat))/nz)*1)
Res_mesh = (resx,resy)
model_top = int(elsl)
inc_t = t_v/t_s

print("grid resolution is: resx: {}, resy: {}".format(resx,resy))



u = GEO.UnitRegistry



half_rate = 25.5 * u.millimeter / u.year
# model_length is the distance between each Res_mesh
model_length  = 1. * u.kilometer
surfaceTemp   = 273.15 * u.degK # 0 * u.degC
baseModelTemp = 1633.15 * u.degK # 1360 * u.degC
bodyforce = (3370 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2) 

KL = model_length
Kt = KL / half_rate
# Kt = 1. * u.year
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT



# Model mesh

Model = GEO.Model(elementRes=Res_mesh,
                  minCoord=(0.0 * scale_coef * u.kilometer, (-1*((resy*nz)-model_top)) * u.kilometer),  
                  maxCoord=(resx*n * scale_coef * u.kilometer, model_top * u.kilometer), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))

print(len(Model.mesh.data))
#print(Model.mesh.data)



Model.outputDir="Output_WUS_2D"



import GenerateInitialStructures2 as GIS



DATA = GIS.MakeStructure(threeD=False,
                         twoD_lat=latitude_input,
                         Model=Model,
                         Moho='Moho_depth_coordinates_36Ma_for_UWG.txt',
                         Topo='elevation_coordinates_36Ma_for_UWG.txt',
                         resolution_long=spacing_long,
                         resolution_lat=spacing_lat,
                         min_lat=min_lat,max_lat=max_lat,min_long=min_long,max_long=max_long,
                         mat=mat,sat=sat,elsl=elsl,ucd=ucd,mcd=mcd,
                         scale_coef=scale_coef,
                         structure_resolution=dfs_2,
                         badlands_resolution=dfs_3)




oc = DATA.oc
air = DATA.air
water = DATA.water
sair = DATA.sair
topo = DATA.topo
uc = DATA.uc
mc = DATA.mc
lc = DATA.lc
ml = DATA.ml
lc_cp = DATA.lc_cp
mc_cp = DATA.mc_cp
surfElevation = DATA.surfElevation
sediment = DATA.sediment



Model.diffusivity = 9e-7 * u.metre**2 / u.second
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)



air.density = 1. * u.kilogram / u.metre**3
air.diffusivity = 1e-6 * u.metre**2 / u.second
air.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)

sair.density = 1. * u.kilogram / u.metre**3
sair.diffusivity = 1e-6 * u.metre**2 / u.second
sair.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)

water.density = 1. * u.kilogram / u.metre**3
water.diffusivity = 1e-6 * u.metre**2 / u.second
water.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)

ml.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)





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



rh = GEO.ViscousCreepRegistry()



Model.minViscosity = 1e+18 * u.pascal * u.second
Model.maxViscosity = 5e+23 * u.pascal * u.second

oc.viscosity = 5e+19 * u.pascal * u.second
air.viscosity = 5e+18 * u.pascal * u.second
sair.viscosity = 5e+18 * u.pascal * u.second
water.viscosity = 8e+18 * u.pascal * u.second

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




topo.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.3,
                                                frictionAfterSoftening=0.03,
                                                epsilon1=0.0, epsilon2=0.5)

uc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.3,
                                                frictionAfterSoftening=0.03,
                                                epsilon1=0.0, epsilon2=0.5)

mc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.3,
                                                frictionAfterSoftening=0.03,
                                                epsilon1=0.0, epsilon2=0.5)

lc.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.3,
                                                frictionAfterSoftening=0.03,
                                                epsilon1=0.0, epsilon2=0.5)

ml.plasticity = GEO.DruckerPrager(name="Upper Mantle",
                                           cohesion=5. * u.megapascal,
                                           cohesionAfterSoftening=1. * u.megapascal,
                                           frictionCoefficient=0.3,
                                           frictionAfterSoftening=0.03,
                                           epsilon1=0.0, epsilon2=0.5)

lc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=15. * u.megapascal,
                                                cohesionAfterSoftening=3. * u.megapascal,
                                                frictionCoefficient=0.44,
                                                frictionAfterSoftening=0.088,
                                                epsilon1=0.0, epsilon2=0.5)

mc_cp.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=15. * u.megapascal,
                                                cohesionAfterSoftening=3. * u.megapascal,
                                                frictionCoefficient=0.44,
                                                frictionAfterSoftening=0.088,
                                                epsilon1=0.0, epsilon2=0.5)




sediment.density = 2300. * u.kilogram / u.metre**3
sediment.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3
sediment.diffusivity = 1.0e-6 * u.metre**2 / u.second
sediment.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
sediment.viscosity = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
pl = GEO.PlasticityRegistry()
sediment.plasticity = pl.Huismans_et_al_2011_Crust
sediment.plasticity.epsilon1 = 0.01
sediment.plasticity.epsilon2 = 1.0




import pandas as pd
import numpy as np

with open("Moho_depth_coordinates_36Ma_for_UWG.txt", "r") as input_crust:
    with open("12.txt", "w") as output_thickness:
        with open("elevation_coordinates_36Ma_for_UWG.txt", "r") as input_elevation:
            with open("13.txt", "w") as output_elevation:    
                crust = pd.read_csv(input_crust, dtype=float, delimiter=",",skiprows=0)
                df_crust = pd.DataFrame(crust)
                df_crust = np.array(df_crust)
                x = df_crust[:,0]
                y = df_crust[:,1]
                z = df_crust[:,2]
                elevation = pd.read_csv(input_elevation, dtype=float, delimiter=",",skiprows=0)
                df_elevation = pd.DataFrame(elevation)
                df_elevation = np.array(df_elevation)
                k = df_elevation[:,2]
                f = len(x)
                q1 = x.reshape(f,1)
                q2 = y.reshape(f,1)
                q3 = z.reshape(f,1)
                q4 = k.reshape(f,1)
                t = np.hstack((q1,q2,q3.round(),q4.round()))
                df2 = pd.DataFrame(t)
                header = " x    y   z   k"
                print("{}".format(header), file=output_thickness)
            
                elevation = ((k/1000).round())
                Moho = ((z/1).round())
                index1 = q1
                print (Moho)
                print (elevation)
                print (index1)
                
                coor_elev = []
                coor_moho = []
                moho_data = []
                topo_data = []
                
                for data in range(len(index1)):
                
                    ele_data = elevation[data]
                    moho_data1 = Moho[data]
                    x_coor = x[data]
                    y_coor = y[data]
                    coor_elev = (x_coor,y_coor,ele_data)
                    topo_data.append(coor_elev)
                    coor_moho = (x_coor,y_coor,moho_data1)
                    moho_data.append(coor_moho)


# Set python variables for model parameters.
dfs = int(n)
spacing_x = 1/spacing_long
spacing_y = 1/spacing_lat
resxx = int(((max_long-min_long) * spacing_x)/n)
resyy = int(((max_lat-min_lat) * spacing_y)/n)
resz = int((((elsl + (-1*dcp) + mat))/nz)*1)
Res_mesh = resx, resy, resz
model_top = int(elsl)
inc_t = t_v/t_s

print("grid resolution is: resxx: {}, resyy: {}, resz: {}".format(resxx,resyy,resz))





# TEMPERATURE Boundary conditions WITHIN THE WHOLE MANTLE ASTHENOSPHERE 
#Boundary conditions
import pandas as pd
import numpy as np
import math

with open("UWG_temp_36.0_Ma.txt", "r") as input_crust:
    with open("30.txt", "w") as output_crust:    
        crust = pd.read_csv(input_crust, dtype=float, delimiter=",",skiprows=0)
        df_crust = pd.DataFrame(crust)
        df_crust = np.array(df_crust)

        x = df_crust[:,0]
        y = df_crust[:,1]
        z = df_crust[:,2]
        f = len(x)
        q1 = x.reshape(f,1)
        q2 = y.reshape(f,1)
        q3 = z.reshape(f,1)
        index1 = q1
        temp = z + 273.15         
                
           
                
        temp_data0 = []
        coor_temp = []
                
        for data in range(len(index1)):
            temp_data = temp[data]

            x_coor = x[data]*scale_coef
            y_coor = y[data]*scale_coef
            coor_temp = (float(x_coor),float(y_coor),float(temp_data))
            temp_data0.append(coor_temp)
    
        

# HERE WE DEFINE WHICH LATITUDE WE WANT TO CHOOSE FOR MAKING 2D PROFILES BY DEFINING LATITUDE VALUE.
# THE RESOLUTION FOR LATITUDES IS 0.1 DEGREE FROM 27N TO 44N.

        lat_int = lat_int_profile
        lat_decimal = lat_decimal_profile 
        latitude = ((lat_int - min_lat) * 10) + (lat_decimal * 10)

            
        max_y = int((resyy*dfs*scale_coef)+(1*dfs*scale_coef))
        max_x = int((resyy*dfs*scale_coef)+(1*dfs*scale_coef))
        
        temp_data00 = []
        for data in range(0,max_y,dfs*scale_coef):
            for node in range(len(temp_data0)):
                val1 = temp_data0[node]
                if round(val1[1]) == data:
                    if round(val1[0]) % (dfs*scale_coef) == 0:
                        temp_data00.append(val1)
        
        
        temp_data1 = []
        for data in range(len(temp_data00)):
            value = temp_data00[data]
                    
            if value[1] == latitude*scale_coef:
                values = value[2]
                temp_data1.append(values)
                        
        print("bottom temperature from left to right is: {}".format(temp_data1))
        
     
    
base = int(resx+1)
nodes_all = []
for node in range(0,base,1):
 
    nodes_all.append(node)

print("nodes to assign temperature to are: {}".format(nodes_all))

print(temp_data1)

temp_at_nodes = []
for data in range(len(nodes_all)):
    temp = ([nodes_all[data]], temp_data1[data] * u.degK)
    temp_at_nodes.append(temp)


Model.set_temperatureBCs(top= 273.15 * u.degK,
                        bottom=1000 * u.degK,
                        materials=[(air, 273.15 * u.degK),
                                  (sair, 273.15 * u.degK)],
                                    nodeSets=(temp_at_nodes))

Model.set_heatFlowBCs(bottom=(-0.022 * u.watt / u.metre**2, ml))



Model.set_velocityBCs(left=[0.0 * u.centimeter / u.year, None],
                      right=[0.0 * u.centimeter / u.year, None],
                      top=[None, 0.0 * u.centimeter / u.year])



import pandas as pd
import numpy as np
import math

with open("traction_vector_uwg_36.0_Ma.txt", "r") as input_crust:
    with open("30.txt", "w") as output_crust:    
        crust = pd.read_csv(input_crust, dtype=float, delimiter=",",skiprows=0)
        df_crust = pd.DataFrame(crust)
        df_crust = np.array(df_crust)

        x = df_crust[:,0]
        y = df_crust[:,1]
        k1 = df_crust[:,2]
        k2 = df_crust[:,3]
        k3 = df_crust[:,4]
        f = len(x)
        q1 = x.reshape(f,1)
        q2 = y.reshape(f,1)
        q3 = k1.reshape(f,1)
        q4 = k2.reshape(f,1)
        q5 = k3.reshape(f,1)
        index1 = q1
                     
        tracton_data0 = []
        coor_traction = []
                
        for data in range(len(index1)):

            x_coor = x[data]*scale_coef
            y_coor = y[data]*scale_coef
            coor_traction = (float(x_coor),float(y_coor),float(k1[data]),float(k2[data]),float(k3[data]))
            tracton_data0.append(coor_traction)
    
# HERE WE DEFINE WHICH LATITUDE WE WANT TO CHOOSE FOR MAKING 2D PROFILES BY DEFINING LATITUDE VALUE.
# THE RESOLUTION FOR LATITUDES IS 0.1 DEGREE FROM 27N TO 44N.

        lat_int = lat_int_profile
        lat_decimal = lat_decimal_profile 
        latitude = ((lat_int - min_lat) * 10) + (lat_decimal * 10)
  
        max_y = int((resyy*dfs*scale_coef)+(1*dfs*scale_coef))
        max_x = int((resyy*dfs*scale_coef)+(1*dfs*scale_coef))
        
        tracton_data1 = []
        
        for data in range(0,max_y,dfs*scale_coef):
            for node in range(len(tracton_data0)):
                val1 = tracton_data0[node]
                if round(val1[1]) == data:
                    if round(val1[0]) % (dfs*scale_coef) == 0:
                        tracton_data1.append(val1)
                                
        values1 = []
        for data in range(len(tracton_data1)):
            value = tracton_data1[data]
            
            if value[1] == latitude*scale_coef:
                values = value[2],value[3],value[4]
                values1.append(values)

base = int(resx+1)
nodes_all = []
for node in range(0,base,1):
 
    nodes_all.append(node)

print("nodes to assign traction to are: {}".format(nodes_all))
    
depth_asthen = int(((0 + ((resy*(nz)) - (elsl-dcp)))/(nz)))
incr_asthen = int(depth_asthen+1)

t_correction = heat_add/incr_asthen

traction_at_nodes = []
for data in range(len(nodes_all)):
    temp = ([nodes_all[data]], [values1[data][0] * u.pascal, values1[data][2] * u.pascal])
    traction_at_nodes.append(temp)

print("bottom tractions from left to right are: {}".format(traction_at_nodes))
            

Model.set_stressBCs(bottom=[0., 0.], nodeSets=(traction_at_nodes))



Model.init_model()




# Set Courant number
GEO.rcParams["CFL"] = 0.5
# Set Solver
Model.solver.set_inner_method("lu")
# Set penalty value. default is 0.0 or None
Model.solver.set_penalty(1e6)
# Set maximal number of Picard iterations (first solve). default is 500.
GEO.rcParams["initial.nonlinear.max.iterations"] = 20
# Set maximal number of Picard iterations. default is 500.
GEO.rcParams["nonlinear.max.iterations"] = 15
# Set nonlinear tolerance for Stokes first solve. default is 1e-2.
# this defines when the convergence happens based on the value for Residual.
GEO.rcParams["initial.nonlinear.tolerance"]= 2e-2 
# Set nonlinear tolerance for solves. default is 1e-2.
# this defines when the convergence happens based on the value for Residual.
GEO.rcParams["nonlinear.tolerance"]= 2e-2
# output units for temperature field. default is u.degK
GEO.rcParams["temperature.SIunits"]= u.degC
# output units for velocity field. default is u.centimetre / u.year
GEO.rcParams["velocityField.SIunits"]= u.millimeter / u.year




Model.run_for(36000000.0 * u.years, checkpoint_interval=100000.0*u.years,dt=100000.0*u.years)


