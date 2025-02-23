import numpy as np
from matplotlib import pyplot as plt


# constants
heater_power = 500  # heater power in watts
heater_diameter = 0.375 * 25.4 / 1000
heater_length = 6 * 25.4 / 1000
heater_flux = heater_power / (np.pi * heater_diameter * heater_length)
k_conductivity = 16.2   # thermal conductivity in W/m/K
rho_density = 7930      # density in kg/m^3
cp_heatCapacity = 500   # heat capacity in J/kg/K
alpha_diffusivity = k_conductivity / rho_density / cp_heatCapacity # thermal diffusivity in m^2/s

dt_timeStep = 0.001      # time step in seconds
dx_cellWidth = 0.001     # cell width in meters
numCells = 65         # number of cells, should be odd
x_start = 0
x_stop = x_start + numCells*dx_cellWidth

# pre-allocate arrays
#x_cellValues = np.arrange(x_start, x_stop, dx_cellWidth)
x_cellValues = np.linspace(x_start, x_stop, numCells)

temperature_initial = 300   # initial temperature in K
T1 = 300                    # temperature at left boundary representing a heat pipe
T2 = 300                    # temperature at right boundary representing a heat pipe
temperature_old = np.ones(numCells) * temperature_initial  # 'old' temperature to be used in loop
temperature_new = np.ones(numCells) * temperature_initial  # 'new' temperature to be used in loop

qdot_heatGenerationRate = np.zeros(numCells)               # volumetric heat generation rate array, will only have one non-zero value
qdot_heatGenerationRate[int(np.ceil(numCells/2))] = heater_flux / dx_cellWidth   # heater power converted to volumetric heat generation rate


# loop constants
tol = 1e-6
err = 100
num_iterations = 0
max_iterations = 1000000


while err > tol and num_iterations < max_iterations:
    for cell in range(numCells):
        if cell == 0:
            temperature_new[cell] = alpha_diffusivity*dt_timeStep/k_conductivity*qdot_heatGenerationRate[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*temperature_old[cell+1] + (1-(2*alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)))*temperature_old[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*T1
        elif cell == numCells-1:
            temperature_new[cell] = alpha_diffusivity*dt_timeStep/k_conductivity*qdot_heatGenerationRate[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*T2 + (1-(2*alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)))*temperature_old[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*temperature_old[cell-1]
        else:
            temperature_new[cell] = alpha_diffusivity*dt_timeStep/k_conductivity*qdot_heatGenerationRate[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*temperature_old[cell+1] + (1-(2*alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)))*temperature_old[cell] + alpha_diffusivity*dt_timeStep/(dx_cellWidth**2)*temperature_old[cell-1]
    err = np.max(np.abs(temperature_new - temperature_old))
#    print(err)
    temperature_old[:] = temperature_new[:]
    num_iterations += 1
#    print(num_iterations)
    if num_iterations >= max_iterations:
        print("Maximum iterations reached: ",num_iterations)

# analytical solution
x = np.linspace(x_start, x_stop/2, numCells)
T = heater_flux / k_conductivity * x + T1


# plot
plt.figure(facecolor='w',edgecolor='k',dpi=300,figsize=(5.75,4.3125))
plt.xlabel(r'Position (m)')
plt.ylabel(r'Temperature (K)')
plt.grid(visible=True,which='major',axis='both')
plt.plot(x_cellValues, temperature_new, '-k',label='numerical')
plt.plot(x,T,'--b',label='analytical')
plt.tight_layout()
plt.savefig('plots/1d_conduction_analytical_plot.pdf',transparent=True)
