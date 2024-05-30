import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from scipy.stats import linregress

# Constants
g = 9.81 # m/s^2

densityCuBlk = 9.08e3 # kg/m^3
densityCuGld = 9.17e3 # kg/m^3
densityRubber = 0.90e3 # kg/m^3

volumeCuBlk = 62.86e-6 # m^3
volumeCuGld = 62.21e-6 # m^3
volumeRubber = 24.43e-6 # m^3

diameterCuBlk = 49.34e-3 # m
diameterCuGld = 49.17e-3 # m
diameterRubber = 36e-3 # m

massCuBlk = 570.7e-3 # kg
massCuGld = 570.3e-3 # kg
massRubber = 22e-3 # kg

specificHeatCuBlk = 0.39e3 # J/kg*K
specificHeatCuGld = 0.39e3 # J/kg*K
specificHeatRubber = 2.01e3 # J/kg*K

thermalConductivityCuBlk = 401 # W/m*K
thermalConductivityCuGld = 401 # W/m*K
thermalConductivityRubber = 0.16 # W/m*K

emissitivity = 0.03

sigma = 5.670374419e-8 # W/(m^2 * K^4)

# Loading in Data
Files = [x for x in listdir('Lab 2') if '.csv' in x]

Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv(f'Lab 2\\{File}')

# Index
# 0 = CopperBlack
# 1 = CopperGold
# 2 = Rubber

# Plotting cooling response
plt.figure(1)
plt.plot(Data[Files[0]]['time'],Data[Files[0]]['specimen'],'b')
plt.xlim([0,10000])
plt.ylim([0,100])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')

plt.figure(2)
plt.plot(Data[Files[1]]['time'],Data[Files[1]]['specimen'],'b')
plt.xlim([0,10000])
plt.ylim([0,100])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')

plt.figure(3)
plt.plot(Data[Files[2]]['time'],Data[Files[2]]['specimen'],'b')
plt.xlim([0,10000])
plt.ylim([0,100])
plt.xlabel('Time (s)')
plt.ylabel('Temperature (Centigrade)')

# Calculating excess temps
for File in Files:
    Data[File]['log(excessTemp)'] = np.log((Data[File]['specimen']-Data[File]['air'])/(Data[File]['water']-Data[File]['air']))

# Linear regression
slopeCuBlk, interceptCuBlk, rCuBlk, p, se = linregress(Data[Files[0]]['time'], Data[Files[0]]['log(excessTemp)'])
slopeCuGld, interceptCuGld, rCuGld, p, se = linregress(Data[Files[1]]['time'], Data[Files[1]]['log(excessTemp)'])
slopeRubber, interceptRubber, rRubber, p, se = linregress(Data[Files[2]]['time'], Data[Files[2]]['log(excessTemp)'])


# Plotting logged excess temp & fitted lines
times = np.linspace(0,10000,10)
interp = lambda slope,intercept,time : slope*time + intercept

plt.figure(4)
plt.scatter(Data[Files[0]]['time'],Data[Files[0]]['log(excessTemp)'],color='b',s=0.5*np.ones(len(Data[Files[0]]['time'])))
plt.plot(times,interp(slopeCuBlk,interceptCuBlk,times),'r')
plt.xlim([0,10000])
plt.ylim([-5,1])
plt.xlabel('Time (s)')
plt.ylabel(r'ln($\frac{\theta}{\theta_i}$)')

plt.figure(5)
plt.scatter(Data[Files[1]]['time'],Data[Files[1]]['log(excessTemp)'],color='b',s=0.5*np.ones(len(Data[Files[1]]['time'])))
plt.plot(times,interp(slopeCuGld,interceptCuGld,times),'r')
plt.xlim([0,10000])
plt.ylim([-5,1])
plt.xlabel('Time (s)')
plt.ylabel(r'ln($\frac{\theta}{\theta_i}$)')

plt.figure(6)
plt.scatter(Data[Files[2]]['time'],Data[Files[2]]['log(excessTemp)'],color='b',s=0.5*np.ones(len(Data[Files[2]]['time'])))
plt.plot(times,interp(slopeRubber,interceptRubber,times),'r')
plt.xlim([0,10000])
plt.ylim([-5,1])
plt.xlabel('Time (s)')
plt.ylabel(r'ln($\frac{\theta}{\theta_i}$)')


# slope and r^2 values
# print(slopeCuBlk)
# print(slopeCuGld)
# print(slopeRubber)

# print(rCuBlk**2)
# print(rCuGld**2)
# print(rRubber**2)

# Calculate sphere surface area
surfaceAreaCuBlk = np.pi*(diameterCuBlk)**2
surfaceAreaCuGld = np.pi*(diameterCuGld)**2
surfaceAreaRubber = np.pi*(diameterRubber)**2


# Calculate h_eff
convectionEffCuBlk = -(slopeCuBlk*densityCuBlk*volumeCuBlk*specificHeatCuBlk)/(surfaceAreaCuBlk)
convectionEffCuGld = -(slopeCuGld*densityCuGld*volumeCuGld*specificHeatCuGld)/(surfaceAreaCuGld)
convectionEffRubber = -(slopeRubber*densityRubber*volumeRubber*specificHeatRubber)/(surfaceAreaRubber)

# print(convectionEffCuBlk)
# print(convectionEffCuGld)
# print(convectionEffRubber)

# Calculate characteristic length
characteristicLengthCuBlk = ((4/3)*np.pi*(diameterCuBlk/2)**3)/surfaceAreaCuBlk
characteristicLengthCuGld = ((4/3)*np.pi*(diameterCuGld/2)**3)/surfaceAreaCuGld
characteristicLengthRubber = ((4/3)*np.pi*(diameterRubber/2)**3)/surfaceAreaRubber

# Calculate Biot Number
BiotCuBlk = convectionEffCuBlk*characteristicLengthCuBlk/thermalConductivityCuBlk
BiotCuGld = convectionEffCuGld*characteristicLengthCuGld/thermalConductivityCuGld
BiotRubber = convectionEffRubber*characteristicLengthRubber/thermalConductivityRubber

# print(BiotCuBlk)
# print(BiotCuGld)
# print(BiotRubber)

# Calculate film temperature
Data[Files[0]]['Film Temperature'] = (Data[Files[0]]['specimen']+Data[Files[0]]['air'])/2 + 273.15
Data[Files[1]]['Film Temperature'] = (Data[Files[1]]['specimen']+Data[Files[1]]['air'])/2 + 273.15

# Linearly itnerpolate the kinematic viscosity & thermal conductivity
Data[Files[0]]['kinematic viscosity'] = (((20.92-15.89)/50)*(Data[Files[0]]['Film Temperature']-300) + 15.89)*1e-6 # m^2/s
Data[Files[1]]['kinematic viscosity'] = (((20.92-15.89)/50)*(Data[Files[1]]['Film Temperature']-300) + 15.89)*1e-6 # m^2/s

Data[Files[0]]['thermal conductivity'] = (((30-26.3)/50)*(Data[Files[0]]['Film Temperature']-300) + 26.3)*1e-3 # m^2/s
Data[Files[1]]['thermal conductivity'] = (((30-26.3)/50)*(Data[Files[1]]['Film Temperature']-300) + 26.3)*1e-3 # m^2/s


# Calculate Beta
Data[Files[0]]['Beta'] = 1/Data[Files[0]]['Film Temperature']
Data[Files[1]]['Beta'] = 1/Data[Files[1]]['Film Temperature']

# Calculate Grashof Number
Data[Files[0]]['Grashof Number'] = g*Data[Files[0]]['Beta']*(Data[Files[0]]['specimen']-Data[Files[0]]['air'])*(diameterCuBlk**3)/(Data[Files[0]]['kinematic viscosity']**2)
Data[Files[1]]['Grashof Number'] = g*Data[Files[1]]['Beta']*(Data[Files[1]]['specimen']-Data[Files[1]]['air'])*(diameterCuGld**3)/(Data[Files[1]]['kinematic viscosity']**2)

# Calculate Prandtl Number
Data[Files[0]]['Prandtl Number'] = -0.00017*Data[Files[0]]['Film Temperature']+0.758
Data[Files[1]]['Prandtl Number'] = -0.00017*Data[Files[1]]['Film Temperature']+0.758

# Calculate Rayleigh Number
Data[Files[0]]['Rayleigh Number'] = Data[Files[0]]['Grashof Number']*Data[Files[0]]['Prandtl Number']
Data[Files[1]]['Rayleigh Number'] = Data[Files[1]]['Grashof Number']*Data[Files[1]]['Prandtl Number']

# Calculate Nusselt Number
Data[Files[0]]['Nusselt Number'] = 2 + (0.589 * Data[Files[0]]['Rayleigh Number']**(1/4))/((1 + (0.469/(Data[Files[0]]['Prandtl Number']**(9/16))))**(4/9))
Data[Files[1]]['Nusselt Number'] = 2 + (0.589 * Data[Files[1]]['Rayleigh Number']**(1/4))/((1 + (0.469/(Data[Files[0]]['Prandtl Number']**(9/16))))**(4/9))

# Calculate convective bar coefficient
Data[Files[0]]['Convective coefficient'] = Data[Files[0]]['Nusselt Number']*Data[Files[0]]['thermal conductivity']/diameterCuBlk
Data[Files[1]]['Convective coefficient'] = Data[Files[0]]['Nusselt Number']*Data[Files[1]]['thermal conductivity']/diameterCuGld

# Calculate convective heat transfer
Data[Files[0]]['Convective HT'] = Data[Files[0]]['Convective coefficient']*surfaceAreaCuBlk*(Data[Files[0]]['specimen']-Data[Files[0]]['air'])
Data[Files[1]]['Convective HT'] = Data[Files[1]]['Convective coefficient']*surfaceAreaCuBlk*(Data[Files[1]]['specimen']-Data[Files[1]]['air'])

# Calculate radiative heat transfer
Data[Files[0]]['Radiative HT'] = emissitivity*sigma*(Data[Files[0]]['specimen']**4 - Data[Files[0]]['air']**4)
Data[Files[1]]['Radiative HT'] = emissitivity*sigma*(Data[Files[1]]['specimen']**4 - Data[Files[1]]['air']**4)

# Calculate total heat transfer using sum of components
Data[Files[0]]['Sum HT'] = Data[Files[0]]['Convective HT'] + Data[Files[0]]['Radiative HT']
Data[Files[1]]['Sum HT'] = Data[Files[1]]['Convective HT'] + Data[Files[1]]['Radiative HT']

# Calculate total heat transfer using effective convective coefficient
Data[Files[0]]['Effective HT'] = convectionEffCuBlk*surfaceAreaCuBlk*(Data[Files[0]]['specimen']-Data[Files[0]]['air'])
Data[Files[1]]['Effective HT'] = convectionEffCuGld*surfaceAreaCuGld*(Data[Files[1]]['specimen']-Data[Files[1]]['air'])

plt.figure(7)
plt.plot(Data[Files[0]]['Film Temperature'],Data[Files[0]]['Convective coefficient'],'b')
plt.xlim([295,335])
plt.ylim([0,10])
plt.xlabel('Film Temperature (K)')
plt.ylabel(r'Convective Coefficient $\left(\frac{W}{m^2\dot K}\right)$')

plt.figure(8)
plt.plot(Data[Files[1]]['Film Temperature'],Data[Files[1]]['Convective coefficient'],'b')
plt.xlim([295,335])
plt.ylim([0,10])
plt.xlabel('Film Temperature (K)')
plt.ylabel(r'Convective Coefficient $\left(\frac{W}{m^2\dot K}\right)$')


plt.figure(9)
plt.plot(Data[Files[0]]['Film Temperature'],Data[Files[0]]['Convective HT'],'b')
plt.plot(Data[Files[0]]['Film Temperature'],Data[Files[0]]['Radiative HT'],'r')
plt.plot(Data[Files[0]]['Film Temperature'],Data[Files[0]]['Sum HT'],'m')
plt.plot(Data[Files[0]]['Film Temperature'],Data[Files[0]]['Effective HT'],'g')
plt.xlim([295,335])
plt.ylim([0,5])
plt.xlabel('Film Temperature (K)')
plt.ylabel(r'Heat transfer rate $\left(\frac{W}{m^2}\right)$')

plt.figure(10)
plt.plot(Data[Files[1]]['Film Temperature'],Data[Files[1]]['Convective HT'],'b')
plt.plot(Data[Files[1]]['Film Temperature'],Data[Files[1]]['Radiative HT'],'r')
plt.plot(Data[Files[1]]['Film Temperature'],Data[Files[1]]['Sum HT'],'m')
plt.plot(Data[Files[1]]['Film Temperature'],Data[Files[1]]['Effective HT'],'g')
plt.xlim([295,335])
plt.ylim([0,5])
plt.xlabel('Film Temperature (K)')
plt.ylabel(r'Heat transfer rate $\left(\frac{W}{m^2}\right)$')

plt.show()