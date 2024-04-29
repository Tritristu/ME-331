import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir

# Constants
# Squares are [length,width,height] and circles are [length,diameter]
dimensionsCuSq = 1e-2*np.array([28.5,1.27,1.27]) # m
dimensionsCuCirc = 1e-2*np.array([17.2,1.27]) # m
dimensionsFeCirc = 1e-2*np.array([28.5,0.95]) # m
dimensionsAlCirc = 1e-2*np.array([28.5,1.27]) # m

conductivityCu = 388 # W/(m*K)
conductivityFe = 16 # W/(m*K)
conductivityAl = 167 # W/(m*K)

tempBaseCuSq = 68 # Centigrade
tempBaseCuCirc = 81 # Centigrade
tempBaseFeCirc = 87 # Centigrade
tempBaseAlCirc = 91 # Centigrade

# Distance from the base
positionsCuSq = 1e-2*np.array([3.06,5.5,9.82,17.15,22.07,28.07]) # m
positionsCuCirc = 1e-2*np.array([3.09,5.56,9.8,17.15]) # m
positionsFeCirc = 1e-2*np.array([3.1,5.44,9.81,17.22,22.07,28.45]) # m
positionsAlCirc = 1e-2*np.array([4.35,5.41,9.72,17.15,22.07,28.07]) # m

# Loading in Data
Files = [x for x in listdir() if '.csv' in x]

Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv(File)

# Calculate the Steady-State Temperatures
AlCircTemps = []
CuCircTemps = []
CuSquareTemps = []
FeCircTemps = []

def tempAvgs(data,avgs):
    positions = data.keys()
    for position in positions:
        if (position != 'Time'):
            avgs.append(data[position].mean())

tempAvgs(Data[Files[0]],AlCircTemps)
tempAvgs(Data[Files[1]],CuCircTemps)
tempAvgs(Data[Files[2]],CuSquareTemps)
tempAvgs(Data[Files[3]],FeCircTemps)

# calculate convective coefficient using least squares error fitting
def normExcessTempAna(position,convectiveCoeff,conductivity,dimensions,circle):
    if circle:
        perimeter = np.pi*dimensions[1]
        areaCS = np.pi*(dimensions[1]/2)**2
    else:
        perimeter = 2*(dimensions[1] + dimensions[2])
        areaCS = dimensions[1]*dimensions[2]

    m = np.sqrt(convectiveCoeff*perimeter/(conductivity*areaCS))
    numerator = np.cosh(m*(dimensions[0]-position)) + (convectiveCoeff/(m*conductivity))*np.sinh(m*(dimensions[0]-position))
    denominator = np.cosh(m*dimensions[0]) + (convectiveCoeff/(m*conductivity))*np.sinh(m*dimensions[0])

    return numerator/denominator

def convectiveCalc(positions,exTemps,range,tempBase,conductivity,dimensions,circle):
    bestError = 1e9
    convectiveCoeff = -1
    baseExcessTemp = tempBase - exTemps[len(exTemps)-1]
    excessTemp = np.array(exTemps) - exTemps[len(exTemps)-1]
    normExcessTemp = excessTemp/baseExcessTemp

    for convection in np.arange(range[0],range[1],range[2]):
        errorLeastSquares = 0
        for i,position in enumerate(positions):
            error = (normExcessTemp[i] - normExcessTempAna(position,convection,conductivity,dimensions,circle))**2
            errorLeastSquares = errorLeastSquares + error
        if errorLeastSquares < bestError:
            bestError = errorLeastSquares
            convectiveCoeff = convection
    
    return convectiveCoeff

convectionCuSq = convectiveCalc(positionsCuSq,CuSquareTemps,[0.1,100,1],tempBaseCuSq,conductivityCu,dimensionsCuSq,circle=False)
convectionCuCirc = convectiveCalc(positionsCuCirc,CuCircTemps,[0.1,100,1],tempBaseCuCirc,conductivityCu,dimensionsCuCirc,circle=True)
convectionFeCirc = convectiveCalc(positionsFeCirc,FeCircTemps,[0.1,100,1],tempBaseFeCirc,conductivityFe,dimensionsFeCirc,circle=True)
convectionAlCirc = convectiveCalc(positionsAlCirc,AlCircTemps,[0.1,100,1],tempBaseAlCirc,conductivityAl,dimensionsAlCirc,circle=True)

print(convectionCuSq)
print(convectionCuCirc)
print(convectionFeCirc)
print(convectionAlCirc)

# Fin Heat Transfer Rates
def heatTransferRates(tempBase,exTemps,convectiveCoeff,conductivity,dimensions,circle):
    if circle:
        perimeter = np.pi*dimensions[1]
        areaCS = np.pi*(dimensions[1]/2)**2
    else:
        perimeter = 2*(dimensions[1] + dimensions[2])
        areaCS = dimensions[1]*dimensions[2]

    baseExcessTemp = tempBase - exTemps[len(exTemps)-1]
    m = np.sqrt(convectiveCoeff*perimeter/(conductivity*areaCS))
    M = baseExcessTemp*np.sqrt(convectiveCoeff*perimeter*conductivity*areaCS)
    numerator = np.sinh(m*dimensions[0]) + (convectiveCoeff/(m*conductivity))*np.sinh(m*dimensions[0])
    denominator = np.cosh(m*dimensions[0]) + (convectiveCoeff/(m*conductivity))*np.sinh(m*dimensions[0])

    return M*numerator/denominator

heatTransferCuSq = heatTransferRates(tempBaseCuSq,CuSquareTemps,convectionCuSq,conductivityCu,dimensionsCuSq,circle=False)
heatTransferCuCirc = heatTransferRates(tempBaseCuCirc,CuCircTemps,convectionCuCirc,conductivityCu,dimensionsCuCirc,circle=True)
heatTransferAlCirc = heatTransferRates(tempBaseAlCirc,AlCircTemps,convectionAlCirc,conductivityAl,dimensionsAlCirc,circle=True)
heatTransferFeCirc = heatTransferRates(tempBaseFeCirc,FeCircTemps,convectionFeCirc,conductivityFe,dimensionsFeCirc,circle=True)


plt.figure(1)
plt.plot(np.concatenate((np.array([0]),positionsCuSq)),np.concatenate((np.array([tempBaseCuSq]),CuSquareTemps[0:len(positionsCuSq)])),'b',label='Copper Square')
plt.plot(np.concatenate((np.array([0]),positionsCuCirc)),np.concatenate((np.array([tempBaseCuCirc]),CuCircTemps[0:len(positionsCuCirc)])),'r',label='Copper Round')
plt.plot(np.concatenate((np.array([0]),positionsAlCirc)),np.concatenate((np.array([tempBaseAlCirc]),AlCircTemps[0:len(positionsAlCirc)])),'g',label='Aluminum Round')
plt.plot(np.concatenate((np.array([0]),positionsFeCirc)),np.concatenate((np.array([tempBaseFeCirc]),FeCircTemps[0:len(positionsFeCirc)])),'m',label='Stainless Steel Round')
plt.xlim([0,dimensionsCuSq[0]])
plt.ylim([0,100])
plt.ylabel('Temperature (Centigrade)')
plt.xlabel('Position (m)')
plt.legend()

plt.show()