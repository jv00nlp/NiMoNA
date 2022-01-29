from Const import *
from numpy import exp


def none(ACompartment, ATimeStep):
    return None


def R0_mitigating(ACompartment, ATimeStep, r0=3.35):
    # Why using the incidence?
    # TODO: Change to percentage of r0
    if (ATimeStep > 30) and (ATimeStep < 60):
        R0 = 0.5 * (ACompartment[0][0] / (ACompartment[0][0] + ACompartment[0][1] + ACompartment[0][2]))
    elif (ATimeStep > 240) and (ATimeStep < 300):
        R0 = 0.9 * (ACompartment[0][0] / (ACompartment[0][0] + ACompartment[0][1] + ACompartment[0][2]))
    else:
        R0 = r0 * (ACompartment[0][0] / (ACompartment[0][0] + ACompartment[0][1] + ACompartment[0][2]))

    return round(R0 * 0.0129, 5)
