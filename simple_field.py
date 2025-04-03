import numpy as np 

def generic_field_3(coef, x):
    return coef[0] + coef[1]*x[:,0] + coef[2]*x[:,2] + coef[3]*(2*x[:,0]*x[:,2]) + coef[4]*(x[:,0]**2-x[:,2]**2) + coef[5]*(3*x[:,2]*(x[:,0]**2)-x[:,2]**3) + coef[6]*(x[:,0]**3 - 3*x[:,0]*(x[:,2]**2)) + coef[7](x[:,0]**4 - 6*((x[:,0]*x[:,2])**2)) + coef[8]*(4*(x[:,2]*x[:,0]**3-x[:,0]*x[:,2]**3))

def generic_xderiv_3(coef, x):
    return coef[1] + coef[3]*(2*x[:,2]) + coef[4]*(2*x[:,0]) + coef[5]*(6*x[:,2]*x[:,0]) + coef[6]*(3*x[:,0]**2 - 3*x[:,2]**2) + coef[7]*(4*x[:,0]**3 - 12*x[:,0]*x[:,2]**2) + coef[8]*(4*(3*x[:,2]*x[:,0]**2)-x[:,2]**3)

def generic_zderiv_3(coef, x):
    return coef[2] + coef[3]*(2*x[:,0]) + coef[4]*(-2*x[:,2]) + coef[5]*(3*x[:,0]**2 - 3*x[:,2]**2) + coef[6]*(-6*x[:,0]*x[:,2]) + coef[7]*(4*x[:,2]**3 - 12*(x[:,0]**2)*x[:,2]) + coef[8]*(4*(x[:,0]**3 - 3*x[:,0]*x[:,2]**2))