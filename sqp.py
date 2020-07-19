from scipy.optimize import minimize,shgo
import numpy as np
from xrotor import Xrotor
from model import RotorModel
import os
"""
------------------------------------
設計諸元
------------------------------------
# rotor(tale)
    - diameter : 0.13 m
    - tip radius : 0.065 m
    - hub radius : 0.005 m
    - Thrust : 2Nf
    - rpm : 7000
------------------------------------
最適化の定義
------------------------------------
# definition
    - R : tip radius
    - r : position @ rotor from its root
    - chord : 翼弦長
    - beta : 取付角
    - r_R : r/R [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    - sn : section number
# variables
    - chord[m] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    - beta[deg] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
# objective
    - minimize efficiency
# optimization method
    - SQP
# variables
    - size : sn × 2
    ## content
        - 0 to 5  : chord[m] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        - 6 to 11 : beta[deg] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
# constraint
    - chord : 0.005mm <= chords[i] <= 0.05
    - beta  : 0deg <= betas[i] <= 90deg
# penalty
    - Thrust >= 2Nf
"""


tipr = 0.065
hubr = 0.005
rpm = 10000
velo = 0.1
r_R = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
sn = len(r_R)
aerof = "AG14_Re50000.txt"
T = 3.0

def function(x):
    global r_R, tipr, hubr, rpm, velo, sn, aerof
    chords = x[:int(len(x)/2)]
    betas = x[int(len(x)/2):]
    radii = [a * tipr for a in r_R]
    rotorf = "rotordum1.txt"
    resultf = "resdum1.txt"

    rm = RotorModel("dumrotor", tipr, hubr, sn, radii, chords, betas)
    rm.writefile(rotorf)
    #xrotorコマンド設定
    xr = Xrotor(2, 0.1)
    xr.aero(aerof)
    xr.impo(rotorf)
    xr.velo = velo
    xr.rpm = rpm
    xr.oper(resultf)
    res = xr.call(timeout = 7)
    try:
        if res == None:
            raise Exception("failed to complete xrotor")
        else:
            result = np.loadtxt(resultf, skiprows=3)
            if len(result) != 13:
                eff = 0
                t = 0
            else:
                eff = result[-1]
                t = result[10]
    except Exception as e:
        print(e)
        eff = 0
        t = 0

    try:
        os.remove(rotorf)
    except Exception as e:
        print(e)
    try:
        os.remove(resultf)
    except Exception as e:
        print(e)

    return eff, t

fun = lambda x : -function(x)[0]
cons = ({'type': 'ineq', 'fun': lambda x: function(x)[1]-T})
bnds = ((0.005,0.06),)*sn + ((0,50),)*sn

init = (2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02,
       2.0000e-02, 4.3204e+01, 4.2122e+01, 3.8006e+01, 2.9634e+01,
       2.3030e+01, 2.0777e+01)
#3.0828Nf(2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02,2.0000e-02, 7.3204e+01, 5.2122e+01, 3.8006e+01, 2.9634e+01,2.3030e+01, 2.0777e+01)
#汚い(0.035, 0.036, 0.043, 0.024, 0.022, 0.035, 32.931, 36.454, 19.340, 35.691, 8.890, 52.345)
#(0.009,0.012,0.0117,0.00964,0.0065,0.006,73.204, 52.122, 38.006, 29.634, 23.030, 20.777)
#(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 7.32040000e+01, 5.21220000e+01,3.80060000e+01, 2.96340000e+01, 2.30300000e+01, 2.07770000e+01)

#res = shgo(fun,bounds=bnds,constraints=cons,options = {'disp': True} )
res = minimize(fun, init, method='SLSQP', bounds=bnds,constraints=cons,options = {'ftol': 1e-9, 'disp': True})
print(res.x,res.fun)
#print(res.x)
