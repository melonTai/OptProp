#author --fujita yuki--
# -*- coding: utf-8 -*-
"""
main.pyの並列処理用プログラム
"""
from xrotor import Xrotor
from model import RotorModel
from scipy import interpolate
import sys,os
import numpy as np
import random

#ディープ
from deap import algorithms
from deap import base
from deap import creator
from deap import tools

#並列処理
from multiprocessing import Pool,Value

#=====================================================
#最適化の定義
#=====================================================
tipr = 0.07
hubr = 0.005
oT1 = 1
oT2 = 3
rpm1 = 6500
rpm2 = 9500
r_R = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
sn = len(r_R)
b = 2
fs = 0.01
velo = 0.01
aerof = "AG14_Re50000.txt"
NOBJ = 1#評価関数の数
NDIM = 12
P = 12
BOUND_LOW = [0.005] * 6 + [0] * 6
BOUND_UP = [0.04] * 6 + [22] * 6
weights = (-1.0,)
MU = 200#人口の数
NGEN = 300#世代数
CXPB = 1.0#交叉の確立(1を100%とする)
MUTPB = 0.7#突然変異の確立(1を100%とする)

# Create uniform reference point
ref_points = tools.uniform_reference_points(NOBJ, P)

# Create classes
creator.create("FitnessMin", base.Fitness, weights=weights)
creator.create("Individual", list, fitness=creator.FitnessMin)
##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#(いじるのここから
#=====================================================
#評価関数の定義
#=====================================================
def evaluate(individual):
    global tipr, hubr, oT1, oT2, rpm1, rpm2, r_R, sn, b, fs, velo, aerof
    #プロセスid取得(並列処理用)
    id = os.getpid()

    #rotorモデル作成
    chords = individual[:int(len(individual)/2)]
    betas = individual[int(len(individual)/2):]
    radii = [a * tipr for a in r_R]
    rotorf = "rotor" + str(id)
    resultf = "res" + str(id)

    rm = RotorModel("dumrotor" + str(id), tipr, hubr, sn, radii, chords, betas)
    rm.writefile(rotorf)
    #xrotorコマンド設定
    xr = Xrotor(b, fs)
    xr.aero(aerof)
    xr.impo(rotorf)
    xr.velo = velo
    xr.rpm = rpm1
    xr.oper()
    xr.rpm = rpm2
    xr.oper()
    xr.cput(resultf)
    res = xr.call(timeout = 7)
    penalty = 0
    try:
        if res == None:
            raise Exception("failed to complete xrotor")
        else:
            result = np.loadtxt(resultf, skiprows=3)
            eff = result[0][-1]
            T1 = result[0][10]
            T2 = result[1][10]
            #ペナルティ
            if T1 < oT1:
                penalty += (oT1 - T1)
            if T2 < oT2:
                penalty += (oT2 - T2)
    except Exception as e:
        print(e)
        eff = 0
        T = -1
        penalty += 10

    try:
        os.remove(rotorf)
    except Exception as e:
        print(e)
    try:
        os.remove(resultf)
    except Exception as e:
        print(e)

    #目的値
    obj1 = -eff + penalty
    return (obj1,)
#いじるのここまで)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Toolbox initialization
def uniform(low, up, size=None):
    try:
        return [random.uniform(a, b) for a, b in zip(low, up)]
    except TypeError:
        return [random.uniform(a, b) for a, b in zip([low] * size, [up] * size)]
##
toolbox = base.Toolbox()
toolbox.register("attr_float", uniform, BOUND_LOW, BOUND_UP, NDIM)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate",evaluate)
toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=BOUND_LOW, up=BOUND_UP, eta=10.0)
toolbox.register("mutate", tools.mutPolynomialBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0, indpb=1/NDIM)
toolbox.register("select", tools.selTournament,  tournsize = 10)


#=====================================================
#最適化アルゴリズム本体
#=====================================================
def main():
    global CXPB, MUTPB, MU, NGEN, tipr, hubr, sn, r_R, toolbox

    #同時並列数(空白にすると最大数になる)
    pool = Pool(4)
    toolbox.register("map", pool.map)
    # Initialize statistics object
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    logbook = tools.Logbook()
    logbook.header = "gen", "evals", "std", "min", "avg", "max"

    #初期化(個体生成のこと)
    pop = toolbox.population(n=MU)

    #進化の始まり
    # Begin the generational process
    for gen in range(NGEN):

        if(gen == 0):
            #0世代目の評価
            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in pop if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

        else:
            offspring = algorithms.varAnd(pop, toolbox, CXPB, MUTPB)
            #評価
            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            #淘汰
            # Select the next generation population from parents and offspring
            pop = toolbox.select(pop + offspring, MU)

        #評価
        pop_fit = np.array([ind.fitness.values for ind in pop])
        chords = pop[0][:int(len(pop[0])/2)]
        betas = pop[0][int(len(pop[0])/2):]
        radii = [a * tipr for a in r_R]
        rotorf = "bestRotorgen"+ str(gen)+".txt"
        rm = RotorModel("bestRotor", tipr, hubr, sn, radii, chords, betas)
        rm.writefile(rotorf)
        record = stats.compile(pop)
        k = 0
        for ind in pop:
            k += 1
            chords = ind[:int(len(ind)/2)]
            betas = ind[int(len(ind)/2):]
            radii = [a * tipr for a in r_R]
            rotorf = "bestRotor" + str(k) + ".txt"
            rm = RotorModel("bestRotor", tipr, hubr, sn, radii, chords, betas)
            rm.writefile(rotorf)
        # Compile statistics about the new population
        logbook.record(gen=gen, evals=len(invalid_ind), **record)
        print(logbook.stream)

    return pop, logbook


if __name__ == "__main__":
    pop, stats = main()
    pop_fit = np.array([ind.fitness.values for ind in pop])
    try:
        k = 0
        for ind in pop:
            k += 1
            chords = ind[:int(len(ind)/2)]
            betas = ind[int(len(ind)/2):]
            radii = [a * tipr for a in r_R]
            rotorf = "bestRotor" + str(k) + ".txt"
            rm = RotorModel("bestRotor", tipr, hubr, sn, radii, chords, betas)
            rm.writefile(rotorf)
        k = 0
        for ind in pop_fit:
            k+=1
            print("fit" + str(k) + ":" + str(ind))
    except Exception as e:
        print("message:{0}".format(e))
