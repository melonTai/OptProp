#author --fujita yuki--
# -*- coding: utf-8 -*-
from nsga3_base import nsga3
from xrotor import Xrotor
from model import RotorModel
from scipy import interpolate
import sys,os
import numpy as np

#ディープ
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
"""
------------------------------------
設計諸元
------------------------------------
# rotor(tale)
    - diameter : 0.3 m
    - tip radius : 0.15 m
    - hub radius : 0.03 m
    - Thrust : 3Nf
    - rpm : 2000
    - density : 1.226e-2
# aerofoil
    - the_best_v4
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
    - GA
# individual
    - size : sn × 2
    ## content
        - 0 to 5  : chord[m] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        - 6 to 11 : beta[deg] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
# constraint
    - chord : 0.005m <= chords[i] <= 0.1m
    - beta  : 0deg <= betas[i] <= 90deg
# penalty
    - Thrust >= 3Nf
------------------------------------
プログラム概要
------------------------------------
設計諸元において効率が最大となるプロペラの形状モデルファイルを作成するプログラム
最適化は遺伝的アルゴリズムを用いている。
最適化のアルゴリズムはnsga3_base.pyから継承している。
各世代ごとの最適解をbestRotorgen[世代数].txtとして出力し、
各世代において暫定最適解複数候補をbestRotor[0~199].txtとして出力する。

プロペラの性能計算はcrotorを使用している。
http://www.esotec.org/sw/crotor.html

モデルファイルの形式は
http://www.esotec.org/sw/dl/Espara_doc.txt
のIMPO欄を参照

self.aerofに指定した翼型モデルファイルは
http://web.mit.edu/drela/Public/web/xrotor/xrotor_doc.txt
のAERO欄を参照
"""
class newnsga3(nsga3):
    def __init__(self):
        super().__init__()
        self.tipr = 0.15
        self.hubr = 0.01
        self.T1 = 3
        self.rpm1 = 2000
        self.r_R = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        self.sn = len(self.r_R)
        self.b = 2
        self.fs = 0.0
        self.velo = 0.1
        self.aerof = "the_best_v4_Re5000.txt"
        self.dens = 1.226e-2

        # individual
        ## size
        self.NDIM = 12
        ## constraint
        self.BOUND_LOW = [0.005] * 6 + [0] * 6
        self.BOUND_UP = [0.1] * 6 + [90] * 6

        self.weights = (-1.0,)
        self.NOBJ = len(self.weights)#評価関数の数
        self.MU = 200#人口の数
        self.NGEN = 300#世代数
        self.CXPB = 0.7#交叉の確立(1を100%とする)
        self.MUTPB = 0.5#突然変異の確立(1を100%とする)
        self.cx_eta = 20
        self.mut_eta = 20
        self.thread = 1

    def setup(self):
        super().setup()
        self.toolbox.register("select", tools.selTournament,  tournsize = 10)

    def beauty(self,x,y):
        yd = [(y[i+1] - y[i])/(x[i+1] - x[i]) for i in range(len(x)-1)]
        ydd = [(yd[i+1] - yd[i])/(x[i+1] - x[i]) for i in range(len(yd)-1)]
        beauty = sum([ydd[i]**2/((1 + yd[i]**2)**2.5)*(x[i + 1] - x[i]) for i in range(len(ydd))])
        return abs(beauty)

    def spline(self,x,y,num = 100):
        """
        #x:controll points of x axis
        #y:controll points of y axis
        """
        tck,u = interpolate.splprep([x,y], k=3, s=0)
        u = np.linspace(0.0, 1.0,num=num,endpoint=True)
        x, y = interpolate.splev(u,tck)
        return x, y

    def evaluate(self,individual):
        #プロセスid取得(並列処理用)
        id = os.getpid()

        #rotorモデル作成
        chords = individual[:int(len(individual)/2)]
        betas = individual[int(len(individual)/2):]
        radii = [a * self.tipr for a in self.r_R]
        rotorf = "rotor" + str(id)
        resultf = "res" + str(id)

        rm = RotorModel("dumrotor" + str(id), self.tipr, self.hubr, self.sn, radii, chords, betas)
        rm.writefile(rotorf)
        #xrotorコマンド設定
        xr = Xrotor(self.b, self.fs)
        xr.aero(self.aerof)
        xr.impo(rotorf)
        xr.dens = self.dens
        xr.velo = self.velo
        xr.rpm = self.rpm1
        xr.oper()
        xr.cput(resultf)
        res = xr.call(timeout = 7)
        penalty = 0
        try:
            if res == None:
                raise Exception("failed to complete xrotor")
            else:
                result = np.loadtxt(resultf, skiprows=3)
                eff = result[-1]
                T1 = result[10]
                #ペナルティ
                if T1 < self.T1:
                    penalty += (self.T1 - T1)
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
        #x, y = self.spline(radii,chords,100)
        #obj2 = self.beauty(x, y)
        #x, y = self.spline(radii,betas, 100)
        #obj3 = self.beauty(x, y)


        return (obj1,)

    def main(self,seed=None):
        self.setup()
        # Initialize statistics object
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("std", np.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        #初期化(個体生成のこと)
        pop = self.toolbox.population(n=self.MU)

        #進化の始まり
        # Begin the generational process
        for gen in range(self.NGEN):

            if(gen == 0):
                #0世代目の評価
                # Evaluate the individuals with an invalid fitness
                invalid_ind = [ind for ind in pop if not ind.fitness.valid]
                fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
                for ind, fit in zip(invalid_ind, fitnesses):
                    ind.fitness.values = fit

            else:
                offspring = algorithms.varAnd(pop, self.toolbox, self.CXPB, self.MUTPB)
                #評価
                # Evaluate the individuals with an invalid fitness
                invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
                fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
                for ind, fit in zip(invalid_ind, fitnesses):
                    ind.fitness.values = fit

                #淘汰
                # Select the next generation population from parents and offspring
                pop = self.toolbox.select(pop + offspring, self.MU)

            #評価
            pop_fit = np.array([ind.fitness.values for ind in pop])
            chords = pop[0][:int(len(pop[0])/2)]
            betas = pop[0][int(len(pop[0])/2):]
            radii = [a * self.tipr for a in self.r_R]
            rotorf = "bestRotorgen"+ str(gen)+".txt"
            rm = RotorModel("bestRotor", self.tipr, self.hubr, self.sn, radii, chords, betas)
            rm.writefile(rotorf)
            record = stats.compile(pop)
            k = 0
            for ind in pop:
                k += 1
                chords = ind[:int(len(ind)/2)]
                betas = ind[int(len(ind)/2):]
                radii = [a * self.tipr for a in self.r_R]
                rotorf = "bestRotor" + str(k) + ".txt"
                rm = RotorModel("bestRotor", self.tipr, self.hubr, self.sn, radii, chords, betas)
                rm.writefile(rotorf)
            # Compile statistics about the new population
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            print(logbook.stream)

        return pop, logbook

if __name__ == "__main__":

    ng = newnsga3()
    pop, stats = ng.main()
    pop_fit = np.array([ind.fitness.values for ind in pop])
    try:
        k = 0
        for ind in pop:
            k += 1
            chords = ind[:int(len(ind)/2)]
            betas = ind[int(len(ind)/2):]
            radii = [a * ng.tipr for a in ng.r_R]
            rotorf = "bestRotor" + str(k) + ".txt"
            rm = RotorModel("bestRotor", ng.tipr, ng.hubr, ng.sn, radii, chords, betas)
            rm.writefile(rotorf)
        k = 0
        for ind in pop_fit:
            k+=1
            print("fit" + str(k) + ":" + str(ind))
    except Exception as e:
        print("message:{0}".format(e))
