#!python3.6
import random
from math import factorial

import numpy as np
#import pymop.factory

#ディープ
from deap import algorithms
from deap import base
#from deap.benchmarks.tools import igd
from deap import creator
from deap import tools

#プロット
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # ３Dグラフ作成のため

#id取得
import os
import sys
import traceback

#並列処理
from multiprocessing.dummy import Pool,Value

class nsga3(object):
    """
    遺伝的アルゴリズムのテンプレートクラス
    デフォルトでnsga3のアルゴリズムが定義されている。
    各メソッドについて継承、オーバーライドして使いまわすことを目的としている。
    nsga3のアルゴリズムはdeapライブラリを用いてプログラム化している。
    https://deap.readthedocs.io/en/master/examples/nsga3.html


    # attributes
        - NOBJ(int)
            評価関数の数
        - NDIM(int)
            遺伝子の数
        - BOUD_LOW(float)
            遺伝子の下限
        - BOUD_UP(float)
            遺伝子の上限
        - MU(int)
            人口の数
        - NGEN(int)
            世代数
        - CXPB(float)
            交叉の確立
        - MUTPB(float)
            突然変異の確立
        - cx_eta(float)
            交叉の学習率
            交叉にcxSimulatedBinaryBoundedを指定したときに使用される
        - mut_eta(float)
            突然変異の学習率
            突然変異にmutPolynomialBoundedを指定したときに使用される
        - thread(int)
            並列処理数
        - weights(tuple)
            (-1.0 or 1.0,)*NOBJ
            評価関数(evaluateメソッドの戻り値)について、
            最大化する関数は1.0、
            最小化する関数は-1.0を指定する。

            例えば、
            evaluateの戻り値が(obj1,obj2,obj3)であるとき、
            obj1、obj3を最大化、obj2を最小化する場合は
            weights = (1.0, -1.0, 1.0)
    # method
        - setup
            deapライブラリのtoolboxに必要な関数を登録する
            詳細は以下リンクを確認
            https://deap.readthedocs.io/en/master/examples/nsga3.html
        - evaluate
            評価関数
            - argument
                - individual(float list)
            - return(array like)
                サイズNOBJのリストやタプルなど
                評価値
        - main
            遺伝的アルゴリズムを実行するメソッド

    """
    def __init__(self):
        #=====================================================
        #最適化の定義
        #=====================================================
        self.NOBJ = 3#評価関数の数
        self.NDIM = 12
        self.BOUND_LOW, self.BOUND_UP = -5.0, 5.0#遺伝子定義域
        self.MU = 100#人口の数
        self.NGEN = 200#世代数
        self.CXPB = 1.0#交叉の確立(1を100%とする)
        self.MUTPB = 0.7#突然変異の確立(1を100%とする)
        self.cx_eta = 10
        self.mut_eta = 20
        self.thread = 4
        self.weights = (-1.0)*self.NOBJ

    def setup(self):
        # Create uniform reference point
        self.ref_points = tools.uniform_reference_points(self.NOBJ, self.P)

        # Create classes
        creator.create("FitnessMin", base.Fitness, weights=self.weights)
        creator.create("Individual", list, fitness=creator.FitnessMin)
        ##

        self.toolbox = base.Toolbox()
        self.toolbox.register("attr_float", self.uniform, self.BOUND_LOW, self.BOUND_UP, self.NDIM)
        self.toolbox.register("individual", tools.initIterate, creator.Individual, self.toolbox.attr_float)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

        self.toolbox.register("evaluate",self.evaluate)
        self.toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=self.BOUND_LOW, up=self.BOUND_UP, eta=self.cx_eta)
        self.toolbox.register("mutate", tools.mutPolynomialBounded, low=self.BOUND_LOW, up=self.BOUND_UP, eta=self.mut_eta, indpb=1/self.NDIM)
        self.toolbox.register("select", tools.selNSGA3, ref_points=self.ref_points)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #(いじるのここから
    #=====================================================
    #評価関数の定義
    #=====================================================
    def evaluate(self,individual):
        pass
    # Toolbox initialization
    def uniform(self,low, up, size=None):
        try:
            return [random.uniform(a, b) for a, b in zip(low, up)]
        except TypeError:
            return [random.uniform(a, b) for a, b in zip([low] * size, [up] * size)]
    ##


    #=====================================================
    #最適化アルゴリズム本体
    #=====================================================
    def main(self,seed=None):
        self.setup()
        #同時並列数(空白にすると最大数になる)
        self.pool = Pool(self.thread)
        self.toolbox.register("map", self.pool.map)

        random.seed(seed)
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

                #0世代目の統計
                # Compile statistics about the population
                record = stats.compile(pop)
                logbook.record(gen=0, evals=len(invalid_ind), **record)

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

            record = stats.compile(pop)
            # Compile statistics about the new population
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            print(logbook.stream)

            return pop, logbook

if __name__ == "__main__":
    ng = nsga3()
    #翼型最適化開始
    try:
        pop, stats = ng.main()
    except KeyboardInterrupt:
        print("Ctrl + Cで停止しました")
        pass

    #最終世代の評価値取得
    pop_fit = np.array([ind.fitness.values for ind in pop])
