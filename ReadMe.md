------------------------------------
# プログラム概要
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


------------------------------------
# 設計諸元
------------------------------------
## rotor(tale)
- diameter : 0.14 m
- tip radius : 0.065 m
- hub radius : 0.005 m
- Thrust : 1Nf
- rpm : 6500
## aerofoil
- AG14
------------------------------------
# 最適化の定義
------------------------------------
## definition
- R : tip radius
- r : position @ rotor from its root
- chord : 翼弦長
- beta : 取付角
- r_R : r/R [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
- sn : section number
## variables
- chord[m] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
- beta[deg] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
## objective
- minimize efficiency
## optimization method
- GA
## individual
- size : sn × 2
- content
    - 0 to 5  : chord[m] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    - 6 to 11 : beta[deg] @ r/R[0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
## constraint
- chord : 0.005mm <= chords[i] <= 0.04
- beta  : 0deg <= betas[i] <= 22deg
## penalty
- Thrust >= 2Nf
