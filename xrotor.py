import numpy as np
import subprocess
import os

class Xrotor(object):
    """
    xrotorを操作するクラス
    各メソッドはxrotorにおけるコマンドと対応する
    # メソッド
    ## コマンド
    実行コマンドを__commandに登録する
    - impo(モデルファイル名)
        rotorモデル作成
    - oper(出力ファイル名)
        各種解析条件は属性のセッターで設定
        (velo,rpmなど)
    ## その他
    - call
        __commandに登録したコマンド実行
        このメソッドは最後に必ず呼び出すこと
        さもなければ何も起きない


    #属性
    - __command(str)
        コマンドを記憶する
    - velo
        解析する流入速度
    - rpm
        解析するrpm
    - nb
        ブレード枚数
    - fs
        設計進行速度
    - dens
        大気密度
    """
    #oper\nn\n100\n\n
    def __init__(self, nb, fs):
        self.__n = 30
        self.__command = "plop\ng\n\n"
        self.__velo = 1.0
        self.__rpm = 160
        self.__nb = nb
        self.__fs = fs
        self.__dens = 1.226

    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, value):
        self.__n = value
        self.__command += "oper\nn\n{n}\n\n".format(n = self.__n)

    @property
    def command(self):
        return self.__command
    @property
    def nb(self):
        return self.__nb

    @nb.setter
    def nb(self,value):
        self.__nb = value

    @property
    def fs(self):
        return self.__fs

    @fs.setter
    def fs(self,value):
        self.__fs = value

    @property
    def dens(self):
        return self.__dens

    @dens.setter
    def dens(self, value):
        self.__dens = value
        self.__command += "dens\n{dens}\n".format(dens = self.__dens)

    def call(self,timeout=7):
        #---xfoilの呼び出し---
        ps = subprocess.Popen(['xrotor.exe'],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        try:
            self.__command += "quit\n"
            res = ps.communicate(bytes(self.__command,"ascii"), timeout=timeout)
        #発散などによってxfoilが無限ループに陥った際の対応
        #タイムアウトによって実現している
        except subprocess.TimeoutExpired:
            res = None
            ps.kill()
        return res

    def aero(self,fname):
        """
        aeroファイル読み込み
        # argument
            - fname(str)
                aeroファイルパス
                ファイル内容は下記リンクのAERO欄を参照
                http://web.mit.edu/drela/Public/web/xrotor/xrotor_doc.txt

                例
                 Sect# =   1 r/R =     0.0000
                 ========================================================================
                 1) Zero-lift alpha (deg):   0.00       7) Minimum Cd           : 0.0070
                 2) d(Cl)/d(alpha)       :  6.280       8) Cl at minimum Cd     : 0.150
                 3) d(Cl)/d(alpha)@stall :  0.100       9) d(Cd)/d(Cl**2)       : 0.0040
                 4) Maximum Cl           :  2.00       10) Reference Re number  : 2000000.
                 5) Minimum Cl           : -1.50       11) Re scaling exponent  : -0.2000
                 6) Cl increment to stall:  0.200      12) Cm                   : -0.100
                                                       13) Mcrit                :  0.620
                 ========================================================================
        """
        pipe = "aero\nread\n{fname}\n\n"\
        .format(
            fname = fname,
        )
        self.__command += pipe
        return pipe

    def impo(self, fname):
        """
        xrotorのimpoコマンドを登録
        # argument
            - fname(str)
                プロペラの各種形状データを記憶するファイルのパス
                ファイル内容は下記リンク
                http://www.esotec.org/sw/dl/Espara_doc.txt
                のIMPOの欄を参照
                同ディレクトリ内のmodel.pyでも生成可能
            - nb(int)
                ブレード枚数
            - fs(float)
                flight speed 対気速度[m/s]
        # return
            - pipe(str)
                xrotorに入力するコマンド
        """

        pipe = "impo\n{fname}\n{nb}\n{fs}\n"\
        .format(
            fname = fname,
            nb = self.nb,
            fs = self.fs
        )
        self.__command += pipe
        return pipe

    @property
    def velo(self):
        return self.__velo

    @velo.setter
    def velo(self,value):
        self.__velo = value


    @property
    def rpm(self):
        return self.__rpm

    @rpm.setter
    def rpm(self, value):
        self.__rpm = value

    def oper(self):
        pipe = "oper\nvelo {velo}\nrpm {rpm}\naddc\nvseq\n\n\n\n\n"\
        .format(
            velo = self.velo,
            rpm = self.rpm
        )
        self.__command += pipe

    def cput(self,fname):
        """
        output analysys data as file
        """
        pipe = "oper\ncput\n{fname}\n\n"\
        .format(
            fname = fname
        )
        self.__command += pipe

if __name__ == "__main__":
    import numpy as np
    xr = Xrotor(2, 1.0)
    xr.velo = 0.1
    xr.rpm = 7000
    xr.aero("HQ_0_10.txt")
    xr.impo("rotor.txt")
    xr.oper()
    xr.velo = 0.1
    xr.rpm = 6000
    xr.oper()
    xr.cput("test3.txt")
    xr.call()
    res = np.loadtxt("test3.txt", skiprows=3)
    print(res)
