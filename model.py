import  numpy as np
class RotorModel(object):
    """
    ロータ形状モデルを定義するクラス。
    writefileメソッドによって生成されるファイルには、
    http://www.esotec.org/sw/dl/Espara_doc.txt
    のIMPOの欄に示される形式で、
    ロータ形状を定義する情報が書き出される。
    このファイルはxrotorのimpoコマンドで読み込むことができる。

    # attributes
        - output(str)
            ファイルに書き出す文字列
        - name(str)
            ロータの名前
        - tipr(float)
            プロペラ半径
        - hubr(float)
            ハブ半径
        - sn(int)
            半径方向のロータ分割数
            radii,chords,betas配列のサイズに等しい
        - radii(float list)
            回転中心からの半径方向の距離[mm]
        - chords(float list)
            radii[i]における翼弦長[mm]がchords[i]
        - betas(float list)
            radii[i]における翼弦長[mm]がbetas[i]
    # method
        - writefile
            ## argument
                - fname
                    ロータ形状を記憶するファイルパス
    
    """
    def __init__(self, rotorname, tipr, hubr, sn, radii, chords, betas):
        self.__output = ""
        self.__name = rotorname
        self.__tipr = tipr
        self.__hubr = hubr
        self.__sn = sn#セクションナンバー
        self.__betas = betas# array like
        self.__chords = chords# array like
        self.__radii = radii# array like

    @property
    def rotorname(self):
        return self.__name

    @rotorname.setter
    def rotorname(self, value):
        self.__name = value

    @property
    def tipr(self):
        return self.__tipr

    @tipr.setter
    def tipr(self, value):
        self.__tipr = value

    @property
    def hubr(self):
        return self.__hubr

    @hubr.setter
    def hubr(self, value):
        self.__hubr = value

    @property
    def secnumber(self):
        return self.__sn

    @secnumber.setter
    def secnumber(self, value):
        self.__sn = value

    @property
    def betas(self):
        return self.__betas

    @betas.setter
    def betas(self, value):
        if isinstance(value,list):
            value = np.ndarray(value)
        if isinstance(value, np.ndarray):
            self.__betas = value
            self.secnumber = len(value)
        else:
            raise Exception("cannot accept without array like object")


    @property
    def chords(self):
        return self.__chords

    @chords.setter
    def chords(self, value):
        if isinstance(value,list):
            value = np.ndarray(value)
        if isinstance(value, np.ndarray):
            self.__chords = value
            self.secnumber = len(value)
        else:
            raise Exception("cannot accept without array like object")

    @property
    def radii(self):
        return self.__radii

    @radii.setter
    def radii(self, value):
        if isinstance(value, list):
            value = np.ndarray(value)
        if isinstance(value, np.ndarray):
            self.__radii = value
            self.secnumber = len(value)
        else:
            raise Exception("cannot accept without array like object")

    @property
    def output(self):
        if len(self.radii) != len(self.betas) or len(self.betas) != len(self.chords) or len(self.chords) != len(self.radii):
            raise Exception("array sizes between radii, betas and chords must be the same")
        self.__output +="  " + str(self.rotorname) + "\n"
        self.__output += "  {tipr:<15.3f}    {hubr:<15.3f}    {secnumber:<15}\n".format(tipr = self.tipr, hubr = self.hubr, secnumber = self.secnumber)
        for r,c,b in zip(self.radii, self.chords, self.betas):
            self.__output += "  {r:<15.3f}    {c:<15.3f}    {b:<15.3f}\n".format(r = r, c = c, b = b)
        return self.__output

    def writefile(self, fname):
        with open(fname, mode='w') as f:
            f.write(self.output)

if __name__ == "__main__":
    r_R = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    init = [2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02, 2.0000e-02,
           2.0000e-02, 7.3204e+01, 5.2122e+01, 3.8006e+01, 2.9634e+01,
           2.3030e+01, 2.0777e+01]
    radii = [a*0.065 for a in r_R]
    chords = init[:int(len(init)/2)]
    betas = init[int(len(init)/2):]
    model = RotorModel("rotor",0.065,0.005,6,radii, chords, betas)
    model.writefile("rotor.txt")
