# 21/12/06, satoko, トランスのみをsimできるようにコード修正
# 真空管のみをsim      : Opt = 0, TubeFlag = 1
# トランスのみをsim    : Opt = 1, TubeFlag = 0
# 真空管+トランスをsim : Opt = 1, TubeFlag = 1

# 21/12/08, satoko, とりあえず山中さんが測定したBHのCSVを読み込んでプロットできるようにした
# %%
import sys

import warnings

import numpy as np
import soundfile as sf

import itertools

from lib_WDF import *
from lib_WDF_Triode import *
from lib_WDF_Hysteresis_JA3_model import *

import csv

import matplotlib.pyplot as plt


#WDF_C800G_version = "1.0.7.2021.12.08"
WDF_C800G_version = "1.1.0.2023.7.11 from 1.0.7.2021.12.08"
WDF_C800G_version = "1.2.0.2023.7.12 from 1.0.7.2021.12.08"
USE_DEBUG = False   # True
USE_NEW_PARAM = True

def get_WDF_C800G_version():
    return WDF_C800G_version

def Usage():
    print('py sim_WDF_Trans.py [Input wav/csv/-sin:freq@Fs] [Tube:0/1/2/3/4/Ig0] [Trans:0/1/11] [gain:G/-db:] [Ltp_shift] [Hys_Amp1] [Hys_amp2] [If_Rp]')
    exit()
    return

# get_TubeSelect(args[2])
# [Tube]    : 0/1/2/3/4/Ig0 {0: withont Tube}, {Ig0: Ig always equal 0}, {1/2/3/4: Tube type (use Type 2)}
#           1---Nakabayashi, 2---Nakabayashi Rev1, 3---Nakabayashi Rev2, 4---NEC(Curve Tracer)
def get_TubeSelect(text):
    if text == 'Ig0':
        f_Ig0 = True
        TubeFlag = 1
        Tube_Select = 2
    else:
        f_Ig0 = False
        Tube_Select = int(text)
        if Tube_Select == 0:
            TubeFlag = 0
            Tube_Select = 2
        else:
            TubeFlag = 1
    return TubeFlag, Tube_Select, f_Ig0

# get_TransSelect(args[3])
# [Opt] : 0/1/2/10 {0: without Trans}, {1: Hysterisis}, {2: Linear Trans w/o Hysterisis}, {10: Hysterisis by Yoneda-san paramaeter}
def get_TransSelect(text):
    USE_YONEDASAN_PARAM = False
    USE_HYSTERESIS = False
    Ltp_shift = 132.0      # 新パラメータ検討結果　微分回路のHPFを補正するためのインダクタ
    Ltp_shift = 11.0      # 暫定 2023.7.12
    Hys_Amp1 = 84
    Hys_Amp2 = 0.05         # 新定数(10Hz) comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
    #If_Rp = 1e3             # 微分回路の内部抵抗 1e3 →　1e5, ここの電流源の値がΦになる
    id = int(text)  # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    if id >= 10:
        USE_YONEDASAN_PARAM = True
        USE_HYSTERESIS = True
        Ltp_shift = 4000.0      # 米田さんER　微分回路のHPFを補正するためのインダクタ
        Hys_Amp1 = 84
        Hys_Amp2 = 0.05         # 新定数(10Hz) comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
        Opt = 1
    elif id == 1:
        USE_HYSTERESIS = True
        Opt = 1
    elif id == 2:
        Opt = 1
    elif id == 0:
        Opt = 0
    else:
        Opt = -1
    return Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2

# get_gain(args[4])
# Input gain 1.414213562=0[dBV] or {"-db:": input by dBV}
def get_gain(text):
    if text[0:4] == '-db:':
        sDB = text[4:]
        fDB = float(sDB)
        gain = 1.4142135623731 * (10 ** (fDB/20.0))
    else:
        sDB = ''
        gain = float(text)
        if gain > 0:
            fDB = 20.0 * np.log10(gain/1.4142135623731)
            sDB = '{:.1f}'.format(fDB)
    return gain,sDB

Fs = 96000  # sample rate (Hz)

Nt1 = 8.69 # トランスの巻き線比, 1次側/2次側

CH = 2 # 入力ステレオ対応
WAVoutbit = "F32"
#WAVoutbit = "24"
#Normalize = 1
Normalize = 0
WAVoutvol   = float(30.0)                # Wavout Fs Voltage(0-p) : if Normalize =0, Voltage(0-p) at 0dBfs 
RL          = float(6.2e3)               # 2022/9/27変更　Load Resistance(orm), 最終段の負荷抵抗

args = sys.argv

#                      args[1]     2                3              4             5           6          7          8
# py sim_WDF_Trans.py [Input wav] [Tube:0/1/2/3/4] [Trans,0/1/11] [gain:'-db:'] [Ltp_shift] [Hys_Amp1] [Hys_amp2] [If_Rp]
if len(sys.argv) == 4:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
elif len(sys.argv) == 5:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
    gain,sDB = get_gain(args[4]) # Input Voltage(0-p)
elif len(sys.argv) == 6:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
    gain,sDB = get_gain(args[4]) # Input Voltage(0-p)
    Ltp_shift   = float(args[5])   # 微分回路のHPFを補正するためのインダクタ
elif len(sys.argv) == 7:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
    gain,sDB = get_gain(args[4]) # Input Voltage(0-p)
    Ltp_shift   = float(args[5])   # 微分回路のHPFを補正するためのインダクタ
    Hys_Amp1    = float(args[6])   # 新定数(10Hz) comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
elif len(sys.argv) == 8:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
    gain,sDB = get_gain(args[4]) # Input Voltage(0-p)
    Ltp_shift   = float(args[5])   # 微分回路のHPFを補正するためのインダクタ
    Hys_Amp1    = float(args[6])   # 新定数(10Hz) comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
    Hys_Amp2    = float(args[7])   # comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
elif len(sys.argv) == 9:
    fname       = args[1]        # 入力wavファイル名
    TubeFlag, Tube_Select, f_Ig0    = get_TubeSelect(args[2])   # 1,2,3 : w/ tube, Ig0 : Ig=0, 0 : w/o Tube, only trans WDF 
    # Opt : 1 Hysteresis Trance, 2 Liner Trans, 0 --- trans less, 1x Use yonedasan parameter
    Opt, USE_HYSTERESIS, USE_YONEDASAN_PARAM, Ltp_shift, Hys_Amp1, Hys_Amp2 = get_TransSelect(args[3])
    gain,sDB = get_gain(args[4]) # Input Voltage(0-p)
    Ltp_shift   = float(args[5])   # 微分回路のHPFを補正するためのインダクタ
    Hys_Amp1    = float(args[6])   # 新定数(10Hz) comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
    Hys_Amp2    = float(args[7])   # comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
    If_Rp       = float(args[8])   # 微分回路の内部抵抗 1e3 →　1e5, ここの電流源の値がΦになる
else:
    Usage()

'''
print(len(args))
if len(sys.argv) == 10:
    fname       = args[1]         # 入力wavファイル名
    Tube_Select = int(args[2])    # 1---Nakabayashi, 2---Nakabayashi Rev1, 3---Nakabayashi Rev2, 4---NEC(Curve Tracer)
    Opt         = int(args[3])    # Opt : 1 --- Output trance 0 --- trans less
    gain        = float(args[4])  # Input Voltage(0-p)
    RL          = float(args[5])  # Load Resistance(orm), 最終段の負荷抵抗
    WAVoutbit   = str(args[6])    # Wavout bit : 16/24/32/F32
    Normalize   = int(args[7])    # Wavout Normalize : Normalize =1, No =0
    WAVoutvol   = float(args[8])  # Wavout Fs Voltage(0-p) : if Normalize =0, Voltage(0-p) at 0dBfs 
    TubeFlag    = int(args[9])    # 1 : w/ tube, 0 : w/o Tube, only trans WDF 
elif len(sys.argv) == 1: # 引数無しの場合は固定パラメータ, jupyterなどでセル毎に実行する場合は、print(len(args)) の数字を入力してください
    # 真空管のみをsim      : Opt = 0, TubeFlag = 1
    # トランスのみをsim    : Opt = 1, TubeFlag = 0
    # 真空管+トランスをsim : Opt = 1, TubeFlag = 1
    print("Warning: Without argument, fixed parameter is used in this script. usage : WDF_C-800G_6AU6A_opt_wav.py <Wav File> <6AU6> <Opt> <Input Voltage(0-p)> <Load Resistance(orm)> <Wavout bit> <Wavout Normalize> <Wavout Fs Voltage(0-p)> <TubeFlag>\n")
    fname       = "sin1kHz_96k_24bit_5s.wav" # 入力wavファイル名
    #fname       = "sin10Hz_96k_24bit_5s.wav" # 入力wavファイル名
    Tube_Select = int(2)                     # 1---Nakabayashi, 2---Nakabayashi Rev1, 3---Nakabayashi Rev2, 4---NEC(Curve Tracer)
    Opt         = int(1)                     # Opt : 1 --- Output trance 0 --- trans less
    gain        = float(0.071046877/2)       # -32dBV (0.071046877/2)[V 0-p] 16dBV=17.84/2 V(0 to peak) # 0dBV -> 1.414V(0-peak)
    RL          = float(12.4e3)              # Load Resistance(orm), 最終段の負荷抵抗
    RL          = float(6.2e3)               # 2022/9/27変更　Load Resistance(orm), 最終段の負荷抵抗
    WAVoutbit   = "F32"                      # Wavout bit : 16/24/32/F32
    Normalize   = int(1)                     # Wavout Normalize : Normalize =1, No =0
    WAVoutvol   = float(30.0)                # Wavout Fs Voltage(0-p) : if Normalize =0, Voltage(0-p) at 0dBfs 
    TubeFlag    = int(0)                     # 1 : w/ tube, 0 : w/o Tube, only trans WDF 
else:
    print("usage : WDF_C-800G_6AU6A_opt_wav.py <Wav File> <6AU6> <Opt> <Input Voltage(0-p)> <Load Resistance(orm)> <Wavout bit> <Wavout Normalize> <Wavout Fs Voltage(0-p)> <TubeFlag>\n")
    sys.exit()
'''

print("Aurgument list")
print("fname       : ", fname)
print("Tube_select : ", Tube_Select)
print("Opt         : ", Opt)
print("gain        : ", gain)
print("RL          : ", RL)
print("WAVooutbit  : ",WAVoutbit)
print("Normalize   : ",Normalize)
print("WAVoutvol   : ",WAVoutvol)
print("TubeFlag    : ",TubeFlag)
print("Ltp_shift =",Ltp_shift)
#print("If_Rp[ohm] =",If_Rp)
print("Hys_Amp1 =",Hys_Amp1)
print("Hys_Amp2 =",Hys_Amp2)
print("\n")

print("\nWDF_C800G Ver.     :",get_WDF_C800G_version())
print("lib_WDF Ver.       :",get_WDF_version())
print("lib_WDF_Triode Ver.:",get_WDF_Triode_version())
print("lib_WDF_Hystereys_JA3_model Ver.:",get_WDF_Hystereys_JA_model_version())
print("\n")

# 6AU6A Triode Connection


# Nakabayashi Tube_Select =1
# G = 0.0020681021, muc = 27.099343, alpha = 0.58944101, Ego = 0.24107953, Cgp = 2.2p, Cgk = 3.3p, Cpk = 5.0p

# Nakabayashi Rev1 Tube_Select =2
# G = 0.0020681021, muc = 27.099343/1.1, alpha = 0.58944101/1.1, Ego = 0.24107953, Cgp = 2.5p, Cgk =3.3p, Cpk = 5.0p

# Nakabayashi Rev2 Tube_Select =3
# G = 0.0020681021, muc = 27.099343/1.1, alpha = 0.58944101/1.1, Ego = 0.24107953, Cgp = 3.2p, Cgk =3.3p, Cpk = 5.0p

# NEC Curve Tracer Tube_Select =4
# G = 0.0024870748, muc = 28.7, alpha = 0.607, Ego = 0.52272826, Cgp = 2.2p, Cgk = 3.3p, Cpk = 5.0p




# 6AU6 Parameter Load
try:
    with open('C800G_6AU6_Param.txt', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=',')
        l = [row for row in reader]

except:
    print("6AU6 Parameter: defult")
    G_data = 0.0020681021
    muc_data = 27.099343
    alpha_data = 0.58944101
    Ego_data = 0.24107953
    Cgp_data = 2.2e-12
    Cgk_data = 3.3e-12
    Cpk_data = 5.0e-12

else:
    print("6AU6 Parameter:")
    G_data     = float(l[1][Tube_Select])
    muc_data   = float(l[2][Tube_Select])
    alpha_data = float(l[3][Tube_Select])
    Ego_data   = float(l[4][Tube_Select])
    Cgp_data   = float(l[5][Tube_Select])
    Cgk_data   = float(l[6][Tube_Select])
    Cpk_data   = float(l[7][Tube_Select])
 
print("G =", G_data)          # パービアンス, 電流の流れやすさを示す(http://ayumi.cava.jp/audio/pctube/node3.html)
print("muc =", muc_data)      # 増幅率
print("alpha =", alpha_data)  # 特性曲線近似パラメータ
print("Ego =", Ego_data)      # 特性曲線近似パラメータ, グリッドのオフセット電圧（ダイオードの0.6Vのようなもの）
print("Cgp =", Cgp_data)      # グリッドプレート間容量
print("Cgk =", Cgk_data)      # グリッドカソード間容量
print("Cpk =", Cpk_data)      # プレートカソード間容量
print("\n\n")

Tube = Triode_model(G_data, muc_data, alpha_data, Ego_data, Cgp_data, Cgk_data, Cpk_data) # 真空管モデル生成

# JA Model parameter
# https://www.comsol.jp/model/download/735431/models.acdc.vector_hysteresis_modeling.pdf
Hys_Ms = 1.31e6      # Saturation magnetization
Hys_a = 233          # Domain wall density
Hys_kp = 500 #374    # Pinning loss
Hys_alpha = 0.000562 # Inter-domain coupling
Hys_cr= 0.9 #0.736   # Magnetization reversibility

Hys_lambda=-30       # https://doc.comsol.com/5.5/doc/com.comsol.help.acdc/acdc_ug_theory.05.14.html, 元の論文探す, -30は米田さんが決めた値

Hys_Ns = 4000.0   # 1次側の巻線数, 本当は3700くらいかも？ヒステリシスは1次側で計算してる
Hys_Ld = 0.074    # 11/25 0.16から0.074へ変更 磁路長 H = NI/Ld
Hys_Sd = 92.0e-6  # 11/25 1.5e-4から92.0e-6へ変更　断面積 phi = BS


Hys_Hs = 500      # 飽和する上限, 磁界の強さの上限, 今回のコードでは使用してない

Hys_Eddy_COEF = 1e-12 # 渦電流係数
#Hys_Eddy_COEF = 1e-11 ～ 1e-14 # 渦電流係数 Takahiro Watanabe 2023.7.11
#Hys_Eddy_COEF = 10000

#Hys_Amp1 = 84
#Hys_Amp2 = 0.1
#Ltp_shift = 6000.0

if True:
    Hys_Amp1 = Hys_Amp1 * (0.074/0.16)   # comsol JAモデルをつかっているため、800Gに合うように入力側の辻褄合わせ
    Hys_Amp2 = Hys_Amp2 * (1.5 / 0.92) # 11/24 0.07から0.05へ変更,comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
else:
    Hys_Amp1 = 84 * (0.074/0.16)   # comsol JAモデルをつかっているため、800Gに合うように入力側の辻褄合わせ
    Hys_Amp2 = 0.05 * (1.5 / 0.92) # 11/24 0.07から0.05へ変更,comsol JAモデルをつかっているため、800Gに合うように出力側の辻褄合わせ
    Ltp_shift = 4000.0             # 11/24 7000.0から4000.0へ変更, 微分回路のHPFを補正するためのインダクタ
'''
if USE_NEW_PARAM:
    Ltp_shift = 132.0             # 11/24 7000.0から4000.0へ変更, 微分回路のHPFを補正するためのインダクタ
    Ltp_shift = 11.0             # 11/24 7000.0から4000.0へ変更, 微分回路のHPFを補正するためのインダクタ
'''

# Hys_Amp1 = Cal_Gain1, Hys_Amp2 = Cal_Gain2, Nsはここでは1にしてる
# ヒステリシスモデルの中で計算してるのではなく、渦電流モデルが登場してるため、
# ヒステリシス＋渦電流の等価回路の下の方で定義してるため
BH_model = Hysteresis_JA_model(Hys_Ms, Hys_a, Hys_kp, Hys_alpha, Hys_cr, Hys_lambda, 1, 1, Hys_Hs, 1, Hys_Ld, Hys_Sd)


# Transformer Parameter Load
if USE_NEW_PARAM:
    C800G_Opt_Param = 'C800G_Opt_Param_new.txt'
else:
    C800G_Opt_Param = 'C800G_Opt_Param_yonedasan.txt'
try:
    with open(C800G_Opt_Param, encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=',')
        l = [row for row in reader]

except:
    print("Transformer Parameter: defult")
    Lp_data = 5.68e5
    Ls_data = 7.52e3
    Ri_data = 3.92e6
    R1_data = 319.98
    Cp_data = 6.3e-11
    R2_data = 29.91
    Cs_data = 6.08e-10
    K_data  = 0.99999

else:
    print("Transformer Parameter:")
    Lp_data = float(l[0][1])
    Ls_data = float(l[1][1])
    Ri_data = float(l[2][1])
    R1_data = float(l[3][1])
    Cp_data = float(l[4][1])
    R2_data = float(l[5][1])
    Cs_data = float(l[6][1])
    K_data  = float(l[7][1])

# トランスの等価回路
print("Lp =", Lp_data) # 1次側インダクタンス
print("Ls =", Ls_data) # 2次側インダクタンス
print("Ri =", Ri_data) # 鉄損
print("R1 =", R1_data) # 1次巻線抵抗
print("Cp =", Cp_data) # 1次浮遊容量
print("R2 =", R2_data) # 2次巻線抵抗
print("Cs =", Cs_data) # 2次浮遊容量
print("K =" , K_data)  # 結合係数

print("\n")

data, Fs = sf.read(fname)

#START_OFFSET = Fs
START_OFFSET = 0


#print(data.shape)

N = data.shape[0]
if N == data.size:
    CH = 1
else:
    CH = data.shape[1]

#print("CH =", CH)
#print("N =", N)

if CH > 2:
    sys.exit()


if CH >1 :
    data= list(itertools.chain.from_iterable(data))
else :
    data  = data.tolist()


if CH == 2:
    data_L = data[0::2]  # L-ch
    data_R = data[1::2]  # R-ch

    data_L = np.r_[np.zeros(START_OFFSET), gain*np.array(data_L)]
    data_R = np.r_[np.zeros(START_OFFSET), gain*np.array(data_R)]

else:
    data_L = data[0::1]  # L-ch
    data_L = np.r_[np.zeros(START_OFFSET), gain*np.array(data_L)]

N = N+START_OFFSET 

output = np.zeros((N, CH))

Vpp = 230  # プレート電圧
Rp = 100e3 # プレート抵抗

# WDFの各素子
R0 = []
Ct = []
Ct2 = []
T1 = []
Rt2 = []
T2 = []
Lt2 = []
T3 = []
T4 = []
Ltp = []
Vl = []
T41 = []
T5 = []
Rti = []
T6 = []
Lt1 = []
T7 = []
Rt1 = []
T8 = []
Ctp = []
T9 = []

C0 = []
A1 = []
A2 = []
V = []
A3 = []
A4 = []
Rk = []
Ck = []
A5 = []


Vi = []
V2 = []
Ci = []
Ri = []
Ccap = []
Rg = []


A7 = []
A8 = []
A9 = []
A10 = []
A11 = []
A12 = []
A13 = []


Vmg = []
Lmg = []
MG1 = []


If = []
Lm = []
M1 = []


amg=np.zeros(CH)
bmg=np.zeros(CH)

# 米田さんER トランス編p.7, p.36を参照
for chnum in range(CH):
    R0.append(Resistor(RL))       # 最終段の負荷抵抗, R0 = RL

###### Output Transformer C-800G

    Lp = Lp_data
    Ls = Ls_data

# k coupling factor
    k = K_data

# Ct2 Secondary ray capacitance
    Ct2.append(Capacitor(Cs_data, Fs))         # トランス2次浮遊容量

    T1.append(Parallel(R0[chnum], Ct2[chnum])) # 最終段の負荷抵抗とトランス2次浮遊容量を並列接続

# Rt2 Secondary DC resistance
    Rt2.append(Resistor(R2_data))              # トランス2次巻線抵抗

    T2.append(Series(Rt2[chnum], T1[chnum]))   # トランス2次巻線抵抗とT1を直列接続

# Lt2 Secondary Leakage inductance
    Lt2.append(Inductor(Ls*(1-k), Fs))         # 2次漏れインダクタンス, (1-k)倍することで漏れを表現してる

    T3.append(Series(Lt2[chnum], T2[chnum]))   # 2次漏れインダクタンスとT2が直列接続

# Nt = 
    Nt = (Lp/Ls)**(0.5)                        # sqrt(1次インダクタンス/2次インダクタンス), 1次2次の比率

    T4.append(IdealTransformer(T3[chnum],Nt))  # T3とNtを理想トランスで接続

# Ltp Primary inductance
    Ltp.append(Inductor(k*Lp/Ltp_shift, Fs))   # ERだとLd = k*Lp/Ltp_shiftの部分, 線形だとk*Lpになる

    Vl.append(TerminatedVs(0, 1e-3))           # ERだとVfの部分, 内部抵抗1e-3の0V, 非線形インダクタを電圧源で表現してる
    T41.append(Series(Ltp[chnum], Vl[chnum]))  # ERだとLdとVfの直列接続部分,       非線形インダクタを電圧源で表現してる

    T5.append(Parallel(T41[chnum], T4[chnum])) # T41とT4を並列接続, 線形トランスの場合は T5.append(Parallel(Ltp[chnum], T4[chnum]))になる

# Rti Iron loss
    Rti.append(Resistor(Ri_data))              # 鉄損

    T6.append(Parallel(Rti[chnum], T5[chnum])) # 鉄損とT5を並列接続

# Lt1 Primary Leakage inductance
    Lt1.append(Inductor(Lp*(1-k), Fs))         # 1次漏れインダクタンス, (1-k)倍することで漏れを表現してる

    T7.append(Series(Lt1[chnum], T6[chnum]))   # 1次漏れインダクタンスとT6が直列接続


# Rt1 Primary DC resistance
    Rt1.append(Resistor(R1_data))              # トランス1次巻線抵抗

    T8.append(Series(Rt1[chnum], T7[chnum]))   # トランス1次巻線抵抗とT7が直列接続

#Ctp Primary stray capacitance
    Ctp.append(Capacitor(Cp_data, Fs))         # トランス1次浮遊容量

    T9.append(Parallel(Ctp[chnum], T8[chnum])) # トランス1次浮遊容量とT8が並列, T9がトランスの入力, トランスだけ検討したい場合はWDFの結合をここで止めればOK

###### Output Transformer end

    if TubeFlag == 1: # 真空管あり
        C0.append(Capacitor(1e-6, Fs))             # プレート側のコンデンサ, DCカット

        if Opt == 1 :
            A1.append(Series(C0[chnum],T9[chnum])) # トランスありの場合はプレート側のコンデンサとT9(トランス部分)が直列接続

        else:
            R0.append(Resistor(RL))
            A1.append(Series(C0[chnum],R0[chnum]))  # トランスなしの場合はコンデンサと出力負荷抵抗が直列接続

        A2.append(PolarityInverter(A1[chnum]))     # 極性反転

        V.append(TerminatedVs(Vpp,Rp))             # プレート電圧源, Vpp = 230V, Rp = 100k

        A3.append(Parallel(V[chnum], A2[chnum]))   # プレート電圧源とA2が並列接続


        Rk.append(Resistor(1e3))                   # カソード抵抗

        Ck.append(Capacitor(330e-6,Fs))            # カソード側のコンデンサ

        A5.append(Parallel(Ck[chnum],Rk[chnum]))   # カソード抵抗とカソード側のコンデンサが並列接続

 
        Vi.append(TerminatedVs(0,0.001))           # グリッド側回路の入力電圧、マイクへの入力信号
        Ci.append(Capacitor(47e-6,Fs))             # グリッド側回路の入力電圧のDCカット用コンデンサ
        Ri.append(Resistor(47e3))                  # グリッド側回路の抵抗
        Ccap.append(Capacitor(50e-12,Fs))          # マイクの容量
        Rg.append(Resistor(500e6))                 # グリッド入力抵抗

        V2.append(TerminatedVs(Vpp,180e3))         # マイクバイアス

        A7.append(Parallel(Ri[chnum], Ci[chnum]))  # グリッド側回路のRiとCiが並列接続, 米田さんER真空管編p.31参照

        A8.append(Parallel(A7[chnum], V2[chnum]))  # A7とマイクバイアスが並列接続

        A9.append(PolarityInverter(A8[chnum]))     # 極性反転

        A10.append(Series(Vi[chnum],A9[chnum]))    # A9とグリッド側入力電圧(マイクへの入力電圧)が直列接続
    
#    A11.append(PolarityInverter(A10[chnum]))

        A12.append(Series(Ccap[chnum],A10[chnum])) # A10とCcapが直列接続


        A13.append(Parallel(Rg[chnum], A12[chnum])) # A12とグリッド側回路の抵抗Rgと並列接続

        # LC Initial Setting
        # 真空管を無音状態で電源いれたときの状態の数値
        # 本当はトランスの中のインダクタの初期化もした方がいいかも？？？
        # 全部初期化しないと立ち上がり時のブオーンが直らないかも
        # コンデンサの両端にかかってる電圧, 実測値に合うようにpythonでパラメータを調整して、
        # 定常状態になったときのコンデンサの両端電圧
        C0[chnum].setWD(69.18780034567018)
        Ck[chnum].setWD(1.6089142338604283)
        Ccap[chnum].setWD(-47.69362504028323)
        Ci[chnum].setWD(47.572451753767815)
    
    else:
        Vi.append(TerminatedVs(0,0.001))           # トランスへの入力信号
        A1 = Parallel(Vi[chnum], T9[chnum])        # トランスと入力信号を並列接続

    # 渦電流に関する部分, Transformer_Magnetic circuit_model_20211019.pptx, p.3の右図参照
    Vmg.append(TerminatedVs(0,1e-12))           # 電圧源で内部抵付き, 12/03だと1e-12から1e-6に変更, 1e-12より小さいとNG
    Lmg.append(Inductor(Hys_Eddy_COEF, Fs))     # 渦電流の係数
    MG1.append(Series(Vmg[chnum],Lmg[chnum]))   # 渦電流を電圧と電流の直列接続で表現してる？
    amg[chnum] = 0 # ルートノードの入射波
    bmg[chnum] = 0 # ルートノードの反射波, p.3の右図の右上のヒステリシスのところがルートノード参照

    If.append(TerminatedIs(0,1e5)) # 11/25 微分回路の内部抵抗 1e3 →　1e5, ここの電流源の値がΦになる
    Lm.append(Inductor(1, Fs))
    M1.append(Parallel(If[chnum],Lm[chnum]))    #  微分ダミー回路 dΦ/dtを求める


# WDF_Jiles-Atherton_COMSOL_Model_Sin_Check_ind_circuit2_trans_60_60_eddy_211022.pyより追加
Hys_Vol_Vi  = np.zeros((N, CH)) # 11/24 佐藤追加 
Hys_Vol_Vl  = np.zeros((N, CH)) # 11/24 佐藤追加
Hys_Vol_R0  = np.zeros((N, CH)) # 11/24 佐藤追加
Hys_B       = np.zeros((N, CH)) # 11/24 佐藤追加
Hys_H       = np.zeros((N, CH)) # 11/24 佐藤追加
Hys_Vol_Phi = np.zeros((N, CH)) # 11/24 佐藤追加
for chnum in range(CH): 

    Vgk = -1.4877409505416912 # グリッドカソード間電圧,  436行目と同じでずれてるかも
    Vpk = 67.51847510219086   # プレートカソード間電圧 , 436行目と同じずれてるかもずれてるかも
    Vgo = 0                   # グリッドのオフセット値, 使われてない

    for num in range(N):
        if chnum == 0 :
            if TubeFlag == 1:
                Vi[chnum].E = Vgo+data_L[num] # 真空管回路への入力信号
            else:
                Vi[chnum].E = data_L[num]     # トランス回路への入力信号
        else:
            if TubeFlag == 1:
                Vi[chnum].E = Vgo+data_R[num]
            else:
                Vi[chnum].E = data_R[num]

        if TubeFlag == 1:
            ag = A13[chnum].WaveUp() # グリッド入射波, 末端からルートまで計算
            ap = A3[chnum].WaveUp()  # プレート入射波, 末端からルートまで計算
            ak = A5[chnum].WaveUp()  # カソード入射波, 末端からルートまで計算

            # Pentode calculations
            # 非線形処理, bg, bp, bk反射波
            [bg, bp, bk, Vgk, Vpk] = Tube.triodeNL(ap, A3[chnum].Rp, ag, A13[chnum].Rp, ak, A5[chnum].Rp, Vgk, Vpk, Fs)

            A13[chnum].WaveDown(bg) # グリッド側の回路, ルートから末端までを計算
            A3[chnum].WaveDown(bp)  # プレート側の回路, ルートから末端までを計算
            A5[chnum].WaveDown(bk)  # カソード側の回路, ルートから末端までを計算
        else:
            WU = A1.WaveUp() # 末端からルートまでを計算
            A1.WaveDown(WU)  # ルートから末端までを計算

        Vmg[chnum].E = -Vl[chnum].Current()*Hys_Ns*Hys_Amp1 # 渦電流回路の電圧源
        amg[chnum] = MG1[chnum].WaveUp()                    # ヒステリシスブロックの入射波
    
        BH_model.getHysteresis((amg[chnum]+bmg[chnum])/2, Fs) # 磁化Mは前の状態を覚えている, 磁化Mを計算

        def fb(bmg1):
            # ここでは磁化MをHとBを動かして収束させてる
            # ここでMを動かすとおかしいことになる
            return  BH_model.getHysteresis2((amg[chnum]+bmg1)/2, Fs) - (amg[chnum]-bmg1)/(2*MG1[chnum].Rp)

        if USE_DEBUG:
            prev_bmg = bmg[chnum]
            print('num:',num,'amg:',amg[chnum],'bmg:',prev_bmg,'MG1.Rp:',MG1[chnum].Rp,'fb(0):',fb(0))
            if num < 0:
                x = np.linspace(-2000000,200000,1000)
                y1 = BH_model.getHysteresis2((amg[chnum]+x)/2, Fs)
                y2 = (amg[chnum]-x)/(2*MG1[chnum].Rp)
                plt.plot(x,y1)
                plt.plot(x,y2)
                plt.plot(x,y1-y2)
                plt.show()
        bmg[chnum] = newton(fb, bmg[chnum], tol=1e-4, maxiter=100) # ヒステリシスブロックの反射波
        if USE_DEBUG:
            print('                     bmg =>',bmg[chnum],'fb:',fb(bmg[chnum]))
            gosa = fb(bmg[chnum])
            #if gosa > 1e-3 or gosa < -1e-3:
            print('fb(bmg[chnum]):',gosa)

        MG1[chnum].WaveDown(bmg[chnum]) # 渦電流回路のところ

        Hys_Phi = -Lmg[chnum].Current()*Hys_Amp2            # 磁束Φ = BS, 渦電流回路でだと電流がΦ, Amp2で出力の調整をしている
        Hys_B[num][chnum] = Hys_Phi / (Hys_Sd*Hys_Amp2)     # 11/24 B-Hカーブをプロットできるように追加, ΦからBを計算している
        Hys_H[num][chnum] = Vmg[chnum].Voltage() / Hys_Ld   # 11/24 B-Hカーブをプロットできるように追加
        If[chnum].Is = -Hys_Phi # モニターするときの辻褄合わせ


        WUM = M1[chnum].WaveUp() # 微分回路部分
        M1[chnum].WaveDown(WUM)  # 微分回路部分

        Hys_Vol_Phi[num][chnum] = Lm[chnum].Voltage() # 微分回路を通った値のためΦを微分してる→dΦ/dtでVになってる
        Vl[chnum].E = Hys_Vol_Phi[num][chnum]*Hys_Ns  # dΦ/dt*巻き数 = V電圧

        Hys_Vol_R0[num][chnum] = R0[chnum].Voltage() # 11/24 佐藤追加 最終段の負荷抵抗
        Hys_Vol_Vi[num][chnum] = T9[chnum].Voltage() # 11/24 佐藤追加 トランスの入力電圧
        Hys_Vol_Vl[num][chnum] = Vl[chnum].Voltage() # 11/24 佐藤追加 非線形ヒステリシスの電圧源の電圧

        output[num][chnum] = R0[chnum].Voltage()

if Opt != 1 :
    output = output/Nt1 # 巻き線比で割ってる


fsname = "WDF_C800G_6AU6_"

fsname  = fsname + '{}'.format(Tube_Select) +"_"


if Opt == 1 :
    fsname  = fsname + "Opt_Hysteresis_JA_"
else:
    fsname = fsname +"Optless_"


fsname  = fsname +'{:.2f}'.format(gain)+"V_"+'{:.2e}'.format(RL)+ "orm_"


fsname  = fsname + "Out_" +  WAVoutbit + "bit_"

if Normalize == 1:
    fsname = fsname + "Normalize_On_"
else :
    fsname = fsname + "Fs_"+ '{:.0f}'.format(WAVoutvol) +"V_"

if TubeFlag == 1:
    fsname = fsname + "w_Tube_"
else:
    fsname = fsname + "wo_Tube_"

fsname = fsname + fname


output = output[START_OFFSET:,:]

output_max = np.amax(np.abs(output))

if Normalize == 1:
    output = output/output_max
else:
    output = output/WAVoutvol
    output = np.clip(output, -1.0, 1.0)

_format = "WAV"

if WAVoutbit=="F32" :
    subtype = 'FLOAT'

elif WAVoutbit=="32" :
    subtype = 'PCM_32'

elif WAVoutbit=="24" :
    subtype = 'PCM_24'

else:
    subtype = 'PCM_16'

sf.write(fsname, output, Fs, format=_format, subtype=subtype)

exit()
## 以下、追加
# WDF_Jiles-Atherton_COMSOL_Model_Sin_Check_ind_circuit2_trans_60_60_eddy_211022.pyより追加
CSV_OUT = np.c_[Hys_H[0:len(data),:], Hys_B[0:len(data),:], data, Hys_Vol_R0[0:len(data),:]]
np.savetxt('JA_model_Out.csv', CSV_OUT, delimiter=',')

print("Processing finished. Ready to plot")

'''
# 全サンプルのプロット
plt.subplot(6,1,1)
#plt.plot(Hys_Vol_Vi) # トランスの入力電圧 プレート出力からDCカットコンデンサの後の電圧
plt.plot(data)       # 引数で与えたwavファイル
plt.title('Input') 
plt.xlim(0, len(data))

plt.subplot(6,1,2)
plt.plot(Hys_H)
plt.title('H')
plt.xlim(0, len(data))

plt.subplot(6,1,3)
plt.plot(Hys_B)
plt.title('B')
plt.xlim(0, len(data))

plt.subplot(6,1,4)
plt.plot(Hys_H,Hys_B)
plt.title('B-H')

plt.subplot(6,1,5)
plt.plot(Hys_Vol_Phi) # phi = BS
plt.title('Vol_Phi')
plt.xlim(0, len(data))

plt.subplot(6,1,6)
plt.plot(Hys_Vol_R0) # 負荷抵抗の両端電圧（最終出力）
plt.title('R0')
plt.xlim(0, len(data))

plt.show()
'''

# %%
# 1秒～1kHz 2周期分
Hz_forplot=1000
fig2 = plt.figure
plt.plot(Hys_H[(Fs*4):(Fs*4+int(Fs*2/Hz_forplot-1)),0],Hys_B[(Fs*1):(Fs*1+int(Fs*2/Hz_forplot-1)),0])
plt.title('B-H zoom')
plt.show()


fig3 = plt.figure
plt.subplot(6,1,1)
#plt.plot(Hys_Vol_Vi) # トランスの入力電圧 プレート出力からDCカットコンデンサの後の電圧
plt.plot(data)       # 引数で与えたwavファイル
plt.title('Input') 
plt.xlim((Fs*4), (Fs*4+int(Fs*2/Hz_forplot-1)))

plt.subplot(6,1,2)
plt.plot(Hys_H)
plt.title('H')
plt.xlim((Fs*4), (Fs*4+int(Fs*2/Hz_forplot-1)))

plt.subplot(6,1,3)
plt.plot(Hys_B)
plt.title('B')
plt.xlim((Fs*4), (Fs*4+int(Fs*2/Hz_forplot-1)))

plt.subplot(6,1,4)
plt.plot(Hys_H[(Fs*1):(Fs*1+int(Fs*2/Hz_forplot-1)),0],Hys_B[(Fs*1):(Fs*1+int(Fs*2/Hz_forplot-1)),0])
plt.title('B-H')

plt.subplot(6,1,5)
plt.plot(Hys_Vol_Phi) # phi = BS
plt.title('Vol_Phi')
plt.xlim((Fs*4), (Fs*4+int(Fs*2/Hz_forplot-1)))

plt.subplot(6,1,6)
plt.plot(Hys_Vol_R0) # 負荷抵抗の両端電圧（最終出力）
plt.title('R0')
plt.xlim((Fs*4), (Fs*4+int(Fs*2/Hz_forplot-1)))

plt.show()
# %%
stopper=0

# %%
# CSVのread
csv_file1 = open("10Hz_14r33dBV_Ch1.csv", "r", encoding="utf-8", errors="", newline="" )
meas_ch1_obj = csv.reader(csv_file1, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

csv_file2 = open("10Hz_14r33dBV_Ch2.csv", "r", encoding="utf-8", errors="", newline="" )
meas_ch2_obj = csv.reader(csv_file2, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

meas_H = np.zeros(10000) # 山中さんのCSVのデータ数10000のため
ii = 0
for row in meas_ch1_obj: # rowはリスト
    meas_H[ii] = row[4]  # CSVのE列のデータを取得
    ii += 1

meas_B = np.zeros(10000) # 山中さんのCSVのデータ数10000のため
ii = 0
for row in meas_ch2_obj: # rowはリスト
    meas_B[ii] = row[4]  # CSVのE列のデータを取得
    ii += 1
   
# BHカーブのプロット
plt.plot(meas_H, meas_B)
plt.title('B-H from measurement')
plt.show()


# %%
# WDFのBHと測定したBHをプロット
plt.plot(Hys_H[(Fs*4):(Fs*4+int(Fs*2/Hz_forplot-1)),0],Hys_B[(Fs*1):(Fs*1+int(Fs*2/Hz_forplot-1)),0], label = "WDF")
plt.plot(meas_H, meas_B, label = "Measurement")
plt.xlabel("H")
plt.ylabel("B")
plt.title('B-H, WDF VS Measurement')
plt.legend()
plt.show()