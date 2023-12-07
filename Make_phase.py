import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from one_spin.one_spin_LLG import A_spin
from one_spin.STO_gif import STO_gif
from one_spin.STO_2Dgraph import one_spin_2D
from chaos.STO_Lyapunov import STO_lyapunov
import pandas as pd

S0 = [1, 0, 0]

t = [0, 500]  # t(時間)が0〜100まで動き、その時のfを求める。
t_eval = np.linspace(*t, 70000)

plotB = [[0, 0, -1.2], [0, 0, 2.4]]

gamma = 0.17
mu_0 = 1.26
mu_h_div_2e = [0.824, -21]
sta_M = [1.2, 0]  # 飽和磁化(T)で入れる
d = [2, -9]


jdc = [2.2 , 10]
Hdcn = mu_h_div_2e[0] * jdc[0] / (sta_M[0] * d[0])
Hdco = mu_h_div_2e[1] + jdc[1] - (sta_M[1] + d[1])
Hdc = Hdcn * (10 ** Hdco) * 1000 / gamma  # 最後の1000はmTにするため

Kx = 0
Ky = 0
Kz = 1481 - 1448

f = open("Lyapunov.txt","w")


for j in [0.042,17.70,26.55,30.97,33.58,33.81,35.4,44.25]:
    for omega in [1.5,2,2.5,3,3.5,4]:

        print(j)
        jac = [j, 11]
        Hacn = mu_h_div_2e[0] * jac[0] / (sta_M[0] * d[0])
        Haco = mu_h_div_2e[1] + jac[1] - (sta_M[1] + d[1])
        Hac = Hacn * (10 ** Haco) * 1000 / gamma  # 最後の1000はmTにするため

        t = [0, 50]  # t(時間)が0〜100まで動き、その時のfを求める。
        t_eval = np.linspace(*t, 1000000)
        spin_graph = one_spin_2D(0.005, gamma, [0, 0, mu_0 * 159], S0, t, t_eval, [0, 0, 0], [omega, omega, 0], [0, 0, 0], Kx, Ky,
                    mu_0 * Kz, 0, 0, 50, [0, Hac, 0], [0, Hdc, 0], 0.288, 0.537)
        spin_graph.get_graph()

        t = [0, 10]  # t(時間)が0〜100まで動き、その時のfを求める。
        t_eval = np.linspace(*t, 10000)
        spin_gif = STO_gif(0.005, gamma, [0, 0, mu_0 * 159], S0, t, t_eval, [0, 0, 0], [omega, omega, 0], [0, 0, 0], Kx, Ky,
                    mu_0 * Kz, 0, 0, 10, [0, Hac, 0], [0, Hdc, 0], 0.288, 0.537)
        #spin_gif.make_gif()

        t = [0, 20]  # t(時間)が0〜100まで動き、その時のfを求める。
        t_eval = np.linspace(*t, 200000)
        Lya_expo = 0
        for lya_st in range(1,181):
            Lyap = STO_lyapunov(0.005, gamma, [0, 0, mu_0 * 159], S0, t, t_eval, 0,omega,[0,0,0],[0, 1, 0], Kx,Ky, mu_0 * Kz, 0, 0, 20,[0,Hac,0],[0,Hdc,0],0.288,0.537,10000 * lya_st,20,0.1)
            tem_lya = Lyap.make_trajec()
            Lya_expo += tem_lya

        Lya_expo = Lya_expo/180
        f = open('Lyapunov.txt' , "a")
        f.write(f'   {j} Am^-2  {omega} GHz  {Lya_expo}   ')
        f.close()




