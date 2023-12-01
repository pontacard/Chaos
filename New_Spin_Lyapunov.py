import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from one_spin.SOT_sin import Sin_SOT

class Lyapunov():           #このメソッドはPhysRevB.100.224422を参考にしている。
    def __init__(self,alpha,gamma,B,S0,t,t_eval,Amp,omega,theta,unit,Kx,Ky,Kz,beta,spin_start,spin_stop,lya_start_step,lya_cycle,delta_t1):
        self.alpha = alpha  # 緩和項をどれくらい大きく入れるか
        self.gamma = gamma  # LLg方程式のγ
        self.B = B  # 外部磁場(tで時間変化させてもいいよ)
        self.S0 = S0  # 初期のスピンの向き
        self.t = t  # どれくらいの時間でやるか
        self.t_eval = t_eval  # ステップ数はどうするか
        self.Amp = Amp
        self.omega = omega
        self.unit = np.array(unit) #振動外場の単位ベクトル(振動外場がx方向を持っていたら[1,0,0]xyだったら[1,1,0]
        self.Kx = Kx
        self.Ky = Ky
        self.Kz = Kz
        self.beta = beta
        self.spin_start = spin_start
        self.spin_stop = spin_stop
        self.lya_start = lya_start_step
        self.lya_cycle = lya_cycle
        self.delta_t1 = np.array(delta_t1)
        self.theta = theta



    def make_trajec(self):
        Amp_v = self.Amp * self.unit
        omega_v = self.omega * self.unit
        theta_v = self.theta * self.unit

        delta_theta = self.omega * self.delta_t1
        per_theta = (self.theta + delta_theta) * self.unit
        print(per_theta)
        perturb_spin = {}

        spin = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v, theta_v, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        perturb_spin[0] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v,
                                  per_theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        spin.history()
        perturb_spin[0].history()
        spin0_log = spin.S.T
        per_spin0_log = perturb_spin[0].S.T

        S_t0 = spin0_log[self.lya_start]            #リヤプノフ指数を図るときの最初のt(最初から(t=0)距離を測り始めてしまうと最初の距離が0になってしまうため。)
        perS_t0 = per_spin0_log[self.lya_start]

        distance_t0 = self.distance(S_t0, perS_t0,delta_theta)
        epsi = distance_t0
        print("D",distance_t0)
        S_t0_delt = spin0_log[self.lya_start + 1]
        perS_t0_delt = per_spin0_log[self.lya_start + 1]
        p1 = self.distance(S_t0_delt, perS_t0_delt, delta_theta)

        lya_expo = [p1]
        d_theta = self.delta_t1
        ex_spin_log = spin0_log
        sum_log_ly = np.log(p1)



        for k in range(1 , self.lya_cycle):
            #print(k)

            d_theta = d_theta * epsi / (lya_expo[k-1]  * self.omega)
            #print(n_omega)

            ptheta = (self.theta + d_theta) * self.unit
            #print(ptheta)
            perturb_spin[k] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v ,ptheta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
            perturb_spin[k].history()
            spin_log = perturb_spin[k].S.T

            distance_kdash = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k],d_theta)
            #print(d_theta)
            distance = distance_kdash/lya_expo[k-1]

            pk = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k],d_theta)
            #print("pk",pk)
            lya_expo.append(pk)
            #print(lya_expo)
            print(np.log(pk))

            sum_log_ly += np.log(pk)


        step_width = (self.t[1] - self.t[0]) / len(self.t_eval)


        Lyapnov_exponent = sum_log_ly / (self.lya_cycle) - np.log(epsi)

        print("here",Lyapnov_exponent)








    def distance(self,S0,S1,del_t):
        S0 = np.array(S0)
        S1 = np.array(S1)
        print(S0,S1)
        #print(np.dot(S0,S1))

        Sdot = np.dot(S0,S1)
        Sdot = np.around(Sdot,12)
        l_cos_dis = np.arccos(Sdot)
        l_cross_dis = (1-Sdot)/2
        #print(l_cos_dis)
        #print( (np.linalg.norm(self.omega) * np.linalg.norm(self.delta_t1))**2)
        #print(del_t)
        print(l_cos_dis**2 , del_t**2)
        dist = np.sqrt(l_cross_dis + del_t**2)
        return dist


if __name__ == '__main__':
    S0 = [0.001, 0, 1]

    t = [0, 5]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 12000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    Amp = 6000
    omega = 3
    theta = [0,0,0]
    unit = [0,0,1]
    Kx = 0
    Ky = 0
    Kz = 5
    delta_t1 = 0.01

    Lyap = Lyapunov(0.01, gamma, [0, 0, -4], S0, t, t_eval, Amp, omega, theta, unit,  Kx,Ky, Kz, 0, 0, 5, 100, 10, delta_t1)
    Lyap.make_trajec()
