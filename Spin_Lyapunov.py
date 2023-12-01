import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from one_spin.SOT_sin import Sin_SOT

class Lyapunov():           #このメソッドはPhysRevB.100.224422を参考にしている。
    def __init__(self,alpha,gamma,B,S0,t,t_eval,Amp,omega,theta,Kx,Ky,Kz,beta,spin_start,spin_stop,lya_start_step,lya_cycle,delta_t1):
        self.alpha = alpha  # 緩和項をどれくらい大きく入れるか
        self.gamma = gamma  # LLg方程式のγ
        self.B = B  # 外部磁場(tで時間変化させてもいいよ)
        self.S0 = S0  # 初期のスピンの向き
        self.t = t  # どれくらいの時間でやるか
        self.t_eval = t_eval  # ステップ数はどうするか
        self.Amp = Amp
        self.omega = omega
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
        delta_theta = self.omega * self.delta_t1
        per_theta = np.array(self.theta) + delta_theta
        perturb_spin = {}

        spin = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega, self.theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        perturb_spin[0] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega,
                                  per_theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        spin.history()
        perturb_spin[0].history()
        spin0_log = spin.S.T
        per_spin0_log = perturb_spin[0].S.T

        S_t0 = spin0_log[self.lya_start]            #リヤプノフ指数を図るときの最初のt(最初から(t=0)距離を測り始めてしまうと最初の距離が0になってしまうため。)
        perS_t0 = per_spin0_log[self.lya_start]

        distance_t0 = self.distance(S_t0, perS_t0,self.delta_t1)
        a = 1/distance_t0
        print("D",distance_t0)
        S_t0_delt = spin0_log[self.lya_start + 1]
        perS_t0_delt = per_spin0_log[self.lya_start + 1]
        p1 = (self.distance(S_t0_delt, perS_t0_delt, self.delta_t1)) / distance_t0

        lya_expo = [p1 * a]
        dt = self.delta_t1
        ex_spin_log = spin0_log
        sum_log_ly = 0



        for k in range(1 , self.lya_cycle):
            #print(k)
            #print(lya_expo[k-1])
            n_omega = np.linalg.norm(self.omega)
            dt = dt / (lya_expo[k-1] * n_omega)
            #print(n_omega)
            d_theta = n_omega * np.linalg.norm(dt)
            ptheta = np.array(self.theta) + d_theta
            #print(ptheta)
            perturb_spin[k] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega ,ptheta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
            perturb_spin[k].history()
            spin_log = perturb_spin[k].S.T

            distance_kdash = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k],dt)
            #print(d_theta)
            distance = distance_kdash/lya_expo[k-1]

            pk = np.sqrt(distance_kdash**2 + d_theta**2)/(distance_t0)
            #print("pk",pk)
            lya_expo.append(pk)
            #print(lya_expo)
            print(np.log(pk))

            sum_log_ly += np.log(pk)


        step_width = (self.t[1] - self.t[0]) / len(self.t_eval)


        Lyapnov_exponent = sum_log_ly / (self.lya_cycle) + np.log(a)

        print("here",Lyapnov_exponent)








    def distance(self,S0,S1,del_t):
        S0 = np.array(S0)
        S1 = np.array(S1)
        print(S0,S1)
        #print(np.dot(S0,S1))

        Sdot = np.dot(S0,S1)
        Sdot = np.around(Sdot,12)
        l_cos_dis = np.arccos(Sdot)
        #print(l_cos_dis)
        #print( (np.linalg.norm(self.omega) * np.linalg.norm(self.delta_t1))**2)
        print(del_t)
        print(l_cos_dis, (np.linalg.norm(del_t)/2))
        dist = np.sqrt(l_cos_dis**2 + del_t[0]**2)
        return dist


if __name__ == '__main__':
    S0 = [0, 0, 1]

    t = [0, 5]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 12000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    mu_h_div_2e = [0.824, -21]
    sta_M = [1.4, 0]  # 飽和磁化(T)で入れる
    theta = [-2, -1]
    j = [2.5, 3]
    d = [1, -9]
    Hsn = mu_h_div_2e[0] * theta[0] * j[0] / (sta_M[0] * d[0])
    Hso = mu_h_div_2e[1] + theta[1] + j[1] - (sta_M[1] + d[1])
    Hs = Hsn * (10 ** Hso) * 1000 / gamma  # 最後の1000はmTにするため
    print(Hs)

    Lyap = Lyapunov(0.01, gamma, [0, 0, -4], S0, t, t_eval, [0, 0, 60],[0,0,6.8],[0,0,0], 0, 0, 12, 0, 0, 5, 100, 10, [0,0,0.01])
    Lyap.make_trajec()

