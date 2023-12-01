import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from New_Spin_Lyapunov import Lyapunov
from one_spin.FMR import FMR_spin
from one_spin.FMR_gif import FMR_gif

class FMR_lyapunov(Lyapunov):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,Amp,omega,theta ,unit ,Kx,Ky,Kz,beta,spin_start,spin_stop,lya_start_step,lya_cycle,delta_t1):
        super().__init__(alpha,gamma,B,S0,t,t_eval,Amp,omega,theta,unit ,Kx,Ky,Kz,beta,spin_start,spin_stop, lya_start_step, lya_cycle, delta_t1)


    def make_trajec(self):
        Amp_v = self.Amp * self.unit
        omega_v = self.omega * self.unit
        theta_v = self.theta * self.unit

        delta_theta = self.omega * self.delta_t1
        per_theta = (self.theta + delta_theta) * self.unit
        print(per_theta)
        perturb_spin = {}

        spin = FMR_spin(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v, theta_v, self.Kx,
                       self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        perturb_spin[0] = FMR_spin(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v,
                                  per_theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        plotB = [[0, 0, -1.2], [0, 0, 2.4]]

        spin_gif = FMR_gif(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v, theta_v, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop, plotB)
        spin_gif.make_gif()

        spin.history()
        perturb_spin[0].history()
        spin0_log = spin.S.T
        per_spin0_log = perturb_spin[0].S.T

        S_t0 = spin0_log[self.lya_start]            #リヤプノフ指数を図るときの最初のt(最初から(t=0)距離を測り始めてしまうと最初の距離が0になってしまうため。)
        perS_t0 = per_spin0_log[self.lya_start]

        print(spin0_log,per_spin0_log)

        #print(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega, self.theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        print(per_theta)

        distance_t0 = self.distance(S_t0, perS_t0, delta_theta)
        epsi = distance_t0
        print("log epsi", np.log(epsi))
        S_t0_delt = spin0_log[self.lya_start + 1]
        perS_t0_delt = per_spin0_log[self.lya_start + 1]
        p1 = self.distance(S_t0_delt, perS_t0_delt, delta_theta)

        lya_expo = [p1]
        d_theta = self.delta_t1
        ex_spin_log = spin0_log
        sum_log_ly = np.log(p1)



        for k in range(1 , self.lya_cycle):
            d_theta = d_theta * epsi / (lya_expo[k - 1] * self.omega)
            # print(n_omega)

            ptheta = (self.theta + d_theta) * self.unit
            # print(ptheta)
            perturb_spin[k] = FMR_spin(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v,
                                      ptheta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
            perturb_spin[k].history()
            spin_log = perturb_spin[k].S.T

            distance_kdash = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k], d_theta)
            # print(d_theta)
            distance = distance_kdash / lya_expo[k - 1]

            pk = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k], d_theta)
            # print("pk",pk)
            lya_expo.append(pk)
            # print(lya_expo)
            print(np.log(pk))

            sum_log_ly += np.log(pk)


        step_width = (self.t[1] - self.t[0]) / len(self.t_eval)


        Lyapnov_exponent = sum_log_ly / (self.lya_cycle) - np.log(epsi)

        print("here",Lyapnov_exponent)

if __name__ == '__main__':
    S0 = [4/5, 0, 3/5]

    t = [0, 50]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 1200000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    Amp = 12000
    omega = 0.3
    theta = [0, 0, 0.02]
    unit = [0, 0, 1]
    Kx = 0
    Ky = 0
    Kz = 4000
    delta_t1 = 0.01

    Lyap = FMR_lyapunov(0.005, gamma, [0, 0, 3300], S0, t, t_eval, Amp, omega, theta, unit,  Kx,Ky, Kz, 0, 0, 50, 300, 10, 0.1)
    Lyap.make_trajec()
