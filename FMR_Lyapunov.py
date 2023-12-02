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
        dt = (self.t[1] - self.t[0]) / len(self.t_eval)
        Amp_v = self.Amp * self.unit
        omega_v = self.omega * self.unit
        theta_v = self.theta * self.unit

        per_theta = (self.theta + self.delta_t1) * self.unit
        print(per_theta)
        perturb_spin = {}

        spin = FMR_spin(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v, theta_v, self.Kx,
                       self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        perturb_spin[0] = FMR_spin(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v,
                                  per_theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        plotB = [[0, 0, -1.2], [0, 0, 2.4]]

        #spin_gif = FMR_gif(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, Amp_v, omega_v, theta_v, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop, plotB)
        #spin_gif.make_gif()

        spin.history()
        perturb_spin[0].history()
        spin0_log = spin.S.T
        per_spin0_log = perturb_spin[0].S.T

        S_t0 = spin0_log[self.lya_start]            #リヤプノフ指数を図るときの最初のt(最初から(t=0)距離を測り始めてしまうと最初の距離が0になってしまうため。)
        perS_t0 = per_spin0_log[self.lya_start]

        print(spin0_log,per_spin0_log)

        #print(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega, self.theta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
        print(per_theta)

        distance_t0 = self.distance(S_t0, perS_t0, self.delta_t1)
        epsi = distance_t0
        print("log epsi", np.log(epsi))
        S_t0_delt = spin0_log[self.lya_start + 1]
        perS_t0_delt = per_spin0_log[self.lya_start + 1]
        p1 = self.distance(S_t0_delt, perS_t0_delt, self.delta_t1) /epsi

        lya_expo = [p1]
        d_theta = self.delta_t1
        ex_spin_log = spin0_log
        sum_log_ly = np.log(p1)
        mcrom = np.cross(S_t0_delt, perS_t0_delt)
        nmcro = mcrom / np.linalg.norm(mcrom)
        l_kdash = self.l(S_t0_delt, perS_t0_delt)
        l_k = l_kdash / p1




        for k in range(self.lya_cycle):
            d_theta = d_theta / (lya_expo[k - 1])
            # print(n_omega)
            k_start = dt * (self.lya_start + k)
            k_start_theta = self.omega * k_start
            ptheta = (self.theta + d_theta + k_start_theta) * self.unit





            cos = np.cos(l_k)
            sin = np.sin(l_k)
            R = np.array([[(nmcro[0] ** 2) * (1 - cos) + cos, (nmcro[0] * nmcro[1]) * (1 - cos) - nmcro[2] * sin,
                           (nmcro[0] * nmcro[2]) * (1 - cos) + nmcro[1] * sin],
                          [(nmcro[0] * nmcro[1]) * (1 - cos) + nmcro[2] * sin, (nmcro[1] ** 2) * (1 - cos) + cos,
                           (nmcro[1] * nmcro[2]) * (1 - cos) - nmcro[0] * sin],
                          [(nmcro[0] * nmcro[2]) * (1 - cos) - nmcro[1] * sin,
                           (nmcro[1] * nmcro[2]) * (1 - cos) + nmcro[0] * sin, (nmcro[2] ** 2) * (1 - cos) + cos]])
            m_kdelt = np.dot(R,spin0_log[self.lya_start + k-1])
            # print(ptheta)
            print("eになるはず",epsi,self.distance(spin0_log[self.lya_start + k - 1], m_kdelt, d_theta))

            perturb_spin[k] = FMR_spin(self.alpha, self.gamma, self.B, m_kdelt, self.t, self.t_eval, Amp_v, omega_v,
                                      ptheta, self.Kx, self.Ky, self.Kz, self.beta, self.spin_start, self.spin_stop)
            perturb_spin[k].history()
            spin_log = perturb_spin[k].S.T

            distance_kdash = self.distance(spin0_log[self.lya_start + k], spin_log[1], d_theta)
            # print(d_theta)
            distance = distance_kdash / lya_expo[k - 1]

            pk = self.distance(spin0_log[self.lya_start + k], spin_log[1], d_theta)/epsi
            # print("pk",pk)
            lya_expo.append(pk)
            print("lpk",np.log(pk))

            sum_log_ly += np.log(pk)

            mcrom = np.cross(spin0_log[self.lya_start + k], spin_log[1])
            nmcro = mcrom / np.linalg.norm(mcrom)
            l_kdash = self.l(spin0_log[self.lya_start + k], spin_log[1])
            l_k = l_kdash  / pk

            ex_m = spin_log[self.lya_start + k]


        step_width = (self.t[1] - self.t[0]) / len(self.t_eval)


        Lyapnov_exponent = sum_log_ly / (self.lya_cycle)
        print("here",Lyapnov_exponent)

if __name__ == '__main__':
    S0 = [4/5, 0, 3/5]

    t = [0, 1]  # t(時間)が0〜100まで動き、その時のfを求める。
    t_eval = np.linspace(*t, 10000)

    plotB = [[0, 0, -1.2], [0, 0, 2.4]]

    gamma = 0.17
    Amp = 8000
    omega = 3
    theta = [0, 0, 0]
    unit = [0, 0, 1]
    Kx = 0
    Ky = 0
    Kz = 4000
    delta_t1 = 0.01

    Lyap = FMR_lyapunov(0.005, gamma, [0, 0, 3300], S0, t, t_eval, Amp, omega, theta, unit,  Kx,Ky, Kz, 0, 0, 50, 300, 30, 0.01)
    Lyap.make_trajec()
