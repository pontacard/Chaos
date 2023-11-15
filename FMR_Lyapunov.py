import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
from Spin_Lyapunov import Lyapunov
from one_spin.FMR import FMR_spin

class FMR_lyapunov(Lyapunov):
    def __init__(self,alpha,gamma,B,S0,t,t_eval,Amp,omega,theta,Kx,Ky,Kz,beta,SOT_start,SOT_stop,lya_start_step,lya_cycle,delta_t1):
        super().__init__(self,alpha,gamma,B,S0,t,t_eval,Amp,omega,theta,Kx,Ky,Kz,beta,SOT_start,SOT_stop, lya_start_step, lya_cycle, delta_t1)


    def make_trajec(self):
        delta_theta = self.omega * self.delta_t1
        per_theta = np.array(self.theta) + delta_theta
        perturb_spin = {}

        spin = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega, self.theta, self.Kx, self.Ky, self.Kz, self.beta, self.SOT_start, self.SOT_stop)
        perturb_spin[0] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega,
                                  per_theta, self.Kx, self.Ky, self.Kz, self.beta, self.SOT_start, self.SOT_stop)
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
        p1 = (self.distance(S_t0_delt, perS_t0_delt, self.delta_t1))

        lya_expo = [p1]
        dt = self.delta_t1
        ex_spin_log = spin0_log
        sum_log_ly = 0



        for k in range(1 , self.lya_cycle):
            #print(k)
            #print(lya_expo[k-1])
            n_omega = np.linalg.norm(self.omega)
            dt = dt / (lya_expo[k-1] * n_omega * a)
            #print(n_omega)
            d_theta = n_omega * np.linalg.norm(dt)
            ptheta = np.array(self.theta) + d_theta
            perturb_spin[k] = Sin_SOT(self.alpha, self.gamma, self.B, self.S0, self.t, self.t_eval, self.Amp, self.omega ,ptheta, self.Kx, self.Ky, self.Kz, self.beta, self.SOT_start, self.SOT_stop)
            perturb_spin[k].history()
            spin_log = perturb_spin[k].S.T

            distance_kdash = self.distance(spin0_log[self.lya_start + k], spin_log[self.lya_start + k],dt)
            #print(d_theta)
            distance = distance_kdash/lya_expo[k-1]

            pk = np.sqrt(distance_kdash**2 + d_theta**2)
            #print("pk",pk)
            lya_expo.append(pk)
            #print(lya_expo)
            print(np.log(pk))

            sum_log_ly += np.log(pk)


        step_width = (self.t[1] - self.t[0]) / len(self.t_eval)


        Lyapnov_exponent = sum_log_ly / (self.lya_cycle) - np.log(a)

        print("here",Lyapnov_exponent)

