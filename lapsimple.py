import numpy as np
from racetrack import racetrack
from typing import TypedDict

class car(TypedDict):
    mu_x: float
    mu_y: float
    mass: float # kg
    power: float # PS
    cLA: float
    cDA: float
    rho: float


class lapsimple:
    maxvel = 1000

    def __init__(self, track: racetrack, car: car):
        self.track = track
        self.mu_x = car["mu_x"]
        self.mu_y = car["mu_y"]
        self.mass = car["mass"]
        self.power = car["power"]
        self.cLA = car["cLA"]
        self.cDA = car["cDA"]
        self.rho = car["rho"]


    def interpolate_ax(self, ay, v, dir):
        if dir in ["fwd", "f"]:
            ax_mx = self.mu_x*9.81 + (self.mu_x*self.cLA - self.cDA)*0.5*self.rho*v**2/self.mass
        elif dir in ["bwd", "b"]:
            ax_mx = self.mu_x*9.81 + (self.mu_x*self.cLA + self.cDA)*0.5*self.rho*v**2/self.mass

        ay_mx = self.get_ay_mx(v)
        # todo: catch inf by inf division error
        t = (1-ay**2/ay_mx**2)*ax_mx**2
        if t < 0:
            return 0
        ax = np.sqrt(t)
        #limit engine
        ax_mx_engine = self.power*1360/(v*self.mass)
        if np.isnan(ax):
            ax = 0
        if ax>ax_mx_engine and dir not in ["bwd", "b"]:
            ax = ax_mx_engine
        return ax
    
    def get_ay_mx(self, v):
        ay = self.mu_y*9.81 + 0.5*self.rho*v**2*self.cLA/self.mass
        return ay

    def get_maxV(self, k):
        k = abs(k)
        t = self.mu_y*9.81/(k-0.5*self.rho*self.cLA/self.mass)
        t[np.where(t<=0)] = np.inf
        v = np.sqrt(t)
        return v


    def intTrack(self, dir="fwd"):
        #INTTRACK integration along the track
        #   Detailed explanation goes here

        maxV = self.get_maxV(self.track.k)

        # initialize velocity vector as infinty
        v       = np.ones_like(self.track.k)*np.inf
        a_y     = np.ones_like(self.track.k)*np.inf
        a_x     = np.zeros_like(self.track.k)
        corners = self.track.corners
        # get maximum corner speed from vehicle model
        v[corners] = maxV[corners]
        is_inf = np.where(np.isinf(v))
        corners = np.setdiff1d(corners, is_inf)
        grad_u = self.track.u

        if dir in ["fwd", "f"]:
            print("forward integration")
            # forward integration
            s = corners[0]
            # start from first 'apex' and integrate forward
            
            for kk in range(s, len(v)-1): # start from last 'apex' and integrate forward
                a_y[kk] = v[kk]**2 * self.track.k[kk]
                a_x[kk] = self.interpolate_ax(a_y[kk],v[kk], dir=dir) # LOOKUP SHOULD NOT RETURN NEGATIVE VALUES!

                v_tmp = v[kk] + grad_u[kk]/v[kk]*a_x[kk]; # could be wrong
                            
                if v_tmp <= min(v[kk+1],maxV[kk+1]):
                    v[kk+1] = v_tmp
                else:
                    v[kk+1] = maxV[kk+1]
            # if circuit is continuous, the first velocity is also the last  
            v[0] = v[-1];
            for kk in range(0, s):
                a_y[kk] = v[kk]**2 * self.track.k[kk]
                a_x[kk] = self.interpolate_ax(a_y[kk],v[kk], dir=dir) # LOOKUP SHOULD NOT RETURN NEGATIVE VALUES!

                v_tmp = v[kk] + grad_u[kk]/v[kk]*a_x[kk]; #could be wrong
                            
                if v_tmp <= min(v[kk+1],maxV[kk+1]):
                    v[kk+1] = v_tmp
                else:
                    v[kk+1] = maxV[kk+1]
        
        #a_x(numel(v)) = 0; # not sure why this is done
        #  lot of code cuplication here, could probably be eliminated
            return v
        
        # backward integration
        elif dir in ["bwd", "b"]:
            print("backward integration")
            s = self.track.corners[-1]     
            # start from last 'apex' and integrate backwards
            for kk in range(s,1,-1): 
                a_y[kk] = v[kk]**2 * self.track.k[kk]
                a_x[kk] = self.interpolate_ax(a_y[kk], v[kk], dir=dir) # SHOULD NOT RETURN NEGATIVE VALUES! Not sure about negative ax_beta here

                v_tmp = v[kk] + grad_u[kk]/v[kk]*a_x[kk] #could be wrong
                
                if v_tmp <= min(v[kk-1],maxV[kk-1]):
                    v[kk-1] = v_tmp
                else:
                    v[kk-1] = maxV[kk-1]
            
            
            v[-1] = v[0]
            # start from last 'apex' and integrate backwards
            for kk in range(len(v)-1, s+1, -1): 
                a_y[kk] = v[kk]**2 * self.track.k[kk]
                a_x[kk] = self.interpolate_ax(a_y[kk], v[kk], dir=dir) # SHOULD NOT RETURN NEGATIVE VALUES! Not sure about negative ax_beta here
                v_tmp = v[kk] + grad_u[kk]/v[kk]*a_x[kk] #could be wrong
                            
                if v_tmp <= min(v[kk-1],maxV[kk-1]):
                    v[kk-1] = v_tmp
                else:
                    v[kk-1] = maxV[kk-1]
            return v
        
        else:
            print("Dir must be either fwd or bwd, {dir} not accepted")
            return 0