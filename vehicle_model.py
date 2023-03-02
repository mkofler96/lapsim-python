import numpy as np
import configparser
import math
from TMSimple import TMSimple
from scipy import optimize

def dict2float(dict):
    new_dict = {}
    for k, v in dict.items():
        new_dict[k] = float(v)
    return new_dict

class vehicle_model:
    def __init__(self, vehicle_model_file):
        config = configparser.ConfigParser()
        config.read(vehicle_model_file)
        self.env = dict2float(config["environment"])
        self.general = dict2float(config["general"])
        self.aero = dict2float(config["aero"])
        self.suspension = dict2float(config["suspension"])
        self.powertrain = dict2float(config["powertrain"])
        self.tire = TMSimple(config["meta"]["tire_parameter_file"])

    def get_vehicle_forces(self, ax, ay, psi_dot, state):
        ## Einspurmodelle 10      S2016/DAQ/Literatur/10_Einspurmodelle
        lf = self.general["wheelbase"] * (1-self.general["weight_distr"]) 
        lr = self.general["wheelbase"] * self.general["weight_distr"] 

        h_w = self.general["rc_z"] 
        h_s = self.general["cog_z"] 

        # Rollen bercksichtigt, Nicken nicht
        # bentigte variablen
        S_fv =  self.general["track_f"]/2      # Normalabstand Feder Vorne
        S_fh =  self.general["track_r"]/2      # Normalabstand Feder Hinten
        S_stv = self.general["track_f"]/2      # Normalabstand Stabi Vorne
        S_sth = self.general["track_r"]/2      # Normalabstand Stabi Hinten
        C_fv =  self.suspension["k_spring_f"] /self.suspension["mr_spring_f"]**2       # Federsteifigkeit Vorne
        C_fv3 = self.suspension["k_spring_f3"]  /self.suspension["mr_spring_f3"] **2     # dritte Federsteifigkeit Vorne
        C_fh =  self.suspension["k_spring_r"]  /self.suspension["mr_spring_r"] **2       # Federsteifigkeit Hinten
        C_stv = self.suspension["k_arb_f"] /self.suspension["mr_arb_f"]**2       # Torsionssteifigkeit Vorne
        C_sth = self.suspension["k_arb_r"] /self.suspension["mr_arb_r"]**2       # Torsionssteifigkeit Hinten

        # Glg 10.104        nenner, zhler  um phi zu berechnen
        nn = self.general["mass"]*ay*(h_s-h_w) 
        zl = 2*(S_fv**2*C_fv+(C_stv*S_stv)/2+S_fh**2*C_fh+(C_sth*S_sth)/2) 
        phi = nn /zl        # Phi Berechnung

        # Glg 10.98 bis 7.101
        F_flv = -C_fv*S_fv*phi 
        F_frv = C_fv*S_fv*phi 
        F_flh = -C_fh*S_fh*phi 
        F_frh = C_fh*S_fh*phi 

        # Glg 10.102 und 10.103
        F_stlv = -0.5*C_stv*phi 
        F_strv = 0.5*C_stv*phi 
        F_stlh = -0.5*C_sth*phi 
        F_strh = 0.5*C_sth*phi 

        # Glg 10.105 - 10.108
        A_lv = F_flv + F_stlv 
        A_rv = F_frv + F_strv 
        A_lh = F_flh + F_stlh 
        A_rh = F_frh + F_strh 

        A = [A_lv, A_rv, A_lh, A_rh] 

        # Glg 10.109 und 10.110

        Mav = S_fv*(F_flv - F_frv) + S_stv*(F_stlv - F_strv) 
        Mah = S_fh*(F_flh - F_frh) + S_sth*(F_stlh - F_strh) 

        Fzmav = np.array([-Mav/self.general["track_f"], Mav/self.general["track_f"], -Mah/self.general["track_r"], Mah/self.general["track_r"]])

        ## Aerodynamics
        # cL Calculation depending on Ride Height

        cLA = self.aero["cl"]*self.aero["farea"]
        cDA = self.aero["cd"]*self.aero["farea"]
        AeroDistr = self.aero["distr"]
        # resulting Fz Aero
        AirDen = self.env["rho"]  # in kg/m**3
        FAero = (0.5  * cLA  * AirDen  * state.vx**2)*np.array([AeroDistr/2,AeroDistr/2,(1-AeroDistr)/2,(1-AeroDistr)/2])

        ## combining Aero and 10_Einspurmodelle
        # lr/(Wheelbase*2) wegen 1/4 masse auf rad und 2*hs/Wb wegen Hebel
        Mass = self.general["mass"]
        Wheelbase = self.general["wheelbase"]
        TrackF = self.general["track_f"]
        TrackR = self.general["track_r"]
        Fz_geom = np.array([Mass*(lr*9.81/(Wheelbase*2)-h_w*lr*ay/(TrackF*Wheelbase)),
            Mass*(lr*9.81/(Wheelbase*2)+h_w*lr*ay/(TrackF*Wheelbase)),
            Mass*(lf*9.81/(Wheelbase*2)-h_w*lf*ay/(TrackR*Wheelbase)),
            Mass*(lf*9.81/(Wheelbase*2)+h_w*lf*ay/(TrackR*Wheelbase))])

        # longitudinal weight transfer (w/o anti geometry, pitch center etc.)
        Fz_ax = Mass*0.5*np.array([-h_s*ax/(Wheelbase),-h_s*ax/(Wheelbase),h_s*ax/(Wheelbase),h_s*ax/(Wheelbase)])

        Fz = Fz_geom + Fz_ax + A + FAero  # eventuell + MAV (kontrollieren)

        # Long/Lat Accelerations
        FxAero = sum(FAero)
        
        # Rolling Force
        
        FxRoll = 0.015*sum(Fz,2)
        Fx = Mass*ax-FxAero - FxRoll
        Fy = Mass*ay

        vehicle_forces = {"Fx": Fx, "Fy": Fy, "Fz": Fz}
        return vehicle_forces
    

    def get_slip_angle(self, psi_dot, state):
        #claude skript Sept.2017 page 685
        a = self.general["wheelbase"] * (1-self.general["weight_distr"])
        b = self.general["wheelbase"] * self.general["weight_distr"]
        tf = self.general["track_f"] / 2
        tr = self.general["track_r"] / 2

        # 10 Einspurmodelle
        n = state.vx* np.tan(np.radians(state.beta)) + psi_dot * np.array([a,a,-b,-b])
        z = state.vx + psi_dot * np.array([-tf,tf,-tr,tr])
        betaT = np.arctan2(n,z)
        # static toe, negativ toe = toe in
        static_toe = np.array([self.suspension["toe_fl"], -self.suspension["toe_fr"], self.suspension["toe_rl"], -self.suspension["toe_rr"]])
         #todo add rear wheel steering
        state.delta_tire_rad = np.radians(np.array([state.delta, state.delta, 0, 0]))
        SA = -betaT + state.delta_tire_rad + np.radians(static_toe)

        return SA

    def transform_to_vehicle_system(self, tire_forces, SA, state):
        # Mirror LHS Tire forces
        #tire_forces["FyT"] = tire_forces["FyT"]
        #tire_forces["Mz"]  =  tire_forces["Mz"]
        # todo: Limit Fx to Max Engine Torque and get engine speed
        #     if(self.PowerTrain.diff_setting ==1)
        #         [FxT_mot, nmot] = self.limitFx(vx, SR, FxT, self.tire.r_0) 
        #     else
        #         [FxT_mot, nmot] = self.limitFx_nodiff(vx, SR, FxT, self.tire.r_0) 
        #     end
        #[FxT_mot, nmot] = self.limitFx(vx, SR, FxT, self.tire.r_0) 
        #[FxF_brake,FxR_brake] = self.limitFx_braking(FxT) 
         
        #FxT(FxT(:,1)<0,:) = [FxF_brake,FxF_brake, FxR_brake,FxR_brake] 
        #FxT(FxT(:,3)>=0,:) = FxT_mot(FxT(:,3)>=0,:) 
        
        #ind = FxT>0 
        #FxT(ind(:,1),1)=0 
        #FxT(ind(:,2),2)=0 
       
        # Transform to Vehicle Coordinate System, todo: check if this equation is correct
        Fx = np.cos(state.delta_tire_rad) * tire_forces["FxT"] - np.sin(state.delta_tire_rad) * tire_forces["FyT"]
        Fy = np.sin(state.delta_tire_rad) * tire_forces["FxT"] + np.cos(state.delta_tire_rad) * tire_forces["FyT"]
        return Fx, Fy

    def get_residual_forces(self, ax, ay, state, get_state=False):
        # Calculate Slip Angle (SA)
        vx = state.vx
        psi_dot = ay / vx      # No convergence when psi_dot ~= 0!
        SA = self.get_slip_angle(psi_dot, state) # C. Patton p.42
        state.SA = SA
        IA = np.array([self.suspension["camber_fl"], self.suspension["camber_fr"], self.suspension["camber_rl"], self.suspension["camber_rr"]])
        
        # Save previous Fz
        #                     Fz_k = Fz 
        vehicle_forces = self.get_vehicle_forces(ax, ay, psi_dot, state)
        # Tire Model
        SA_tm = np.degrees(SA)
        state.SA_tm = SA_tm
        IA = IA*0 
        mu_track = (self.suspension["mu_track_x"], self.suspension["mu_track_y"])
        tire_forces = self.tire.gForces(vehicle_forces["Fz"],SA_tm,state.SR,IA,mu_track, self.suspension["tire_temp"]) 
        #print(f"tire forces: {tire_forces}")
        # Yaw Accelerations
        lx = np.array([-self.general["track_f"]/2,self.general["track_f"]/2,-self.general["track_r"]/2,self.general["track_r"]/2])
        ly = np.array([(1-self.general["weight_distr"])*self.general["wheelbase"],
                      (1-self.general["weight_distr"])*self.general["wheelbase"],
                      -self.general["weight_distr"]*self.general["wheelbase"],
                      -self.general["weight_distr"]*self.general["wheelbase"]])
        #self.general["weight_distr"])*self.general["wheelbase"]
        #self.general["weight_distr"])*self.general["wheelbase"]] 
        Fx, Fy = self.transform_to_vehicle_system(tire_forces, SA, state)
        #print(f"after transformation: Fx: {Fx}, Fy: {Fy}")
        Mz = tire_forces["Mz"]
        # neglect Mz for now
        state.yaw_moment = sum(Fx*lx + Fy*ly)
        #print(f"yaw: Fx*lx: {Fx*lx}, Fy*ly: {Fy*ly}, Mz: {Mz}")
        state.psi_dotdot =  state.yaw_moment / self.general["inertia_z"]
        res_x = sum(Fx)- vehicle_forces["Fx"]
        res_y = sum(Fy)- vehicle_forces["Fy"] 
        if get_state:
            return state
        else:
            return np.array([res_x, res_y])
    

class state:
    result = []
    def __init__(self, vehicle, vx = 20, beta=1, SR=0, delta=20):
        self.vx = vx
        self.beta = beta
        self.SR = SR
        self.delta = delta
        self.warnings = []
        self.vehicle = vehicle

    def __str__(self) -> str:
        return f"State vx = {self.vx}, beta = {self.beta}, SR = {self.SR}, delta = {self.delta}, solution = {self.result}"

    def solve(self):
        res, info, ier, mesg = optimize.fsolve(lambda x: self.vehicle.get_residual_forces(x[0], x[1], self), [0,0], full_output=True)
        self.result = np.append(res, self.vx)
        self.ax = self.result[0]
        self.ay = self.result[1]
        self.warnings.append(mesg)
        return self