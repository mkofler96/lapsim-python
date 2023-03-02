import configparser
import numpy as np

class TMSimple:
    def __init__(self, tire_parameter_file):
        #todo read file
        params = configparser.ConfigParser()
        params.read(tire_parameter_file)
        params = params["tire_parameters"]
        self.fz_nom = float(params["fz_nom"]) 

        self.a1x = 2*float(params["fxmx_fzn"]) - 0.5*float(params["fxmx_2fzn"]) 
        self.a2x =  -float(params["fxmx_fzn"]) + 0.5*float(params["fxmx_2fzn"]) 
        self.b1x = 2*float(params["dfx0_fzn"]) - 0.5*float(params["dfx0_2fzn"]) 
        self.b2x =  -float(params["dfx0_fzn"]) + 0.5*float(params["dfx0_2fzn"]) 
        self.c1x = 2*float(params["fxin_fzn"]) - 0.5*float(params["fxin_2fzn"]) 
        self.c2x =  -float(params["fxin_fzn"]) + 0.5*float(params["fxin_2fzn"]) 
            
        self.a1y = 2*float(params["fymx_fzn"]) - 0.5*float(params["fymx_2fzn"])
        self.a2y =  -float(params["fymx_fzn"]) + 0.5*float(params["fymx_2fzn"]) 
        self.b1y = 2*float(params["dfy0_fzn"]) - 0.5*float(params["dfy0_2fzn"]) 
        self.b2y =  -float(params["dfy0_fzn"]) + 0.5*float(params["dfy0_2fzn"]) 
        self.c1y = 2*float(params["fyin_fzn"]) - 0.5*float(params["fyin_2fzn"]) 
        self.c2y =  -float(params["fyin_fzn"]) + 0.5*float(params["fyin_2fzn"]) 
        
        # I don't know which values are correct
        #a1x = (fxmx_2fzn-fxmx_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #a2x = fxmx_fzn - a1x 
        #b1x = (dfx0_2fzn-dfx0_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #b2x = dfx0_fzn - b1x 
        #c1x = (fxin_2fzn-fxin_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #c2x = fxin_fzn - c1x 

        #a1y = (fymx_2fzn-fymx_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #a2y = fymx_fzn - a1y 
        #b1y = (dfy0_2fzn-dfy0_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #b2y = dfy0_fzn - b1y 
        #c1y = (fyin_2fzn-fyin_fzn*((fz_2nom)**2 / (fz_nom)**2)) / ((fz_2nom/fz_nom)-((fz_2nom)**2 / (fz_nom)**2)) 
        #c2y = fyin_fzn - c1y 

        #MZ
        self.nL0_FzNom       = float(params["nL0_FzNom"])
        self.nL0_2FzNom      = float(params["nL0_2FzNom"])     

        self.ts_FzNom        = float(params["ts_FzNom"]) 
        self.ts_2FzNom       = float(params["ts_2FzNom"])    

        self.tinf_FzNom      = float(params["tinf_FzNom"])  
        self.tinf_2FzNom     = float(params["tinf_2FzNom"])   
        
        #
        self.r_0             = float(params["r_0"])
        self.c_v_FzNom       = float(params["c_v_FzNom"])   
        self.c_v_2FzNom      = float(params["c_v_2FzNom"]) 
        
        self.mu_x0           = float(params["mu_x0"]) 
        self.mu_y0           = float(params["mu_y0"]) 
        
        # camber 
        self.Sg4              = float(params["Sg4"])  # camber at 4° (usually 0.57~)
        
        # temperature
        self.K1_t = float(params["K1_t"])
        self.g_t = float(params["g_t"])
        self.B1_t = float(params["B1_t"]) 
        self.A1_t = float(params["A1_t"])
        self.k_t = float(params["k_t"])
        self.Ymax = float(params["Ymax"]) 
        self.Ymax_temp0 = float(params["Ymax_temp0"]) 
        self.T0 = float(params["T0"])
        self.Tnom = float(params["Tnom"]) 
    
    def gForces(self,Fz,SA,SR,IA,mu,T): 
        # todo: assert length is equal   
        fy_mx_t = self.Ymax_temp0 + (self.K1_t*np.sin(self.g_t*self.B1_t*(1-np.exp((-abs(T-self.T0))/self.A1_t)))*np.sign(T-self.T0)) 
                    
        fxmx = (self.a1x * Fz/self.fz_nom + self.a2x * (Fz/self.fz_nom)**2)*(mu[0]/self.mu_x0) 
        dfx0 = self.b1x * Fz/self.fz_nom + self.b2x * (Fz/self.fz_nom)**2 
        fxin = (self.c1x * Fz/self.fz_nom + self.c2x * (Fz/self.fz_nom)**2)*(mu[0]/self.mu_x0) 
        fymx = (self.a1y * Fz/self.fz_nom + self.a2y * (Fz/self.fz_nom)**2)*(mu[1]/self.mu_y0)*fy_mx_t/self.Ymax 
        dfy0 = (self.b1y * Fz/self.fz_nom + self.b2y * (Fz/self.fz_nom)**2)*(self.Tnom*self.k_t/T+1-self.k_t) 
        fyin =(self.c1y * Fz/self.fz_nom + self.c2y * (Fz/self.fz_nom)**2)*(mu[1]/self.mu_y0) 
        
        
        Sg = -self.Sg4/4*IA 
#             Sg = deg2rad(Sg) 
        
        if (np.logical_or(fxin>fxmx, fxin<-fxmx)).any():
            return {"FxT": np.zeros_like(SA), "FyT": np.zeros_like(SA), "Mz": np.zeros_like(SA)}

        if (np.logical_or(fyin>fymx, fyin<-fymx)).any():
            return {"FxT": np.zeros_like(SA), "FyT": np.zeros_like(SA), "Mz": np.zeros_like(SA)}

        Kx = fxmx 
        Bx = np.pi-np.real(np.arcsin(fxin/fxmx)) 
        Ax = Kx * Bx / dfx0 
        
        Ky = fymx 
        By = np.pi-np.real(np.arcsin(fyin/fymx)) 
        Ay = Ky * By / dfy0      
        
        # Weighing Factor
        G = Ay * Kx * Bx / Ax / Ky / By 

        slx     = SR 
        sly     = SA/(G+1e-5)  # wieso 1e-5? 
        
        beta    = np.arctan2(sly,slx) 
        sc      = np.sqrt(slx**2+sly**2) 
        
        # Pure Fx
        Fx = Kx * np.sin(Bx*(1-np.exp(-abs(sc)/Ax))*np.sign(sc))         #maybe sc = SR
        
        # Pure Fy without camber
        Fy0 = Ky * np.sin(By*(1-np.exp(-abs(G*sc)/Ay))*np.sign(G*sc))    # maybe G*sc = SA
        
        # Pure Fy with camber
        Fy = Ky * np.sin(By*(1-np.exp(-abs(G*sc+Sg)/Ay))*np.sign(G*sc+Sg)) 
        
        # delta Fy pure
        dFy = Fy-Fy0 

        
        # Combined Force
        F  = 0.5 * (abs(Fx)+abs(Fy0)+(abs(Fx)-abs(Fy0))*np.cos(2*beta)) 
        
        FxT = F * np.cos(beta) 
        FyT = F * np.sin(beta)+dFy 
        
        # MZ
        nL0 = (self.nL0_2FzNom - self.nL0_FzNom)/self.fz_nom * (Fz-self.fz_nom) + self.nL0_FzNom 
        ts = (self.ts_2FzNom - self.ts_FzNom)/self.fz_nom * (Fz-self.fz_nom) + self.ts_FzNom 
        tinf = (self.tinf_2FzNom - self.tinf_FzNom)/self.fz_nom * (Fz-self.fz_nom) + self.tinf_FzNom 
        
        c_v = (self.c_v_2FzNom - self.c_v_FzNom)/self.fz_nom * (Fz-self.fz_nom) + self.c_v_FzNom 
        
        rstat = self.r_0 - Fz / c_v 
        L = 2*np.sqrt(self.r_0**2 - rstat**2) 
        
        temp2 = -nL0 * (abs(SA)-ts)/ ts *((tinf-abs(SA))/(tinf-ts))**2   
        
        nL = -0.5*(np.sign(abs(SA)-ts)-1)* nL0 * (1-abs(SA)/ts) -0.5*(np.sign(abs(SA)-tinf)-1)*0.5*(np.sign(abs(SA)-ts)+1)*temp2 
            
        Mz = -Fy * nL *L
        forces = {"FxT": FxT, "FyT": FyT, "Mz": Mz}
        return forces