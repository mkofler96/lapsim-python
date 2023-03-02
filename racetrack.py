import numpy as np
import scipy

class racetrack:
    #RACETRACK will hold all necessary geometric info for a laptime sim.
    #   first version will use x,y,(z) information in meters to create a
    #   tracks centerline and use a fixed track width to add left and right
    #   to define track limits. Further version will differ in this part!
    #   Also no information about banked corners will be used at this time.
    #
    #   Once the data is stored in this class the curvature k in each point
    #   is evaluated and stored for further use in the laptime simulation.
    
    # 17/06/22 Markus Ager
    # 18/04/16 Michael Kofler, added Dist and alpha

        # rawx
        # rawy
        # rawz
        # k
        # u,v
        # vlim
        # x,y,z
        # dir,dirN
        # r
        # local_v_min
        # dist
        # gradebeta #[rad], road grade
        # name        #name of track
        # alpha
    
    def __init__(self, xyz_csv_file):
        
        track_data = np.loadtxt(xyz_csv_file, delimiter=",", dtype=float)
        # # store raw data
        self.rawx = track_data[:,0]
        self.rawy = track_data[:,1]
        self.r = track_data[:,[0,1]]
        # self.rawy = rawy
        # self.rawz = rawz

        # # initialize k
        self.k = self.get_k()
        self.u = self.get_u()
        self.dist = np.cumsum(self.u, axis=0)
        self.corners, _ = scipy.signal.find_peaks(self.k)
        # # update xyz
        # self = update_xyz(self);
        # # update k
        # self = update_k(self);
        # # update dist
        # self = updatedist(self);

    
#     def update_xyz(self):            
#         # kk x 3 matrix
#         new_xyz = self.r' + self.dirN'.*[self.v,self.v,self.v];
        
#         self.x = new_xyz(:,1);
#         self.y = new_xyz(:,2);
#         self.z = new_xyz(:,3);
        
#         self.u = sqrt(gradient(self.x(:)).**2+gradient(self.y(:)).**2+gradient(self.z(:)).**2);
#         self.u(isnan(self.u)) = 0;
    def get_u(self):
        return np.linalg.norm(np.gradient(self.r, axis=0), axis=1)
    
    def get_k(self):        
         c_dot = np.gradient(self.r, axis=0)
         c_dotdot = np.gradient(np.gradient(self.r, axis=0), axis=0)
         int1 = np.cross(c_dot, c_dotdot)
         #take euclidean norm
         int2 = np.linalg.norm(np.gradient(self.r, axis=0), axis=1)**3
        
         return int1/int2;          

    
#     function self = update_v(self,v)
#         if (size(self.v) ~= size(v))
#             error('size of v does not match stored v''size');

#         self.v = v;
#         self = update_xyz(self);
#         self = update_r(self);
#         self = update_xyz(self);

    
#     function self = update_gradebeta(self)
#         # self = update_u(self)
        
#         s2D = sqrt(gradient(self.x).**2 + gradient(self.y).**2);
#         s3D = self.u;
        
#         deltah = sqrt(s3D.**2 - s2D.**2).*sign(gradient(self.z));
        
#         self.gradebeta = atan2(deltah,s2D);

    
#     function self = updatedist(self)
        
#         self.dist = cumtrapz(self.u);
    

    def find_local_k_peaks(self):
        return 

#     function self = getalpha(self)
#         self.alpha = atan(diff(self.y)./diff(self.x));
#         self.alpha = [self.alpha; self.alpha(end)];

