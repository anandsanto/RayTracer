##########################################################################
#
#
#
#               TODO : ADD YOUR DOCUMENTATION HERE
#
#
#
#
##########################################################################


import math
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat



#To add the given vectors a and b
def addvector(a,b):
    r=[]
    for (i,j) in zip(a,b):
        r.append(i+j)
    return r

#To perform scalar multiplication where scalar is b and vector is a
def scalmul(a,b):
    r=[]
    for i in a:
        r.append(i*b)
    return r

# To get the modulus of a 3-D vector a
def modvector(a):
    b = (a[0]**2)+(a[1]**2)+(a[2]**2)
    return(b**0.5)

# To get dot-products of 2 vectors a and b
def dotprod(a,b):
    sum = 0
    for (i,j) in zip(a,b):
        sum = sum + (i*j)
    return sum

# To get unit vector of vector a
def unitvec(a):
    mod = modvector(a)
    mod = 1/mod
    return(scalmul(a,mod))

# Implement Snell's law from incident vector, normal vector and refractive indices of mediums 1 and 2
# Inputs : Incident vector, Normal vector, Refarctive index 1, Refractive index 2, curvature
# Output : Direction of refracted ray

def snell(i,n,n1,n2,coc):

    uniti = unitvec(i)
    unitn = unitvec(n)
    #print("unit vectors are",uniti,unitn)
    dotin = dotprod(uniti,unitn)
    #print("Dot prod of i and n is %f"%dotin)
    angi = math.acos(dotin) # Angle of incidence
    if(angi>(math.pi/2)): angi = math.pi - angi
    sinangr = (n1 * math.sin(angi))/n2

    # Check for Total Internal Reflection
    if(sinangr>1):
        print("We experience Total Internal Reflection")
        #TODO Deal with TIR
        return None

    angr = math.asin(sinangr) #Angle of refraction
    b = n1/n2 # Ratio of refractive indices
    cosangr = math.sqrt(1-((sinangr)**2))

    # Component of refracted ray along incident ray
    trans = scalmul(uniti,b)
    term2 = b*(math.cos(angi)) - cosangr

    # Component of refracted ray along normal
    if(coc<0): term2_vec = scalmul(unitn,-term2)
    else: term2_vec = scalmul(unitn,term2)

    # Summed up refracted ray
    trans = addvector(trans, term2_vec)
    return trans

# A method to 3-D plot all the rays in the list Rays
def plot_rays(Rays):
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    axes = list(repeat(0, 100))
    #print Rays[0].vertices()
    count = 0
    for ray in Rays:
        iterator = 0
        axes[count] = fig.gca(projection='3d')
        points,dirs = ray.vertices()
        x = np.array([]).reshape(0,1)
        y = np.array([]).reshape(0,1)
        z = np.array([]).reshape(0,1)
        #print points
        for (i,j) in zip(points,dirs):
            if(j[0]==0 and j[1]==0 and j[2]==0): break
            try:
                k = points[iterator+1]
            except:
                break
            # lambd = (k[0] - i[0])/j[0]
            x = np.append(x,np.linspace(i[0],k[0],3))
            y = np.append(y,np.linspace(i[1],k[1],3))
            z = np.append(z,np.linspace(i[2],k[2],3))
            iterator = iterator + 1
        #print(x.shape,y.shape,z.shape)
        axes[count].plot(z, y, x)
        #axes[count].legend()
        count = count + 1

    plt.show()

# Class for Ray
class Ray:

    # To maintain history of points and directions
    prevpoints = np.array([]).reshape(0,3)
    prevdirs = np.array([]).reshape(0,3)
    terminated = False

    #Constructor definition
    def __init__(self, start, direction):
        self.curpoint = np.array(start)
        self.curdir = np.array(direction)
        self.prevdirs=np.append(self.prevdirs,[self.curdir],axis=0)
        self.prevpoints=np.append(self.prevpoints,[self.curpoint],axis=0)

    # Current Point Getter method
    def get_point(self):
        return self.curpoint

    # Current Direction Getter method
    def get_dir(self):
        return self.curdir

    # Update the current point and direction method
    def append(self, point, direction):
        if(self.terminated == False):
            self.curpoint = np.array(point)
            self.curdir = np.array(direction)
            self.prevdirs=np.append(self.prevdirs,[self.curdir],axis=0)
            self.prevpoints=np.append(self.prevpoints,[self.curpoint],axis=0)
        else:
            print("Sorry, the ray was already terminated..!")
            self.prevdirs=np.append(self.prevdirs,[self.curdir],axis=0)
            self.prevpoints=np.append(self.prevpoints,[self.curpoint],axis=0)

    # To get all the points and directions in the history
    def vertices(self):
        return self.prevpoints,self.prevdirs

    def set_terminated(self):
        self.terminated = True


# A general class for Optical Elements - All the optical elements should be defined as children to this class
class OpticalElement:
    def propagate_ray(self, ray):

        # Mandate the implementation of propagate_ray in child classes

        raise NotImplementedError()

# Spherical refraction Element - subclass of Optical Element
class SphericalRefraction(OpticalElement):

    # Constructor to get 5 inputs and instantiation

    def __init__(self, z0, curvature, n1, n2, aperadius):

        # Pt where lens meets z axis
        self.z0 = z0
        #Signed Inversion of Radius
        self.curvature = curvature
        #Refractive indices
        self.n1 = n1
        self.n2 = n2
        #Aperture radius
        self.aperadius = aperadius

    # Intercept finding method

    def intercept(self,ray):

        #Case where the refracting surface is planar
        if(self.curvature==0):
            posz = ray.get_point()[2]
            dirz = ray.get_dir()[2]
            lamb = (self.z0 - posz)/dirz
            term2 = scalmul(ray.get_dir(),lamb)
            intercept = addvector(term2,ray.get_point())
            levelpoint=[-1*ray.get_point()[0],-1*ray.get_point()[1],-1*self.z0]
            diffvec = addvector(levelpoint,intercept)
            mod = modvector(diffvec)
            if(mod>self.aperadius):
                print "The ray doesnt intersect the medium"
                return None
            else: return intercept

        else:
            # Vector for centre of curvature (COC) from origin
            rz = (self.z0 + (1/self.curvature))
            # Vector between COC and current point of ray
            r = [ray.get_point()[0],ray.get_point()[1],(ray.get_point()[2] - rz)]
            #Unit vector in ray direction
            k=[]
            for i in ray.get_dir():
                k.append(i/modvector(ray.get_dir()))
            modr = modvector(r)
            dotrk = dotprod(r,k)

            #To find distance between intercept and current point of ray
            if(((dotrk**2) - ((modr**2)-((1/self.curvature)**2))) < 0): return None # Root of a negative number : Invalid
            l1 = (-1*dotrk) + ((dotrk**2) - ((modr**2)-((1/self.curvature)**2)))**(0.5)
            l2 = (-1*dotrk) - ((dotrk**2) - ((modr**2)-((1/self.curvature)**2)))**(0.5)
            #print modr,l1,l2, ((dotrk**2) - ((modr**2)-((1/self.curvature)**2)))
            if(l1>=0 and l2>=0): self.l = min(l1,l2)
            elif((l1<0 and l2>=0)): self.l = l2
            elif((l1>=0 and l2<0)): self.l = l1
            else:
                print("Invalid Distance to the intercept from current point..!")
                return None # Invalid Solution

            #To find the intercept
            intercept = []
            for (i,j) in zip(ray.get_point(),ray.get_dir()):
                intercept.append(i + (self.l * j))
            #print("Intercept is ",intercept)
            if(abs(modvector([intercept[0],intercept[1],0]))>self.aperadius):
                print("Please increase the Aperture radius")
                return None
            return intercept

    def propagate_ray(self, ray):
        intercept = self.intercept(ray)
        if(intercept == None):
            print ("There is no Intercept found..!")
            # Dummy append at a different pt in same direction to maintain array length
            ray.append(addvector(ray.get_point(),ray.get_dir()),ray.get_dir())
            return
        #Check if the refracting surface is planar
        if(self.curvature != 0):
            coc = [0,0,(self.z0 + (1/self.curvature))]
            n = addvector(intercept,scalmul(coc,-1))
        # Else it is circular
        else:
            if(self.z0 != 0): n = [0,0,-self.z0]
            else: n = [0,0,-ray.get_point()[2]]
        new_dir = snell(ray.get_dir(),n,self.n1,self.n2,self.curvature)

        #If ray could not be refracted, terminate the ray!
        if(new_dir == None):
            ray.append(intercept,[0,0,0])
            ray.set_terminated()
        ray.append(intercept,new_dir)
        return

# Output plane - OpticalElement subclass to terminate
class OutputPlane(OpticalElement):

    def __init__(self, z0):
        self.z0 = z0

    def intercept(self,ray):

        posz = ray.get_point()[2]
        dirz = ray.get_dir()[2]
        lamb = (self.z0 - posz)/dirz
        term2 = scalmul(ray.get_dir(),lamb)
        intercept = addvector(term2,ray.get_point())
        levelpoint=[-1*ray.get_point()[0],-1*ray.get_point()[1],-1*self.z0]
        diffvec = addvector(levelpoint,intercept)
        mod = modvector(diffvec)
        return intercept

    def propagate_ray(self, ray):
        #print("Ray propagated to output")
        ray.append(self.intercept(ray),[0,0,0])
        ray.set_terminated()
        return

### MAIN PROGRAM ###
if __name__ == "__main__":


    Rays = [] # List registry containg the rays that will be plotted

    # To create a circular bundle of ray9s
    theta = np.linspace(0,2*math.pi,70)

    # A convex surface at z = 100
    b = SphericalRefraction(100,0.03,1,1.5,11)

    # Output plane at z = 350
    c = OutputPlane(350)

    # Varying radius to create a spiral
    r = np.linspace(0.1,14,70)

    # A concave surface at z=160
    d = SphericalRefraction(130,-0.03,1.5,1,11)
    for i in range(69):
        a = Ray([math.cos(theta[i]),math.sin(theta[i]),0],[0,0,1]) # Ray instantiation
        b.propagate_ray(a) # Pass the ray through surface b
        d.propagate_ray(a) # Pass the ray through surface d
        c.propagate_ray(a) # Output surface and terminate the waves
        Rays.append(a) # Each iteration, the ray created is appended to Rays

    plot_rays(Rays)
