# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:38:02 2016

@author: YinCheang
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import matplotlib.tri as tri

class CM:
    
    inlet_Ma   = 1.01
    inlet_leny = 1.0
    outlet_Ma  = 2.0
    
    inlet_n = 10    
    pointx = []
    pointy = []
    pointv = []
    pointa = []
    
    def wall_shape(self,x,a):
        
        f = a*((0.0010*x**4)-(0.0116*x**3)+(0.0133*x**2)+(0.2355*x)+0.9998)
        
        return f
        
    def inipoint(self):
        
        f = lambda x : (0.0040*x**3)-(3*0.0116*x**2)+(2*0.0133*x)+(0.2355)  
        
        temp = f(0)
#        print temp
        temp = temp/self.inlet_n
        
        delta = self.inlet_leny/self.inlet_n
        
        f1 = lambda y,r : np.sqrt(r**2 - (y-0.5)**2)
        
        for i in range(self.inlet_n+1):
            print i
            self.pointy.append((i)*delta )
            self.pointx.append(0.0)
            
            self.pointv.append(self.cal_v(self.inlet_Ma))
            self.pointa.append(i*(temp))
            
#        plt.plot(self.pointa)
        print self.pointx
        plt.plot(self.pointx,self.pointy,'r*')
        plt.show()
#        print self.pointa
        pause = raw_input("WTF")
        return 0
        
    def cal_RIGHTwallpoint(self,v1,s1,x1):
        
#        f = lambda x : ((0.0040*x**3)-(3*0.0116*x**2)+(2*0.0133*x)+(0.2355))*-1  
        
        Q     = v1+s1        
        swall = 0.0 #np.arctan(f(x1))
        vwall = Q-swall 
        
        return vwall,swall
        
    def cal_LEFTwallpoint(self,v1,s1,x1):
        
        f = lambda x : (0.0040*x**3)-(3*0.0116*x**2)+(2*0.0133*x)+(0.2355)        
        
        R     = v1-s1        
        print x1 
        swall = np.arctan(f(x1))
#        print
#        print "wall"
#        print swall
#        print
#        pause = raw_input("HAHAHA")
        vwall = R+swall 
        
        print
        print vwall*180./np.pi
        print R
        print swall*180./np.pi
#        pause = raw_input("76")        
        
        return vwall,swall
        
    def intepoint_posi(self,x1,y1,x2,y2,v1,v2,s1,s2):
        
        M1 = self.cal_Mayerfun(v1)
        M2 = self.cal_Mayerfun(v2)
        
        MaAng1 = self.cal_MaAng(M1)
        MaAng2 = self.cal_MaAng(M2)
        
        m1 = np.tan(MaAng1+s1)
        m2 = np.tan(s2-MaAng2)
        

        
        a = np.array([ [1.0,-m1]
                      ,[1.0,-m2] ] )
        b = np.array([ y1-m1*x1
                      ,y2-m2*x2] )
        
        a = np.linalg.inv(a)
        ans = np.dot(a,b)

#        ans = np.linalg.solve(a,b)
        y   = ans[0]
        x   = ans[1]
        
        print M1 , M2
        print m1 , MaAng1*180/np.pi  , s1*180/np.pi 
        print m2 , MaAng2*180/np.pi  , s2*180/np.pi 
        print x , y
          
#        print ans[1]
      
#        pause = raw_input("WHAT")    
    
        return x,y
        
    def wallLEFTpoint_posi(self,x1,y1,v1,s1):
        
        #x = wallfunction(y)
        
        M1      = self.cal_Mayerfun(v1)
        MaAng1  = self.cal_MaAng(M1)
        m1      = np.tan(s1+MaAng1)
       
#        wallfun = lambda y : 1.0
        f = lambda x : (0.0010*x**4)-(0.0116*x**3)+(0.0133*x**2)+(0.2355*x)+0.9998 - (m1*x + (y1-m1*x1))
       
#        xans = (1.0 - (y1 - m1*x1 ))/m1
        
        xans = sp.fsolve(f,x1)
        yans = self.wall_shape(xans,1.0)
#        yans = wallfun(xans)
        
        print x1,y1
        print xans,yans, m1 ,x1
#        plt.plot(xans,yans,'b*')
#        plt.show()
#        pause =  raw_input("124")

        return  xans,yans
    
        
    def wallRIGHTpoint_posi(self,x1,y1,v1,s1):
        
        #x = wallfunction(y)
#        print x1 ,y1
        
        M1      = self.cal_Mayerfun(v1)
        MaAng1  = self.cal_MaAng(M1)
        m1      = np.tan(s1-MaAng1)
       
        wallfun = lambda y : 0.0
        
        xans = (0.0 - (y1 - m1*x1 ))/m1
        yans = wallfun(xans)
        
#        print xans,yans, m1 

        return  xans,yans    
    
    def cal_MaAng(self,M):
    
        MaAng = np.arcsin(1./M)
        
        return MaAng
        
    def cal_v(self,M):
        
        f = np.sqrt((2.4/0.4))*np.arctan(np.sqrt(0.4/2.4)*np.sqrt(M**2 -1)) - np.arctan(np.sqrt(M**2 -1))
        
        return f
    
    def cal_Mayerfun(self,v):

        Guess = np.array([1.5])        
        
        f = lambda M : v-(np.sqrt((2.4/0.4))*np.arctan(np.sqrt(0.4/2.4)*np.sqrt(M**2 -1))- np.arctan(np.sqrt(M**2 -1)))
                        
        Mans = sp.fsolve(f,Guess)   
        
        return Mans
        
    def calculation(self,leng,toleng):
        
#        print "toleng" , toleng - leng-1
        if leng == 9 :
            toleng = toleng - 1
        
        cco = 0
        v = lambda v1,v2,s1,s2 : 0.5*(v1+v2)+0.5*(s1-s2)
        a = lambda v1,v2,s1,s2 : 0.5*(v1-v2)+0.5*(s1+s2)
        
#        v = lambda v1,v2,s1,s2 : 0.5*(v2+v1)+0.5*(s2-s1)
#        a = lambda v1,v2,s1,s2 : 0.5*(v2-v1)+0.5*(s2+s1)
        
        for i in range(leng):
            i = toleng-leng+i-1
            
            x,y = self.intepoint_posi(self.pointx[i],self.pointy[i],self.pointx[i+1],self.pointy[i+1], 
                                      self.pointv[i],self.pointv[i+1],self.pointa[i],self.pointa[i+1])
                                      
            self.pointx.append(float(x))
            self.pointy.append(float(y))
            self.pointv.append(float(v(self.pointv[i+1],self.pointv[i],self.pointa[i+1],self.pointa[i])))
            self.pointa.append(float(a(self.pointv[i+1],self.pointv[i],self.pointa[i+1],self.pointa[i])))
            

            
#            print self.pointa[len(self.pointa)-1]
            
#            plt.plot(self.pointx,self.pointy,'r*')  
#            plt.show()
#            pause = raw_input("intepoint")
 
#            xk = np.linspace(0,4.836928805312580)
#            yk = self.wall_shape(xk)
#            plt.plot(xk,yk)   
#            plt.plot(xk,-yk)
            cco +=1
#            print "-----181---- ", i 
#            pause = raw_input("WAHAHHAA")
        return cco
        
    def cal_pros(self):
        
        j   = 0
        
        self.inipoint()
        
        cco = len(self.pointx)-1
#        print cco
        
#        print self.pointv 
#        print self.pointa       
        
#        plt.plot(self.pointx,self.pointy,'r*')
        
        
        while (self.pointx[j-1]) < 4.836928805312580:


#            print cco    
            cco = self.calculation(cco,len(self.pointx))
            
            print "----206---- ",cco
#            pause = raw_input("HAHAHA")
            
            
            if cco == 9 :
                i   = len(self.pointx)-1-(cco+1)
#                print i , len(self.pointx) , cco
#                pause = raw_input("210")
                x,y = self.wallLEFTpoint_posi(self.pointx[i],self.pointy[i],self.pointv[i],self.pointa[i])
                v,s = self.cal_LEFTwallpoint(self.pointv[i],self.pointa[i],x)
                self.pointx.append(float(x) )
                self.pointy.append(float(y) )
                self.pointv.append(float(v) )
                self.pointa.append(float(s) )   

                print 
                print x , y , " WALL "
                print "-------v------"
                print v*180/np.pi
                print "-------s------"   
                print s*180/np.pi 
                print
            
                
            if cco == 10 :
                i   = len(self.pointx)-cco   
                x,y = self.wallRIGHTpoint_posi(self.pointx[i],self.pointy[i],self.pointv[i],self.pointa[i])
                v,s = self.cal_RIGHTwallpoint(self.pointv[i],self.pointa[i],x)
                self.pointx.append(float(x))
                self.pointy.append(float(y))
                self.pointv.append(float(v))
                self.pointa.append(float(s))

            
            if cco == 10 :
                cco -= 1
            else:
                cco += 1
                
#            print "237---- ", cco
                
                
#            print v,s
            
            
#        i   = len(self.pointx) -11
#        print i            
#        print self.pointx[i] , self.pointy[i]
#        x,y = self.wallLEFTpoint_posi(self.pointx[i],self.pointy[i],self.pointv[i],self.pointa[i])
#        self.pointx.append(x)
#        self.pointy.append(y)
        
#            plt.plot(self.pointx,self.pointy,'r*') 
#            xk = np.linspace(0,4.836928805312580)
#            yk = self.wall_shape(xk,1.0)
#            plt.plot(xk,yk)   
#            plt.plot(xk,-yk)
#            plt.show()
            
#            pause = raw_input("wait")
            
#            pause = raw_input("Thanks")
            j = len(self.pointx)
#            print j ,j+cco-1, self.pointx[j-1]
        
        return 0
        
    
def main():
    
    qes = CM()
    qes.cal_pros()
    M = []
    
    plt.plot(qes.pointx,qes.pointy,'r*') 
    xk = np.linspace(0,4.836928805312580)
    yk = qes.wall_shape(xk,1.0)
    plt.plot(xk,yk)   
    plt.show()


#    z = np.zeros([len(qes.pointx),len(qes.pointy)])+100
    M = qes.pointv 

    for i in range(len(qes.pointv)):
        M[i] = float(qes.cal_Mayerfun(M[i]))
        
#        qes.pointv[i] = qes.pointv[i]*180.0/np.pi

#    plt.figure(1)

    print len(M) , len(qes.pointx)


    plt.plot(qes.pointx,M,'b*')   
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.tricontourf(qes.pointx,qes.pointy, M)

    
#    plt.plot(qes.pointx,qes.pointa)
#    plt.figure(2)
#    plt.plot(qes.pointx,qes.pointv)
    
#    x = np.linspace(0,4.836928805312580)
#    y = qes.wall_shape(x)
#    plt.plot(x,y)
#    plt.show()
    
    return 0
    
    

if __name__=="__main__":
    main()