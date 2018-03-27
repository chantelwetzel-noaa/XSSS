


source("C:/Users/Chantel.Wetzel/Documents/GitHub/XSSS_AIS/R/Do_XSSS_AIS.R")

SSS.ais.fxn(filepath = "C:/My_Directory/.../", # parent directory to the file.name folder 
            control.name = "ENGL.ctl", 
            dat.name = "ENGL.dat",
            m.in = c(0.28, 0.27,0.40,0.40, 0), 
            h.in = c(0.87, 0.093, 0.20, 1, 1), 
            depl.in = c(0.60, 0.10, 0.01, 0.99, 1) )


SSS.ais.fxn(filepath = "C:/Users/Chantel.Wetzel/Desktop/pop_ais", # parent directory to the file.name folder 
            control.name = "pop.ctl" ,
            dat.name = "pop.dat",
            Niter = 2000, AIS.iter = 2000, final.Niter = 5000,
            m.in = c(0.054, 0.054,0.20,0.20, 0) ,
            h.in = c(0.50, 0.093, 0.20, 1, 1) ,
            depl.in = c(0.20, 0.10, 0.01, 0.99, 1))


filepath = "C:/Users/Chantel.Wetzel/Desktop/pop_ais" # parent directory to the file.name folder 
control.name = "pop.ctl" 
dat.name = "pop.dat"
m.in = c(0.054, 0.054,0.20,0.20, 0) 
h.in = c(0.50, 0.093, 0.20, 1, 1) 
depl.in = c(0.20, 0.10, 0.01, 0.99, 1)
Niter = 50 
AIS.iter = 50
final.Niter = 400

