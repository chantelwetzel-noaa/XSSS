

SSS.ais.fxn(filepath = "C:/My_Directory/.../", # parent directory to the file.name folder 
            control.name = "control_file.ctl", 
            dat.name = "data_file.dat",
            m.in = c(0.28, 0.27,0.40,0.40, 0), 
            h.in = c(0.87, 0.093, 0.20, 1, 1), 
            depl.in = c(0.60, 0.10, 0.01, 0.99, 1) )


SSS.ais.fxn(filepath = "C:/Users/Chantel.Wetzel/Desktop/pop_ais", # parent directory to the file.name folder 
            control.name = "pop.ctl" ,
            dat.name = "pop.dat",
            entropy.level = 0.80,
            Niter = 200, AIS.iter = 200, final.Niter = 200,
            m.in = c(0.054, 0.054,0.20,0.20, 0) ,
            h.in = c(0.50, 0.093, 0.20, 1, 1) ,
            depl.in = c(0.20, 0.10, 0.01, 0.99, 1))


filepath = "C:/Users/Chantel.Wetzel/Desktop/pop_ais" # parent directory to the file.name folder 
control.name = "pop.ctl" 
dat.name = "pop.dat"
m.in = c(0.054, 0.054,0.20,0.20, 0) 
h.in = c(0.50, 0.093, 0.20, 1, 1) 
depl.in = c(0.20, 0.10, 0.01, 0.99, 1)
Niter = 200 
AIS.iter = 200
final.Niter = 200
tantalus=FALSE 
read.seed = FALSE
entropy.level = 0.80
fmsy.m.in = NULL
bmsy.b0.in = NULL

