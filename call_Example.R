


source("F:/NOAA/SSSais/code/XSSS_AIS/XSSS_AIS_v2.R")
SSS.ais.fxn(filepath = "C:/My_Directory/.../", # parent directory to the file.name folder 
            file.name = "EnglishSole" , # folder that contains the model files 
            control.name = "ENGL.ctl", 
            dat.name = "ENGL.dat",
            m.in = c(0.28, 0.27,0.40,0.40), 
            h.in = c(0.87, 0.093, 0.20, 1, 1), 
            depl.in = c(0.60, 0.10, 0.01, 0.99, 1), 
            start.m.equal = FALSE, 
            years = seq(1876,2012,1), 
            ofl_yrs = seq(2013,2024,1),
            depl.yr = 2000 )
