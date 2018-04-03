
# Generic call to function 
# Only the required inputs are specified here
SSS.ais.fxn(filepath = "C:/My_Directory/.../", # parent directory to the file.name folder 
            control.name = "control_file.ctl", # name of the model control file
            dat.name = "data_file.dat", # name of the model data file
            m.in = c(0.28, 0.27,0.40,0.40, 0), # c(female M median value, male M median value, female M sd, male M sd, M equal between sexes (0=No, 1=Yes))
            h.in = c(0.87, 0.093, 0.20, 1, 1), # c(steepness, mean, steepness sd, lower bound, upper bound, distribution (1= beta))
            depl.in = c(0.60, 0.10, 0.01, 0.99, 1) #c(depletion in select year, sd, lower bound, upper bound, distribution (1 = beta, 2 = lognormal))
            )


# Function call with all inputs specified
SSS.ais.fxn(filepath = "C:/My_Directory/.../", 
            control.name = "control_file.ctl", 
            dat.name = "data_file.dat", 
            tantalus = FALSE, 
            read.seed = FALSE, # Option to allow the user to predefine seeds
            entropy.level = 0.92, # Level of entropy to meet before finishing AIS
            Niter = 2000, # Initial model run for MCMC
            AIS.iter = 2000, # Interum run size for AIS
            final.Niter = 5000, # Final model sample size
            m.in = c(0.28, 0.27,0.40,0.40, 0), 
            h.in = c(0.87, 0.093, 0.20, 1, 1), 
            depl.in = c(0.60, 0.10, 0.01, 0.99, 1), 
            fmsy.m.in = NULL, # Not yet fully implemented
            bmsy.b0.in = NULL # Not yet fully implemented
            ) 


