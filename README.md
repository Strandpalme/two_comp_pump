# two_comp_pump
two compatment model of the leech t cell with cell intrinsic plasticity

# Setup
Safe the folder "two_comp_pump" on your computer. Change the path in line 
10 of "model_demo_sc.m" or in line 2 of "model_demo.mlx" to the corresponding
path on your system. 

# Models
Two versions of the two compartment model can be found in the "models" 
directory. 
    - TcellDoublePump.m     
                            This is the model as it was used in my master
                            thesis. It's a good idea to keep it unchanged 
                            as a reference.

    - TcellDoublePumpFitting.m
                            This version of the model gets an array of pump
                            parameter values to change properties of the
                            Na+/K+ pump. This is included as an example how
                            parameters can be changed in the model.

# tools
This directory includes two basic analysis scripts, based on the analysis
in my master thesis. They are added to the MATLAB path in both demo scripts.