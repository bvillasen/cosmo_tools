import numpy as np
import matplotlib.pyplot as plt

OmegaM_h = np.array([ [ 0.2730324, 0.706981 ],
                      [ 0.2808033, 0.70019  ],
                      [ 0.2905362, 0.691785 ],
                      [ 0.300962,  0.683534 ],
                      [ 0.3110418, 0.675953 ] ])
                      
OmegaM_s8 = np.array([ [ 0.282131, 0.79138 ],
                       [ 0.287195, 0.800129 ],
                       [ 0.292056, 0.808553 ], 
                       [ 0.297324, 0.816653 ],
                       [ 0.301376, 0.822809 ],
                       [ 0.305224, 0.829937 ],
                       [ 0.309886, 0.836092 ],
                       [ 0.313937, 0.842896 ],
                       [ 0.317786, 0.849052 ],
                       [ 0.320623, 0.853588 ] ])
                       
OmegaM_s8_line_params = np.polyfit( OmegaM_s8[:,0], OmegaM_s8[:,1], 1)

OmegaM_vals = OmegaM_h[:,0]
h_vals = OmegaM_h[:,1]
s8_vals = OmegaM_s8_line_params[1] + OmegaM_vals*OmegaM_s8_line_params[0]

alpha = 0.45
sigma8_vals = s8_vals / ( OmegaM_vals / 0.3 )**alpha
