import numpy as np 


data_thermal_history_Hiss_2018 = {
'name': "Hiss et al. (2018)",
'z': np.array([ 2, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4 ]),
'T0': np.array([ 13721, 10927, 13334, 16281, 20036, 18371, 16244, 13439 ]),
'T0_sigma_plus': np.array([ 1694, 1961, 1206, 1940, 1416, 1087, 1153, 1542 ]),
'T0_sigma_minus':np.array([ 2152, 3843, 1530, 1601, 1507, 1388, 1135, 2318 ]),
'gamma': np.array([ 1.47, 1.67, 1.56, 1.38, 1.29, 1.12, 1.38, 1.31  ]),
'gamma_sigma_plus': np.array([ 0.12, 0.27, 0.12, 0.08, 0.07, 0.12, 0.13, 0.14 ]),
'gamma_sigma_minus':np.array([ 0.10, 0.14, 0.12, 0.10, 0.07, 0.06, 0.13, 0.10]),
}

data_thermal_history_Bolton_2014 = {
'name': "Bolton et al. (2014)",
'z': np.array([ 2.4,  ]),
'T0': np.array([ 10000 ]),
'T0_sigma_plus': np.array([ 3200 ]),
'T0_sigma_minus':np.array([ 2100 ]),
'gamma': np.array([ 1.54  ]),
'gamma_sigma_plus': np.array([ 0.11 ]),
'gamma_sigma_minus':np.array([ 0.11 ]),
}

data_thermal_history_Walther_2019 = {
'name': 'Walther et al. (2019)',
'z': np.array([ 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4 ]) ,
'T0':            np.array([ 0.768, 0.732, 1.014, 1.165, 1.234, 1.286, 1.289, 1.186, 1.404, 1.038, 1.205, 0.940, 0.890, 0.877, 0.533, 0.599 ])*1e4,
'T0_sigma_plus': np.array([ 0.369, 0.196, 0.250, 0.290, 0.193, 0.191, 0.182, 0.133, 0.165, 0.313, 0.229, 0.220, 0.093, 0.130, 0.122, 0.152 ])*1e4,
'T0_sigma_minus':np.array([ 0.218, 0.091, 0.150, 0.189, 0.139, 0.147, 0.144, 0.115, 0.157, 0.267, 0.194, 0.173, 0.073, 0.106, 0.091, 0.134 ])*1e4,
'gamma':             np.array([ 1.63, 1.88, 1.74, 1.63, 1.67, 1.78, 1.60, 1.75, 1.74, 1.69, 1.41, 1.27, 1.85, 1.84, 1.64, 1.54 ]),
'gamma_sigma_plus':  np.array([ 0.16, 0.20, 0.15, 0.16, 0.13, 0.11, 0.14, 0.11, 0.10, 0.14, 0.20, 0.24, 0.23, 0.23, 0.26, 0.29 ]),
'gamma_sigma_minus': np.array([ 0.25, 0.27, 0.21, 0.19, 0.15, 0.12, 0.16, 0.13, 0.11, 0.25, 0.23, 0.24, 0.33, 0.33, 0.32, 0.29 ])
}

data_thermal_history_Boera_2019 = {
'name': 'Boera et al. (2019)',
'z': np.array([ 4.2, 4.6, 5.0 ]),
'T0': np.array([ 8.13, 7.31, 7.37 ]) * 1e3,
'T0_sigma_plus': np.array([ 1.34, 1.35, 1.67 ]) * 1e3,
'T0_sigma_minus': np.array([0.97, 0.88, 1.39 ]) * 1e3,
'gamma': np.array([ 1.21, 1.29, 1.33 ]),
'gamma_sigma_plus': np.array([ 0.23, 0.19, 0.18 ]),
'gamma_sigma_minus': np.array([ 0.28, 0.26, 0.27 ]),

}


data_thermal_history_Gaikwad_2020a = {
'name': 'Gaikwad et al. (2020a)',
'z': np.array([ 5.4, 5.6, 5.8 ]),
'T0': np.array([ 11000, 10500, 12000 ]),
'T0_sigma_plus':  np.array([ 1600, 2100, 2200 ]) ,
'T0_sigma_minus': np.array([ 1600, 2100, 2200 ]) ,
'gamma': np.array([ 1.20, 1.28, 1.04 ]),
'gamma_sigma_plus': np.array([ 0.18, 0.19, 0.22 ]),
'gamma_sigma_minus': np.array([ 0.18, 0.19, 0.22]),
}


data_thermal_history_Gaikwad_2020b = {
'name': 'Gaikwad et al. (2020b)',
'z': np.array([ 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8 ]),
'T0': np.array([ 9500, 11000, 12750, 13500, 14750, 14750, 12750, 11250, 10250, 9250  ]),
'T0_sigma_plus':  np.array([ 1393, 1028, 1132, 1390, 1341, 1322, 1493, 1125, 1070, 876 ]) ,
'T0_sigma_minus': np.array([ 1393, 1028, 1132, 1390, 1341, 1322, 1493, 1125, 1070, 876 ]) ,
'gamma': np.array([ 1.5, 1.425, 1.325, 1.275, 1.250, 1.225, 1.275, 1.350, 1.400, 1.525 ]),
'gamma_sigma_plus': np.array([ 0.096, .133, .122, .122, .109, .120, .129, .108, .101, .140  ]),
'gamma_sigma_minus': np.array([ 0.096, .133, .122, .122, .109, .120, .129, .108, .101, .140 ]),
}



# 
# data_thermal_history_Lidz_2010 = {
# 'name': "Lidz+ 2010",
# 'z': np.array([ 2.4,  ]),
# 'T0': np.array([ 10000 ]),
# 'T0_sigma_plus': np.array([ 3200 ]),
# 'T0_sigma_minus':np.array([ 2100 ]),
# 'gamma': np.array([ 1.54  ]),
# 'gamma_sigma_plus': np.array([ 0.11 ]),
# 'gamma_sigma_minus':np.array([ 0.11 ]),
# },  