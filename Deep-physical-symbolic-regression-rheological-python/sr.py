# External packages
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import torch
import physo
import physo.learn.monitoring as monitoring
plt.rcParams["text.usetex"] = False
import pandas as pd
import chardet
seed = 2
np.random.seed(seed)
torch.manual_seed(seed)
#Making toy synthetic data

data = np.genfromtxt('Fig.4c-d_data.csv', delimiter=',', skip_header=1)
data = data[np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])]
inputdata1= data[:, 0]
#inputdata2= data[:, 1]
#inputdata3= data[:, 2]
outputdata=data[:, 1]
X1 = inputdata1.reshape(-1,1)
#X2 = inputdata2.reshape(-1,1)
#X3= inputdata3.reshape(-1,1)
X = np.hstack([X1])
X = X.T
y = outputdata
print(f"X1 shape: {X1.shape}")
#print(f"X2 shape: {X2.shape}")
#print(f"X2 shape: {X3.shape}")
print(f"y shape: {y.shape}")
plt.figure()
plt.scatter(X1, y, color='k', marker='.')
  # Set x-axis to log scale
#plt.yscale('log')  # Set y-axis to log scale plt.xlabel("X (log scale)")
#plt.ylabel("y (log scale)")
plt.title("rheology")
plt.tight_layout()
plt.savefig("rheology_input.png", dpi=300)
plt.close()

save_path_training_curves = 'rheology.png'
save_path_log             = 'rheology.log'

run_logger     = lambda : monitoring.RunLogger(save_path = save_path_log,
                                                do_save = True)

run_visualiser = lambda : monitoring.RunVisualiser (epoch_refresh_rate = 1,
                                           save_path = save_path_training_curves,
                                           do_show   = False,
                                           do_prints = True,
                                           do_save   = True, )
# Running SR task
expression, logs = physo.SR(X,y,
                            # Giving names of variables (for display purposes)
                            X_names = [ "x1" ],
                            # Associated physical units (ignore or pass zeroes if irrelevant)
                            X_units = [ [1]],
                            # Giving name of root variable (for display purposes)
                            y_name  = "y",
                            y_units = [1],
                            # Fixed constants
                            fixed_consts       = [  ],
                            fixed_consts_units = [  ],
                            # Free constants names (for display purposes)
                            free_consts_names = [ "a" ,"b","c","d","e"],
                            free_consts_units = [ [1],[1],[1],[1],[1]],
                            # Symbolic operations that can be used to make f
                            op_names = ["mul", "add", "sub", "div" ,"pow","exp","log"],
                            get_run_logger     = run_logger,
                            get_run_visualiser = run_visualiser,
                            # Run config
                            run_config = physo.config.config2.config2,
                            # Parallel mode (only available when running from python scripts, not notebooks)
                            parallel_mode = False,
                            # Number of iterations
                            epochs = int(os.environ.get("PHYSO_EPOCHS", "500"))
)
best_expr = expression
print(best_expr.get_infix_pretty())




