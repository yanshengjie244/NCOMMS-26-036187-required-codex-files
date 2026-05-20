Anomalous rheology symbolic regression code and source data

This folder contains the Python code, local Python environment, and figure source
data used for the anomalous rheology symbolic-regression analysis.

Environment
-----------
The code was checked with the bundled Windows Python executable:

    python.exe

The package requirements are listed in requirements.txt. The active local
environment was verified with:

    Python 3.8.20
    numpy 1.24.3
    pandas 1.5.2
    matplotlib 3.7.2
    sympy 1.10.1
    scikit-learn 1.0.1
    torch 1.12.1
    tqdm 4.66.5
    jupyterlab 2.2.6

To recreate the environment in a new conda installation, use:

    conda create -n anomalous-rheology python=3.8
    conda activate anomalous-rheology
    conda install --file requirements.txt

The submitted folder already contains a working local environment, so the script
can also be run directly with the bundled python.exe.

Running the symbolic regression script
--------------------------------------
From this folder, run:

    python.exe sr.py

By default, sr.py uses Fig.4c-d_data.csv and runs 500 PhySO epochs. For a quick
installation check, run one epoch:

    set PHYSO_EPOCHS=1
    python.exe sr.py

On PowerShell, use:

    $env:PHYSO_EPOCHS="1"
    .\python.exe sr.py

The script writes:

    rheology_input.png
    rheology.png
    rheology.log

Data files
----------
The figure data files are:

    Fig.4a-b_data.csv
    Fig.4a-b_data.xlsx
    Fig.4c-d_data.csv
    Fig.4e-f_data.csv
    Fig.4g-h_data.csv

Fig.4a-b_data.csv was converted from Fig.4a-b_data.xlsx. The Excel file contains
one worksheet with four paired x/y curve groups labeled 0.5, 1, 2, and 3. The CSV
columns are:

    x_0.5, Curve1_0.5, x_1, Curve1_1, x_2, Curve1_2, x_3, Curve1_3

Notes on checks performed
-------------------------
The requirements in requirements.txt were read and the required packages were
verified in the local environment. The sr.py script was checked with a one-epoch
run. During checking, blank trailing entries in the first curve of
Fig.4c-d_data.csv were filtered before training to avoid passing NaN values into
PhySO. Matplotlib was configured to use a non-interactive backend and not require
a LaTeX installation.
