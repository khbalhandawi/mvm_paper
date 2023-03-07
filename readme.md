# Welcome to mvm_paper: the code repository for the design space exploration paper

## Quickstart guide

Create a virtual environment for your project

**MacOS/Linux**

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Windows**
```
python -m venv .env
.env\Scripts\activate
pip install -r requirements.txt
```

run the examples normally using python.

## List of python examples

Here is a brief description of what each python script does:

|**File**                                                   |  **Description** |
|-----------------------------------------------------------|------------------|
|[man_defs.py](man_defs.py)|Includes the structure of the margin analysis network (MAN) for the FEA strut example (`get_man_combined` and the analytical strut example `get_man`.|
|[strut_design_combined.py](strut_design_combined.py)|Performs margin value analyis (MVM) on the FEA strut example.|
|[strut_design_manufacturability.py](strut_design_manufacturability.py)|Performs margin value analyis (MVM) on the analytical strut example. |
|[postprocess_DOE.py](postprocess_DOE.py)|Performs a full-factorial DOE on all the design parameters of the strut FEA problem and processes results into dataframes and a parallel coordinates plot (PCP).|
|[postprocess_results.py](postprocess_results.py)|Performs a full-factorial DOE on the *input* design parameters and *decisions* of the strut FEA problem and processes results into dataframes and a parallel coordinates plot (PCP).|
|[postprocess_concepts_min_excess.py](postprocess_concepts_min_excess.py)|Postprocesses the results of arbitrary runs from [strut_design_combined.py](strut_design_combined.py) into dataframes for PCPs.|

## Installing R libraries and dependancies

To reproduce the plots shown in the paper you will need R. With the `renv` packages installed use the following commands in an R console to download and install all the dependancies you need:

```
library(renv)
renv::restore()
```

## List of R scripts

You may then run all the R files to produce the plots. Here is a description of their content

|**File**                                                   |  **Description** |  **Depends on** |
|-----------------------------------------------------------|------------------|-----------------|
|[plot_funcs.R](plot_funcs.R)|Plot theme related functions and definitions. Executed at the beginning of other scripts.|~|
|[plot_mvp.R](plot_mvp.R)|Plots the margin value map (MVP) of various runs and superimposes them.|[strut_design_combined.py](strut_design_combined.py)|
|[plot_scatter_matrix.R](plot_scatter_matrix.R)|Plots a scatter matrix of a particular run and margin node.|[strut_design_combined.py](strut_design_combined.py)|
|[visualize_impact.R](visualize_impact.R)|Plots a bar chart of the impact on performance for the manufacturability example.|[strut_design_manufacturability.py](strut_design_manufacturability.py)|

## For development

First clone this repository and its submodule ``mvmlib``:

```
git clone https://github.com/khbalhandawi/mvm_paper
git submodule init
git submodule update
cd mvm_paper
```

in the ```mvm_paper`` directory do the following:

**MacOS/Linux**

Create a virtual environment to develop this library

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements_dev.txt
pip install -r mvmlib/requirements_dev.txt
```

**Windows**
```
python -m venv .env
.env\Scripts\activate
pip install -r requirements_dev.txt
pip install -r mvmlib/requirements_dev.txt
```

make changes to ``mvm_paper`` and ``mvmlib`` as necessary. Remember to commit changes to both repositories to avoid breaking the API

First for ``mvmlib`` got to its directory and use git commands as usual
```
git checkout -b <your development branch>
git add -A
git commit -m "your message"
(optional) git push origin <your development branch>
```

Then for ``mvm_paper`` simply
```
git add -A
git commit -m "your message"
git push origin master
```