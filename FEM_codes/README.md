This folder contains the following Python (Jupyter notebook) code:

'GRWpaper_Rubber_FEM.ipynb' is the Python (Jupyter notebook) code to solve the rubber problem using finite element method and compare the numerical results with experimental data. To run this code: you first  launch Jupyter Notebook  and then you should open the .ipynb file and then run it by  clicking the "Run" button (or pressing Shift + Enter). It will save two .mat files for solution data, the first one is 'matlab_data_st_FEM_rubber1.mat' for the position of moving boundary s(t) and the second one is 'matlab_data_con_FEM_rubber1.mat' for the concentration profile m(x,t) at T = 31 minutes. These two solution .mat files are used to compare the solution with GRW and RW to produce Figure 2 in the manuscript (For more details, see the folder: Rubber). 

'GRWpaper_Stefan_FEM.ipynb' is the Python (Jupyter notebook) code to compute the estimated orders of convergence (EOC) for Stefan problem using finite element method. To run this code: you first  launch Jupyter Notebook  and then you should open the .ipynb file and then run it by  clicking the "Run" button (or pressing Shift + Enter). It will give the orders of convergence (EOC) in space and computing time for Stefan problem using FEM.

The file 'matlab_data_st_FEM_rubber1.mat' contains the  finite element method  (FEM) approximation of the position of moving boundary s(t) for the rubber problem. This data is used to produce FEM solution in Figure 2 of the manuscript. 

The file 'matlab_data_con_FEM_rubber1.mat' contains the  finite element method  (FEM) approximation of concentration profile m(x,t) for the rubber problem. This data is used to produce FEM solution in Figure 2 of the manuscript. 

The file 'matlab_data_st_RW_rubber1.mat' contains the  random walk method (RWM) approximation of the position of moving boundary s(t) for the rubber problem. This data is used to produce the RW solution in Figure 2 of the manuscript. 

The file 'matlab_data_con_RW_rubber1.mat' contains the Random Walk Method (RWM) approximation of the concentration profile m(x,t) for the rubber problem. This data is used to produce the RW solution in Figure 2 of the manuscript."
The file 'FigMB.pdf' contains a comparison of the approximated moving boundary results obtained through FEM, GRM, RW, and experimental data from the lab. This is used in Figure 2.

'The file 'FigCon.pdf' presents a comparison of the approximated concentration profiles obtained through FEM, GRM, and RW methods. This is used in Figure 2.