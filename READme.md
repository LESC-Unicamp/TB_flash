
<div style="text-align: center; font-size: 40px; font-weight: bold;">
Tβ-flash code
</div>

*Authors*: **Arley A. Cruz; Estefânia P. Canzian; Luís F. M. Franco**
### The structure of code is divided into:
	
1. **main.py**
2. **EoS.py**
3. **Saturation_point.py**
4. **Tp_flash.py**
5. **Functionalities.py**
6. **Experimental_points.py**
7. **Parameters.py**
8. **Output**
9. **requirements.txt**

### Explanation of each *.py* file:
1. The file *main.py* is the user-executable part of the program. The input parameters necessary for execution of program are:

   - **T**: Temperature (K)
   - **tol**: Algorithm convergence tolerance
   - **n_ite**: Maximum number of iterations
   - **save**: Save the results (1 ---> save; other number ---> don't save)
   - **complete**: Complete binary diagram (1 ---> yes; other number ---> no)
   - **beta_v**: volumetric fraction
   - **case**: Which result to obtain
   - **ftol**: tolerance for pressure optimization
   - **system**: available mixtures dataset

   The output of **main.py**:
   - Phase equilibrium diagram --------------> case = 0
   - Computational efficiency plots ----------> case = 1
   - Table of results ------------------------------> case = 2
   - Deviation (AARD) --------------------------> case = 3
   

2. **EoS.py** is the file which Peng-Robinson EoS remains. The user doesn't need to change this part of program.


3. In the file **Saturation_point.py** the calculation of saturation pressure (bubble and dew points) are accomplished. 
It uses bubble and dew P algorithms for obtaining the range of pressure which is used in the optimization process in **Functionalities.py**.


4. **Tp_flash.py** file is responsible for the equilibrium calculation to a specific pressure and temperature. 
The output of such algorithm is the volumetric fraction in the vapor phase (Beta).


5. **Functionalities.py** file is the heart of the entire code which gives all the results. 
It runs a minimization problem which the objective function is the difference between calculated volumetric fraction (given in **Tp_flash**) 
and specified volumetric fraction (**beta_v**). The bounds of minimization problem are bubble and dew points (obtained in **Saturation_point.py**).
The output value of minimization problem is the pressure of system, in that specific volumetric fraction.
Also, is possible to take phase equilibrium plots for binary and ternary systems and compare with experimental points. 


6. **Experimental_points.py** contains some experimental phase equilibrium points to evaluate the code performance. The systems available are:
   - R32+R125 ---------------------------> system = 0
   - R32+R152 ---------------------------> system = 1
   - Propane+Butane -------------------> system = 2
   - CO2+Butane ------------------------> system = 3
   - Ethylene+CO2 ----------------------> system = 4
   - Methane+CO2 ----------------------> system = 5
   - Methane+Ethane+CO2 -----------> system = 6
   
   OBS: The system 6 doesn't support cases 1 and 2
   

7. **Parameters.py** contains all necessary input parameters of each component for running the code.


8. **Output.py** is responsible for organizing all possible executions (cases and systems).


9. **requirements.txt** is the list of all python extensions needed for running the code.To install such extensions use the command:
    **pip install -r requirements.txt** 