# EcoHydroSlide-Software
The eco-hydro-geotechnical model is written in languages of Fortran and Matlab. The framework is executed by several steps that we will introduce here. Users are advised to read this manual along with the paper we have published. First of all, we'd like to declare the computational hardware that is necessary for modeling.

#### Hardware

- Operating system: Linux
- Random access memory (RAM): > 8G
- Available Storage: > 50G

We have tested the model on [Delft High Performance Computing Centre](https://www.tudelft.nl/dhpc) (DHPC). The operating system version is Red Hat Enterprise Linux release 8.6 (Ootpa). The model is also expected to be executed on other Linux release versions. The MATLAB environment (with [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)) is also required for the modeling workflow. 

#### Data

All the data have been prepared in this repository except for two folders that include the necessary forcing data. They are not provided in the GitHub repository because they could occupy a large amount of storage that exceeds the upload limit. We adopted the secondary storage path to download them: https://drive.google.com/drive/folders/1b6pWfpgdGv6C2NMO_4v7rAvH1WUlN6FS?usp=sharing

If they are rightly found by the users, the folder names will be found as **<u>LandslideBasics</u>** and **<u>DownscalingBasicData</u>**. Users do nothing but only need to put them in the first level directory, i.e., the same path as the folders like **<u>Basics</u>** or **<u>Results</u>**.

#### Run the model step by step

All the preliminary work is done, and users can start running the model. You will easily find there are three files named ***Bare_100yr_linux*, *SVeg_100yr_linux*, *and DVeg_100yr_linux***. They correspond to the three computational scenarios: (i) a bare soil catchment ("**Bare**"); (ii) a static fully vegetated catchment ("**S-Veg**"); and (iii) a dynamic fully vegetated catchment ("**D-Veg**"). If you are a macOS user, use the versions in **<u>Exe_macOS</u>** (copy the files to the main path).

1. Now you directly run the the executable file :

```bash
./Bare_100yr_linux  #(for Bare)
./SVeg_100yr_linux  #(for SVeg)
./SVeg_100yr_linux  #(for DVeg)
```

2. Then the model will ask for a project file and you will input:

```bash
Simu_100yr   #(can be found as Simu_100yr.Project in main path)
```

The model starts running and will output the hydrological series and grid spatial variables. The outputs (e.g., the soil moisture) will be saved in **<u>Results</u>** folder. The typical grid file for the soil moisture is named like **OutDT_SM_20220101**, with a suffix of simulation moments.

3. After the executable has finished running, move all the soil moisture grid files to the specific folder. For example, if you are testing the scenario of **Bare**, move the files to **<u>SM_bare</u>**.

4. Start running the slope stability module. You can easily find the *.m files on the main path: ***StabilityModel_main_CertainCe.m***, ***StabilityModel_main_randomCe_NoRoot.m***, and ***StabilityModel_main_randomCe.m***. They refer to the computational scenarios as: (i) CertainCe: the failure depth ($c_e$) is set to 1.5 m and 5 m, denoting the fully rooted and partially rooted failures, respectively; (ii) randomCe_NoRoot: $c_e$ is randomly generated ranging from 1.5 m to 5 m, but the root cohesion is set to zero manually; and (iii) randomCe: $c_e$ is randomly generated ranging from 1.5 m to 5 m with the existing root cohesion. 

   Then you can run the *.m file within the MATLAB environment.

5. Check all the results in the **<u>Results</u>** folder.

#### References and manuals

1. Manual for using CREST hydrological model

   [Hydrometeorology and Remote Sensing Laboratory » CREST (ou.edu)](http://hydro.ou.edu/research/crest/)

2. Manual for using the coupled hydrological-geotechnical model

   https://github.com/GuodingChen/iHydroSlide3D_v1.0/blob/master/iHydroSlide3D_v1.0_manual.html





