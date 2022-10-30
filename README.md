# EcoHydroSlide-Software
The eco-hydro-geotechnical model is written by languages of Fortran and Matlab. The framework is executed by several steps that we will introduce here. Users are advised to read this manual along with the paper we have published. First of all, we'd like to declare the computational hardware that is necessary for modeling.

#### Hardware

- Operating system: Linux
- Random access memory (RAM): > 8G
- Available Storage: > 50G

We have tested the model on [Delft High Performance Computing Centre](https://www.tudelft.nl/dhpc) (**DHPC**). The operating system version is Red Hat Enterprise Linux release 8.6 (Ootpa). The model is also expected to be executed on other Linux release versions. MATLAB environment (with [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)) is also required for the modeling workflow. 

#### Data

All the data have been prepared in this repository except two folders that including the necessary forcing data. They are not provided in GitHub repository because they could occupy a large storage that exceeds the upload limit. We adopted the secondary storage path to download them: https://drive.google.com/drive/folders/1b6pWfpgdGv6C2NMO_4v7rAvH1WUlN6FS?usp=sharing

If they are rightly found by the users, the folders names will be found as **LandslideBasics** and **DownscalingBasicData**. Users do nothing but only need do put them in the first level directory, i.e., the same path as the folders like **Basics** or **Results**.

#### Run the model step by step

Now, all the preliminary work is done and users can start running the model. You will easily find there are three files named **Bare_100yr_linux, SVeg_100yr_linux, and DVeg_100yr_linux**. 
