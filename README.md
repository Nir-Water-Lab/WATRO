# WATRO

## Weak Acids Transport in Reverse Osmosis

<img src= https://github.com/Nir-Water-Lab/WATRO_1/blob/main/watro_1.2_install/WATRO_1.2_Install/Watro%20pic.PNG>

### WATRO: Installation Instruction for Windows
WATRO (Weak Acids Transport Reverse in Osmosis) 2021 release. Including user 
interface (not graphic) for SWRO 1st and 2nd pass (experimentally validated) and 
brackish water (not experimentally validated yet), Now works with installation of 
the Anaconda Python 3 package.

1. Download and install the Anaconda python package (works only with Python 
2, not 3): https://www.anaconda.com/download.
2. Download and install IPhreeqc COM modules (msi files at the bottom of the 
page): https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.
3. Extract the content of the Installation folder to a chosen location in your hard 
drive.
4. Go to the IPhreeqc folder. Find the file ‘IPhreeqcCOM.dll’, copy it and paste in 
the same folder where you placed the WATRO scripts.
5. Go to the IPhreeqc folder. From the database folder in there copy all the ‘ .dat’ 
files and paste them in the same folder with the WATRO code.
6. Installation of .Net Framework 3.5 from Microsoft may be needed.
7. You are ready to go. Open one of the interfaces using Anaconda, follow the 
instructions and run the code.

The WATRO (Weak Acid Transport Reverse Osmosis) computer code was Initially 
developed by Oded Nir and Ori Lahav at the Technion - Israel Institute of 
Technology, Faculty of Civil and Environmental Engineering. Further development 
is now performed by Oded Nir in the Zuckerberg Institute for Water Research, 
Blaustein Institutes for Desert Research, Ben Gurion University, Israel.


## References
1. **O. Nir**, O. Lahav, 2016, Acid-base dynamics in seawater reverse osmosis: experimental evaluation 
of a reactive-transport algorithm, Environmental Science: Water Research and Technology 2(1), 
107-116.
2. **O. Nir**, When does commercial software fail in predicting scaling tendency in reverse osmosis and 
what can we do better? Desalination and Water Treatment. 2018 Nov 1;131:34-42
