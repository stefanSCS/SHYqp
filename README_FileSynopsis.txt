
-------Materials in article Figures
Correspondence between the SHYqp models shown in the Figures of the article 
and the material input files used to generate them. 
One can then trace the corresponding *Err_and_Coeff.txt file and generated figs 
by the 'name' field in the material input file. 

Note: This mapping applies only to the revised version of the article (not the preprint).

Fig.6-7: matAZ31B_Lou2007.txt
Fig.8-9: matAZ31B_Lou2007_B.txt
Fig.10-11: matAZ31B_Andar2012.txt
Fig.12-13: matTiG4_Raemy2017B.txt
Fig.14-15: matAA5042H2.txt
Fig.16-17: matAA2090T3.txt
Fig.18-19: matDP980_Li2020.txt
Fig.22-23: [matAZ31B_Lou2007_TXT.txt (Model-A), 
            matAZ31B_Lou2007_TXT2.txt (Model-B),
            matAZ31B_Lou2007_TXT2B.txt (Model-C)]
Fig.24-25: [matAZ31B_Andar2012_00.txt (Model-A),
            matAZ31B_Andar2012_002.txt (Model-B),
            matAZ31B_Andar2012_002B.txt (Model-C)]


Input files (with biaxial data):
matAZ31B_Lou2007_TXT2B.txt (Model-C): AZ31B_biaxYS_Data.txt 
matAZ31B_Andar2012_002B.txt (Model-C): AZ31B_Biaxial_Andar.txt


-------AZ31 Damask input/output data

Texture: EulerData_Sim555b_2607.txt
Note: 
The angles were generated using classical (Z-X-Z) Euler angles formulas. 
Damask requires the first and third angles to be in [0,2*pi]. 
When processing, add 2*pi to the negative angles only.

Geometry: AZ31B_2607_16x16x16.vti

Calculated biaxial stresses: AZ31B_biaxYS_Data.txt 
