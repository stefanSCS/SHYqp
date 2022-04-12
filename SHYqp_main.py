'''
--This is the driver script for SHYqpV1
--
'''
import SHYqpV1 as SHYqp
##import numpy as np
from matplotlib import pyplot as plt


###Read mechanical data and other global parameters from text file  
data=SHYqp.readData('mat000File.txt')
### echo data 
SHYqp.prtData(data)  


####---Bezier5YS section    
#######Generate the sequence of points and tangents from directional data 
uaxData=SHYqp.uaxLambda(data)  #;print(uaxData)
#######plot Bezier5YS
savePngProtoModel=True  ##change this to 'False' if saving figures is not desired
SHYqp.protoBez5YS_uaxPlot(uaxData,savePngProtoModel)

###'nPIplaneSections' represents the number of PI-plane sections used to calculate Lambda_max
###(do NOT decrease this number, as Lambda_max is actually a minimum over these sections)
nPIplaneSections=31
lbd,vPatch=SHYqp.protoData(uaxData,nPIplaneSections)
print("sectional shape parameter Lambda_max = ",lbd)
SHYqp.protoBez5YS_Plot(lbd,data['shapeBAX'],vPatch,uaxData,savePngProtoModel)


###---SHYqp section 
####NOTE: 
#####'nEquator' represents the number of locations on the biaxial curve sigma_{xy}=0
#####where convexity constraints are imposed; 
#####This number determines the total number of locations on the unit sphere where convexity constraints are enforced.
#####Increasing 'nEquator' too much has significant impact on the overall run-time.
#####In all tests, 'nEquator=200' seemed more than sufficient. 
####NOTE:
#####'epsilon' represents the lower bound of a convexity constraint;
#####If convexity is not achieved, increase 'epsilon' slightly (rather than 'nEquator') 
if(1):
    ######Select the solver   
    #qpSolver='quadprog'  
    qpSolver='cvxopt'       
    ######Calculate SHYqp parameters     
    if(data['assym']):
        vCoeff,ddMon,nQ,nP=SHYqp.dataFitSHYqp(data,uaxData,lbd,qpSolver,nEquator=200,epsilon=0.01)
    else:
        vCoeff,ddMon,nQ,nP=SHYqp.dataFitSHYqpSymm(data,uaxData,lbd,qpSolver,nEquator=200,epsilon=0.01)
    ######Check convexity     
    cvxCheck=SHYqp.SHYqp_HessGaussCheck(vCoeff,ddMon,nQ,nP)
    ######Calculate overall performance report 
    SHYqp.SHYqp_Predictions(uaxData,vCoeff,ddMon,nQ,nP,qpSolver,cvxCheck)
    ######Calculate  plots 
    savePngSHYqp=True
    SHYqp.SHYqp_uax_Plot(uaxData,vCoeff,ddMon,nQ,nP,savePngSHYqp)
    SHYqp.SHYqp_bax_Plot(vCoeff,ddMon,nQ,nP,data['assym'],uaxData['name'],savePngSHYqp)
    SHYqp.SHYqp_surf_Plot(vCoeff,ddMon,nQ,nP,uaxData['name'],savePngSHYqp)    
###Show all plots 
plt.show()
