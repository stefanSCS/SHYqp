'''version: 'SHYqpV1' (initial)
--This code generates the Bezier5YS and SHYqp-models of a sheet metal based on uniaxial and biaxial mechanical test data;
--Materials can be symmetric/assymetric w.r.t. tension-compression;
--Read the accompanying 'mat001File_Instructions.txt' file for input instructions; 
--Theory in the paper: 
'Bezier5YS and SHYqp: 
A general framework for generating data and for modeling symmetric and asymmetric orthotropic yield surfaces';
--This proof of concept code is released under the MIT licence by its author Stefan C. Soare.
'''

import numpy as np
import quadprog as qpg
import cvxopt
from matplotlib import pyplot as plt
from os import path as osp

figDir='.\\FIGS\\'
lineWidthUax=1.5

def readData(fName):
    try:
        ff=open(fName,'r')
    except IOError as err:
        print(err)
        exit()        
    vLine=[];fVal=0.0
    data={'name':'','assym':'*','DEG':'*',
            'sT0':'*','sT15':'*','sT22.5':'*','sT30':'*','sT45':'*','sT60':'*','sT67.5':'*','sT75':'*','sT90':'*',
            'rT0':'*','rT15':'*','rT22.5':'*','rT30':'*','rT45':'*','rT60':'*','rT67.5':'*','rT75':'*','rT90':'*',
            'sC0':'*','sC15':'*','sC22.5':'*','sC30':'*','sC45':'*','sC60':'*','sC67.5':'*','sC75':'*','sC90':'*',
            'rC0':'*','rC15':'*','rC22.5':'*','rC30':'*','rC45':'*','rC60':'*','rC67.5':'*','rC75':'*','rC90':'*',
            'sTb':'*','sCb':'*','rTb':'*','rCb':'*',
            'wsTb':'*','wsCb':'*','wrTb':'*','wrCb':'*',
            'ww':'*','LTAN':'*','LUAX':'*','LBAX':'*','VLBAX':np.ones(6),'fileData':'*'}
    vCheck={'name':False,'assym':False,'DEG':False,
            'sT0':False,'sT15':False,'sT225':False,'sT30':False,'sT45':False,'sT60':False,'sT675':False,'sT75':False,'sT90':False,
            'rT0':False,'rT15':False,'rT225':False,'rT30':False,'rT45':False,'rT60':False,'rT675':False,'rT75':False,'rT90':False,
            'sC0':False,'sC15':False,'sC225':False,'sC30':False,'sC45':False,'sC60':False,'sC675':False,'sC75':False,'sC90':False,
            'rC0':False,'rC15':False,'rC225':False,'rC30':False,'rC45':False,'rC60':False,'rC675':False,'rC75':False,'rC90':False,
            'sTb':False,'sCb':False,'rTb':False,'rCb':False,
            'ww':False,'LTAN':False,'LUAX':False,'LBAX':False,'VLBAX':True,'fileData':False}   
    for line in ff:
        errLine=line
        line=line.strip()
        if(line=='' or line[0]=='#'):
            continue
        if('=' not in line):
            print('Incorrect data format (equal sign missing) on line:\n'+errLine+'\nCalculations aborted')
            exit()
        vLine=line.split('=')
        vLine[0]=vLine[0].strip();vLine[1]=vLine[1].strip()
        if(vLine[0]=='' or vLine[1]=='' or len(vLine)!=2):
            print('incorrect(A) data format on line:\n'+errLine+'\nCalculations aborted')
            exit()
        if(vLine[0]=='name'):
            if(len(vLine[0])>36):
                print('Warning: name too long (only the first 30 chars are retained):\n')
            data['name']=vLine[1][0:30]; vCheck['name']=True;continue  
        if(vLine[0]=='assym'):
            if((vLine[1][0]).upper()=='Y' or vLine[1]=='*'):data['assym']=True
            elif(vLine[1][0].upper()=='N'):data['assym']=False
            else:print('assym: unknown value\nCalculations aborted');exit()
            vCheck['assym']=True;continue               
        if(vLine[0] not in ['name','assym','fileData']):            
            try:
                fVal=float(vLine[1])
                if(fVal<0.0):
                    print('incorrect zero-value on line:\n'+errLine+'\nCalculations aborted')
                    exit()
                if(fVal>10.0**6):  
                    print('unacceptable large value on line:\n'+errLine+'\nCalculations aborted')
                    exit()                
            except ValueError:
                if(vLine[1]=='*'):
                    fVal='*'
                else:    
                    print('incorrect(B) data format on line:\n'+errLine+'\nCalculations aborted')
                    exit()
        if(vLine[0]=='DEG'):
            if(fVal=='*'):data['DEG']=4
            else:
                if(fVal-int(fVal)):print("DEG: must be an integer");exit()
                if(int(fVal) not in [2*k for k in range(2,13)]):
                    print("DEG: incorrect value (must be an even integer >=4 and <=24)");exit()
                data['DEG']=int(fVal)
            vCheck['DEG']=True;continue                
        if(vLine[0]=='sT0'):
            data['sT0']=fVal; vCheck['sT0']=True;continue
        if(vLine[0]=='rT0'):
            data['rT0']=fVal; vCheck['rT0']=True;continue
        if(vLine[0]=='sT90'):
            data['sT90']=fVal; vCheck['sT90']=True;continue 
        if(vLine[0]=='rT90'):
            data['rT90']=fVal; vCheck['rT90']=True;continue
        if(vLine[0]=='sTb'):
            data['sTb']=fVal; vCheck['sTb']=True;continue 
        if(vLine[0]=='rTb'):
            data['rTb']=fVal; vCheck['rTb']=True;continue            
        if(vLine[0]=='sT15'):
            data['sT15']=fVal; vCheck['sT15']=True;continue
        if(vLine[0]=='sT225'):
            data['sT22.5']=fVal; vCheck['sT225']=True;continue            
        if(vLine[0]=='sT30'):
            data['sT30']=fVal; vCheck['sT30']=True;continue
        if(vLine[0]=='sT45'):
            data['sT45']=fVal; vCheck['sT45']=True;continue
        if(vLine[0]=='sT60'):
            data['sT60']=fVal; vCheck['sT60']=True;continue
        if(vLine[0]=='sT675'):
            data['sT67.5']=fVal; vCheck['sT675']=True;continue    
        if(vLine[0]=='sT75'):
            data['sT75']=fVal; vCheck['sT75']=True;continue
        if(vLine[0]=='rT15'):
            data['rT15']=fVal; vCheck['rT15']=True;continue
        if(vLine[0]=='rT225'):
            data['rT22.5']=fVal; vCheck['rT225']=True;continue            
        if(vLine[0]=='rT30'):
            data['rT30']=fVal; vCheck['rT30']=True;continue
        if(vLine[0]=='rT45'):
            data['rT45']=fVal; vCheck['rT45']=True;continue
        if(vLine[0]=='rT60'):
            data['rT60']=fVal; vCheck['rT60']=True;continue
        if(vLine[0]=='rT675'):
            data['rT67.5']=fVal; vCheck['rT675']=True;continue    
        if(vLine[0]=='rT75'):
            data['rT75']=fVal; vCheck['rT75']=True;continue 
        if(vLine[0]=='sC0'):
            data['sC0']=fVal; vCheck['sC0']=True;continue
        if(vLine[0]=='rC0'):
            data['rC0']=fVal; vCheck['rC0']=True;continue
        if(vLine[0]=='sC90'):
            data['sC90']=fVal; vCheck['sC90']=True;continue 
        if(vLine[0]=='rC90'):
            data['rC90']=fVal; vCheck['rC90']=True;continue
        if(vLine[0]=='sCb'):
            data['sCb']=fVal; vCheck['sCb']=True;continue 
        if(vLine[0]=='rCb'):
            data['rCb']=fVal; vCheck['rCb']=True;continue             
        if(vLine[0]=='sC15'):
            data['sC15']=fVal; vCheck['sC15']=True;continue   
        if(vLine[0]=='sC30'):
            data['sC30']=fVal; vCheck['sC30']=True;continue
        if(vLine[0]=='sC45'):
            data['sC45']=fVal; vCheck['sC45']=True;continue
        if(vLine[0]=='sC60'):
            data['sC60']=fVal; vCheck['sC60']=True;continue
        if(vLine[0]=='sC75'):
            data['sC75']=fVal; vCheck['sC75']=True;continue
        if(vLine[0]=='rC15'):
            data['rC15']=fVal; vCheck['rC15']=True;continue   
        if(vLine[0]=='rC30'):
            data['rC30']=fVal; vCheck['rC30']=True;continue
        if(vLine[0]=='rC45'):
            data['rC45']=fVal; vCheck['rC45']=True;continue
        if(vLine[0]=='rC60'):
            data['rC60']=fVal; vCheck['rC60']=True;continue
        if(vLine[0]=='rC75'):
            data['rC75']=fVal; vCheck['rC75']=True;continue
        if(vLine[0]=='sC225'):
            data['sC22.5']=fVal; vCheck['sC225']=True;continue
        if(vLine[0]=='sC675'):
            data['sC67.5']=fVal; vCheck['sC675']=True;continue
        if(vLine[0]=='rC225'):
            data['rC22.5']=fVal; vCheck['rC225']=True;continue    
        if(vLine[0]=='rC675'):
            data['rC67.5']=fVal; vCheck['rC675']=True;continue
        if(vLine[0]=='ww'):
            data['ww']=fVal
            if(fVal=='*' or fVal>1.0):
                #data['ww']=0.7
                data['ww']=0.9
            vCheck['ww']=True;continue
        if(vLine[0]=='LTAN'):
            data['LTAN']=fVal
            if(fVal=='*' or fVal>1.0):
                data['LTAN']=0.5  
            vCheck['LTAN']=True;continue            
        if(vLine[0]=='LUAX'):
            data['LUAX']=fVal
            if(fVal=='*' or fVal>1.0):
                data['LUAX']=0.6
            vCheck['LUAX']=True;continue
        if(vLine[0]=='LBAX'):
            if(fVal=='*' or fVal>1.0):
                data['LBAX']=0.0  ##one shape parameter
            else:    
                data['LBAX']=fVal    
            vCheck['LBAX']=True;continue
        if(vLine[0]=='LBAX1'):
            if(fVal=='*' or fVal>1.0):
                data['VLBAX'][0]=1.0
            else:
                data['VLBAX'][0]=fVal 
            continue                
        if(vLine[0]=='LBAX2'):
            if(fVal=='*' or fVal>1.0):
                data['VLBAX'][1]=1.0
            else:
                data['VLBAX'][1]=fVal
            continue                
        if(vLine[0]=='LBAX3'):
            if(fVal=='*' or fVal>1.0):
                data['VLBAX'][2]=1.0
            else:
                data['VLBAX'][2]=fVal
            continue                
        if(vLine[0]=='LBAX4'):
            if(fVal=='*' or fVal>1.0):
                data['VLBAX'][3]=1.0
            else:
                data['VLBAX'][3]=fVal
            continue                
        if(vLine[0]=='fileData'):
            if(vLine[1] !='*'):
                try:
                    ffd=open(vLine[1],'r')
                    ysz=[]
                    for line in ffd:
                        line=line.strip().split(',') ##;print(line,len(line))
                        if(len(line)==2):
                            ysz.append([float(line[0]),float(line[1]),0.0])
                        elif(len(line)==3):
                            ysz.append([float(line[0]),float(line[1]),float(line[2])])                            
                    ffd.close()
                    data['fileData']=np.array(ysz)
                    #ysz=np.array(ysz)
                    #idx=np.argwhere(ysz[:,0]>0)
                    #ysRD=ysz[idx][np.argmin(np.abs(ysz[idx,1]))][0][0]
                    #data['fileData']=ysz/ysRD ##normalize by stress along RD
                    #if(0):##check input file 
                    #    ffd=open(vLine[1].split('.')[0]+'_out.txt','w')
                    #    for ys in data['fileData']:
                    #        ffd.write('{:.5f}, {:.5f}\n'.format(ys[0],ys[1]))
                    #    ffd.close()                        
                except IOError as err:
                    print(err)
                    exit()
            else:
                data['fileData']=np.array([])
            vCheck['fileData']=True;continue                
    if('*' in [data['sT0'],data['sT45'],data['sT90'],data['sC0'],data['sC45'],data['sC90']]):
        print('Missing value: sT0,sT45,sT90,sC0,sC45,sC90 must be provided\nCalculations aborted');exit()
    if('*' in [data['rT0'],data['rT45'],data['rT90'],data['rC0'],data['rC45'],data['rC90']]):
        print('Missing value: rT0,rT45,rT90,rC0,rC45,rC90 must be provided\nCalculations aborted');exit()  
    for item in vCheck:
        if(not vCheck[item]):
            print(item+': no data provided\nCalculations aborted')
            exit()     
    ff.close()
    s0=data['sT0']
    for key in data:
        if((key[0]=='s') and (data[key]!='*')):data[key]/=s0
    if(data['sTb']=='*'):
        data['wsTb']=0
        data['sTb']=0.5*(data['sT0']+data['sT90'])
    else:
        data['wsTb']=1    
    if(data['sCb']=='*'):
        data['wsCb']=0
        data['sCb']=0.5*(data['sC0']+data['sC90'])
    else:
        data['wsCb']=1    
    if(data['rTb']=='*'):data['rTb']=1.0;data['wrTb']=0
    else:data['wrTb']=1
    if(data['rCb']=='*'):data['rCb']=1.0;data['wrCb']=0
    else:data['wrCb']=1
    data['fileData']/=s0
    npi=np.pi
    dTheta={'0':0,'15':npi/12,'22.5':npi/8,'30':npi/6,'45':npi/4,'60':npi/3,'67.5':3*npi/8,'75':5*npi/12,'90':npi/2}
    options=[['0','45','90'],['0','22.5','45','67.5','90'],['0','15','30','45','60','75','90']]
    if(data['assym']):
        aData={'name':data['name'],'assym':True,'DEG':data['DEG'],
               'sT':{},'sC':{},'rT':{},'rC':{},'thetaT':{},'thetaC':{},
               'sTb':data['sTb'],'sCb':data['sCb'],'rTb':data['rTb'],'rCb':data['rCb'],
               'wsTb':data['wsTb'],'wsCb':data['wsCb'],'wrTb':data['wrTb'],'wrCb':data['wrCb'],
               'weight':data['ww'],'shapeTAN':data['LTAN'],'shapeUAX':data['LUAX'],'shapeBAX':data['LBAX'],
               #'shapeVBAX':data['VLBAX']}
                'shapeVBAX':np.array([data['VLBAX'][0],data['VLBAX'][1],data['VLBAX'][0],data['VLBAX'][2],data['VLBAX'][3],data['VLBAX'][2]]),
                'fileData':data['fileData']}                 
        for key in data:
            if(key in ['VLBAX','fileData'] or data[key]=='*'):continue
            if(key[0:2]=='sT' and key!='sTb'):aData['sT'][key[2:]]=data[key];aData['thetaT'][key[2:]]=dTheta[key[2:]]
            if(key[0:2]=='rT' and key!='rTb'):aData['rT'][key[2:]]=data[key];aData['thetaT'][key[2:]]=dTheta[key[2:]]
            if(key[0:2]=='sC' and key!='sCb'):aData['sC'][key[2:]]=data[key];aData['thetaC'][key[2:]]=dTheta[key[2:]]
            if(key[0:2]=='rC' and key!='rCb'):aData['rC'][key[2:]]=data[key];aData['thetaC'][key[2:]]=dTheta[key[2:]]
            if(key=='LBAX'):
                if(int(data[key])==0):
                    data['VLBAX']=data['VLBAX'][0]*np.ones(6) ##all shape parameters equal to 'LBAX1'
        if(len(aData['sT'])!=len(aData['thetaT']) or len(aData['sT'])!=len(aData['rT'])):
            print('sT and rT: the numbers of directional angles, stresses and r-values are not the same\nCalculations aborted');exit()
        if(len(aData['sC'])!=len(aData['thetaC']) or len(aData['sC'])!=len(aData['rC'])):
            print('sC and rC: the numbers of directional angles, stresses and r-values are not the same\nCalculations aborted');exit()
        for key in ['sT','rT','sC','rC','thetaT','thetaC']:    
            if([k for k in aData[key]] not in options):
                print("Unknown combination of angles. The options allowed are:")
                print(options[0]);print(options[1]);print(options[2]);exit()                
    else:
        aData={'name':data['name'],'assym':False,'DEG':data['DEG'],
               'sT':{},'sC':{},'rT':{},'rC':{},'thetaT':{},'thetaC':{},
               'sTb':data['sTb'],'sCb':data['sTb'],'rTb':data['rTb'],'rCb':data['rTb'],
               'wsTb':data['wsTb'],'wsCb':data['wsTb'],'wrTb':data['wrTb'],'wrCb':data['wrTb'],
               'weight':data['ww'],'shapeTAN':data['LTAN'],'shapeUAX':data['LUAX'],'shapeBAX':data['LBAX'],
               #'shapeVBAX':np.array([data['VLBAX'][0],data['VLBAX'][1],data['VLBAX'][2],data['VLBAX'][0],data['VLBAX'][1],data['VLBAX'][2]])}
               'shapeVBAX':np.array([data['VLBAX'][0],data['VLBAX'][1],data['VLBAX'][0],data['VLBAX'][0],data['VLBAX'][1],data['VLBAX'][0]]),
               'fileData':data['fileData']}                 
        for key in data:
            if(key in ['VLBAX','fileData'] or data[key]=='*'):continue
            if(key[0:2]=='sT' and key!='sTb'):aData['sT'][key[2:]]=data[key];aData['thetaT'][key[2:]]=dTheta[key[2:]]
            if(key[0:2]=='rT' and key!='rTb'):aData['rT'][key[2:]]=data[key];aData['thetaT'][key[2:]]=dTheta[key[2:]]
            if(key=='LBAX'):
                if(int(data[key])==0):
                    data['VLBAX']=data['VLBAX'][0]*np.ones(6) ##all shape parameters equal to 'LBAX1'
        if(len(aData['sT'])!=len(aData['thetaT']) or len(aData['sT'])!=len(aData['rT'])):
            print('sT and rT: the numbers of directional angles, stresses and r-values are not the same\nCalculations aborted');exit()
        for key in ['sT','rT','thetaT']:    
            if([k for k in aData[key]] not in options):
                print("Unknown combination of angles. The options allowed are:")
                print(options[0]);print(options[1]);print(options[2]);exit() 
        aData['sC']=aData['sT'];aData['rC']=aData['rT'];aData['thetaC']=aData['thetaT']                
    return aData 


'''
function:'prtData'
--Utility for displaying the input data 
'''
def prtData(data):
    for key in data:
        if(type(data[key]) in [str,bool]):
            print(key+':{}'.format(data[key]));continue
        if(type(data[key])==float):
            print(key+':{:.3f}'.format(data[key]));continue
        if(type(data[key])==int):
            print(key+':{}'.format(data[key]));continue
        if(type(data[key])==dict):
            txt=key+':{'
            for kk in data[key]:
                txt+="'"+kk+"'"+':{:.3f}'.format(data[key][kk])+','
            print(txt[:len(txt)-1]+'}');continue        
        if(type(data[key])==list):
            print((key+':['+'{:.3f},'*(len(data[key])-1)+'{:.3f}'+']').format(*data[key]));continue
    return                


def nMonomials():
    vDeg=[]
    for kk in range(2,13):
        deg=2*kk
        jj=0;sOdd=0;sEven=0
        while(jj<deg):
            sOdd+=deg-jj
            sEven+=deg+1-jj
            jj+=2
        sEven+=1
        vDeg.append((deg,sOdd,sEven,sOdd+sEven))
    return vDeg 
def nMonomials2():
    vDeg=[]
    for kk in range(2,13):
        deg=2*kk
        sEven=(kk+1)**2
        sOdd=kk*(kk+1)
        vDeg.append((deg,sOdd,sEven,sOdd+sEven))
    return vDeg    


def nMonoms(degQ):
    deg=degQ//2    
    return (deg+1)**2,deg*(deg+1)


def vPoly(degree):
    nQ=int(degree)
    if(nQ%2):print("'degree' must be an even integer\nCalculations aborted");exit()
    dd={'nQ':nQ}
    vP=[];vPidx=[];jj=0
    while(jj<=nQ-1):
        lv=len(vP)
        vPidx.append([jj,[k for k in range(lv,lv+nQ-jj)]])
        vP+=[(nQ-1-jj-k,k,jj) for k in range(nQ-jj)]
        jj+=2
    dd['vP']=vP ;dd['vPidx']=vPidx   
    dd['vD1P']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in vP]
    dd['vC1P']=np.array([mon[0] if(mon[0]>0) else 0 for mon in vP])
    dd['vD2P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in vP]
    dd['vC2P']=np.array([mon[1] if(mon[1]>0) else 0 for mon in vP])
    dd['vD3P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in vP]
    dd['vC3P']=np.array([mon[2] if(mon[2]>0) else 0 for mon in vP]) 
    dd['vH11P']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH11P']=np.array([cf*mon[0] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH12P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH12P']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH13P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD1P']]
    dd['vCH13P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC1P'],dd['vD1P'])])
    dd['vH22P']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD2P']]
    dd['vCH22P']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC2P'],dd['vD2P'])])
    dd['vH23P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD2P']]
    dd['vCH23P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC2P'],dd['vD2P'])])
    dd['vH33P']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD3P']]
    dd['vCH33P']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC3P'],dd['vD3P'])])    
    vQ=[];vQidx=[];jj=0
    while(jj<=nQ):
        lv=len(vQ)
        vQidx.append([jj,[k for k in range(lv,lv+nQ+1-jj)]])
        vQ+=[(nQ-jj-k,k,jj) for k in range(nQ+1-jj)]
        jj+=2 
    dd['vQ']=vQ;dd['vQidx']=vQidx
    dd['vD1Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in vQ]
    dd['vC1Q']=np.array([mon[0] if(mon[0]>0) else 0 for mon in vQ])
    dd['vD2Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in vQ]
    dd['vC2Q']=np.array([mon[1] if(mon[1]>0) else 0 for mon in vQ])
    dd['vD3Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in vQ]
    dd['vC3Q']=np.array([mon[2] if(mon[2]>0) else 0 for mon in vQ])
    dd['vH11Q']=[(mon[0]-1,mon[1],mon[2]) if(mon[0]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH11Q']=np.array([cf*mon[0] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH12Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH12Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH13Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD1Q']]
    dd['vCH13Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC1Q'],dd['vD1Q'])])
    dd['vH22Q']=[(mon[0],mon[1]-1,mon[2]) if(mon[1]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH22Q']=np.array([cf*mon[1] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH23Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD2Q']]
    dd['vCH23Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC2Q'],dd['vD2Q'])])
    dd['vH33Q']=[(mon[0],mon[1],mon[2]-1) if(mon[2]>0) else (0,0,0) for mon in dd['vD3Q']]
    dd['vCH33Q']=np.array([cf*mon[2] for (cf,mon) in zip(dd['vC3Q'],dd['vD3Q'])])     
    return dd


'''
function: 'prtMonomials'
--Utility for checking that the list of monomials (and derivatives) is correct 
'''
def prtMonomials(dd):
    print("vP = \n",dd['vP'])
    print('derivatives')
    print("vD1P=\n",dd['vD1P']);print(dd['vC1P']);print("vD2P =\n",dd['vD2P']);print(dd['vC2P']);print("vD3P =\n",dd['vD3P']);print(dd['vC3P'])
    print('Hessian')
    print("vH11P =\n",dd['vH11P']);print(dd['vCH11P']);print("vH12P =\n",dd['vH12P']);print(dd['vCH12P']);print("vH13P =\n",dd['vH13P']);print(dd['vCH13P'])
    print("vH22P =\n",dd['vH22P']);print(dd['vCH22P']);print("vH23P =\n",dd['vH23P']);print(dd['vCH23P']);print("vH33P =\n",dd['vH33P']);print(dd['vCH33P'])
    print('nCoeff P = ',len(dd['vP']))
    print("vQ = \n",dd['vQ'])
    print('derivatives')
    print("vD1Q = \n",dd['vD1Q']);print(dd['vC1Q']);print("vD2Q = \n",dd['vD2Q']);print(dd['vC2Q']);print("vD3Q = \n",dd['vD3Q']);print(dd['vC3Q'])
    print('Hessian')
    print("vH11Q = \n",dd['vH11Q']);print(dd['vCH11Q']);print("vH12Q = \n",dd['vH12Q']);print(dd['vCH12Q']);print("vH13Q = \n",dd['vH13Q']);print(dd['vCH13Q'])
    print("vH22Q = \n",dd['vH22Q']);print(dd['vCH22Q']);print("vH23Q = \n",dd['vH23Q']);print(dd['vCH23Q']);print("vH33Q = \n",dd['vH33Q']);print(dd['vCH33Q'])
    print('nCoeff Q = ',len(dd['vQ']))
    print('total number of coeffs = ',len(dd['vP'])+len(dd['vQ']))
    return

'''
def lambdaMax(bOne,tOne,bTwo,tTwo):
    aa=bTwo-bOne
    ddet=tOne[0,0]*tTwo[1,0]-tOne[1,0]*tTwo[0,0]
    mu=(tTwo[1,0]*aa[0,0]-tTwo[0,0]*aa[1,0])/ddet
    lbd=mu
    mu=(tOne[0,0]*aa[1,0]-tOne[1,0]*aa[0,0])/ddet
    return 0.5*min(lbd,mu)
'''
    
def curveSeg(Bs,Ts,Be,Te,Lshape,plot=True,nTT=5):
    '''
    Bs=start point; Ts=start tangent;Be=end point; Te=end tangent
    All must be numpy vectors of shape (2,1)
    '''
    if(plot):
        nTT=101;vTT=np.linspace(0.0,1.0,nTT)
    else:    
        ##nTT=5;
        vTT=np.linspace(0.3,0.7,nTT)
    TTs=Lshape*Ts/np.sqrt(Ts[0,0]**2+Ts[1,0]**2) 
    TTe=Lshape*Te/np.sqrt(Te[0,0]**2+Te[1,0]**2)
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be
    return Bs+vTT*(A1+vTT*(A2+vTT*(A3+vTT*(A4+vTT*A5))))

def curveSeg3(Bs,Ts,Be,Te,LshapeS,LshapeE,plot=True,nTT=5):
    '''
    The same function as 'curveSeg' with one exception:
    No shape restriction is imposed on the input vectors.
    Note: the caller must ensure that tangents are normalized (as unit vectors)    
    '''
    if(plot):
        nTT=101;vTT=np.linspace(0.0,1.0,nTT)
    else:    
        vTT=np.linspace(0.2,0.8,nTT)
    TTs=LshapeS*Ts;TTe=LshapeE*Te
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be
    return Bs+vTT*(A1+vTT*(A2+vTT*(A3+vTT*(A4+vTT*A5))))
    
    
def curveSegF(Bs,Ts,Be,Te,Lshape,vTT):
    '''
    Mainly the same function as 'curveSeg' with two exceptios:
    -No shape restriction is imposed on the input vectors.
    -Calculates interpolated values for a given vector of t-values (vTT)
    Note: the caller must ensure that tangents are normalized (as unit vectors)    
    '''
    TTs=Lshape*Ts;TTe=Lshape*Te
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be   
    return Bs+vTT*(A1+vTT*(A2+vTT*(A3+vTT*(A4+vTT*A5)))) 

def curveSegFDF(Bs,Ts,Be,Te,Lshape,vTT):
    '''
    Mainly the same function as 'curveSeg' with two exceptios:
    -No shape restriction is imposed on the input vectors.
    -Calculates interpolated values and derivatives for a given vector of t-values (vTT)
    Note: the caller must ensure that tangents are normalized (as unit vectors)    
    '''
    TTs=Lshape*Ts;TTe=Lshape*Te
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be
    f=Bs+vTT*(A1+vTT*(A2+vTT*(A3+vTT*(A4+vTT*A5))))
    df=A1+vTT*(2.0*A2+vTT*(3.0*A3+vTT*(4.0*A4+5.0*A5*vTT)))    
    return f,df    

    
'''
function:'curveSegTheta'
--Calculates the list of t-parameters corresponding to a list of theta values 
  for a Bezier segment defined by Bs,Ts,Be,Te,Lshape 
'''
def curveSegTheta(vTheta,Bs,Ts,Be,Te,Lshape):
    TTs=Lshape*Ts;TTe=Lshape*Te
    B1=Bs+TTs;B2=B1+TTs;B4=Be-TTe;B3=B4-TTe
    A1=5.0*(B1-Bs)
    A2=10.0*(Bs+B2-2.0*B1)
    A3=10.0*(3.0*B1-Bs+B3-3.0*B2)
    A4=5.0*(Bs-4.0*B1+6.0*B2-4.0*B3+B4)
    A5=5.0*(B1-2.0*B2+2.0*B3-B4)-Bs+Be
    vtt=[]
    t=0.0;epsF=1.0e-9;epsDF=1.0e-12
    for theta in vTheta:
        nn=0
        while(nn<500):
            f=Bs+t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))-theta
            if(abs(f)<=epsF):break            
            df=A1+t*(2.0*A2+t*(3.0*A3+t*(4.0*A4+5.0*A5*t)))
            if(abs(df)<epsDF):break
            t-=0.75*f/df
            f=Bs+t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))-theta
            nn+=1
        ###print("nIter = ",nn)    
        vtt.append(t)
    return np.array(vtt)    

'''
function:'uaxLambda'
--Calculates the maximum shape parameter lambda for the uniaxial interpolated data
--Input:
----data=the overall data structure of the material
--Output:
--four (or two) maximum lambda values
--lists of precalculated values of tangents
--Note: the angle along the theta-axis of each directional data set is in degrees 
'''
def uaxLambda(data):
    vvPolygon={'lambdaSTmax':0,'lambdaRTmax':0,'pointsST':0,'tanST':0,'pointsRT':0,'tanRT':0,
               'sTb':data['sTb'],'rTb':data['rTb'],'shapeUAX':data['shapeUAX'],'assym':data['assym'],'name':data['name'],'shapeBAX':data['shapeBAX']}
    ##vThetaT=np.array([data['thetaT'][theta] for theta in data['thetaT']])
    vTheta=np.array([float(theta) for theta in data['thetaT']]) ##use degrees
    nTheta=vTheta.shape[0]
    pointsS=np.zeros((2,nTheta))
    pointsS[0,:]=vTheta
    pointsS[1,:]=[data['sT'][k] for k in data['sT']]
    pointsR=np.zeros((2,nTheta))
    pointsR[0,:]=vTheta
    pointsR[1,:]=[data['rT'][k] for k in data['rT']]
    tanS=np.zeros((2,nTheta))
    tanR=np.zeros((2,nTheta))
    ##dx=(0.5*np.pi)/(nTheta-1)
    dx=90.0/(nTheta-1.0)
    Ldata=data['shapeTAN']
    tanS[:,0]=[1.0,0.0];tanR[:,0]=[1.0,0.0]
    jj=1
    lsMax=10**6;lrMax=10**6
    while(jj<nTheta-1):
        ty=(1.0-Ldata)*(pointsS[1,jj]-pointsS[1,jj-1])+Ldata*(pointsS[1,jj+1]-pointsS[1,jj])
        nt=np.sqrt(dx*dx+ty*ty)
        tanS[:,jj]=[dx/nt,ty/nt]
        lsMax=min(lsMax,(0.5*dx)/(tanS[0,jj]+tanS[0,jj-1]))
        ty=(1.0-Ldata)*(pointsR[1,jj]-pointsR[1,jj-1])+Ldata*(pointsR[1,jj+1]-pointsR[1,jj])
        nt=np.sqrt(dx*dx+ty*ty)
        tanR[:,jj]=[dx/nt,ty/nt]
        lrMax=min(lrMax,(0.5*dx)/(tanR[0,jj]+tanR[0,jj-1]))
        jj+=1    
    tanS[:,-1]=[1.0,0.0];tanR[:,-1]=[1.0,0.0]
    ##scale=180.0/np.pi
    ##vvPolygon['lambdaSTmax']=scale*min(lsMax,(0.5*dx)/(tanST[0,jj]+tanST[0,jj-1]))    
    ##vvPolygon['lambdaRTmax']=scale*min(lrMax,(0.5*dx)/(tanRT[0,jj]+tanRT[0,jj-1]))
    vvPolygon['lambdaSTmax']=min(lsMax,(0.5*dx)/(tanS[0,jj]+tanS[0,jj-1]))    
    vvPolygon['lambdaRTmax']=min(lrMax,(0.5*dx)/(tanR[0,jj]+tanR[0,jj-1]))
    vvPolygon['pointsST']=pointsS;vvPolygon['tanST']=tanS
    vvPolygon['pointsRT']=pointsR;vvPolygon['tanRT']=tanR
    if(not data['assym']):
        pass
        ##return vvPolygon
    vTheta=np.array([float(theta) for theta in data['thetaC']]) ##use degrees
    nTheta=vTheta.shape[0]
    pointsS=np.zeros((2,nTheta))
    pointsS[0,:]=vTheta
    pointsS[1,:]=[data['sC'][k] for k in data['sC']]
    pointsR=np.zeros((2,nTheta))
    pointsR[0,:]=vTheta
    pointsR[1,:]=[data['rC'][k] for k in data['rC']]
    tanS=np.zeros((2,nTheta))
    tanR=np.zeros((2,nTheta))
    dx=90.0/(nTheta-1.0)
    tanS[:,0]=[1.0,0.0];tanR[:,0]=[1.0,0.0]
    jj=1
    lsMax=10**6;lrMax=10**6
    while(jj<nTheta-1):
        ty=(1.0-Ldata)*(pointsS[1,jj]-pointsS[1,jj-1])+Ldata*(pointsS[1,jj+1]-pointsS[1,jj])
        nt=np.sqrt(dx*dx+ty*ty)
        tanS[:,jj]=[dx/nt,ty/nt]
        lsMax=min(lsMax,(0.5*dx)/(tanS[0,jj]+tanS[0,jj-1]))
        ty=(1.0-Ldata)*(pointsR[1,jj]-pointsR[1,jj-1])+Ldata*(pointsR[1,jj+1]-pointsR[1,jj])
        nt=np.sqrt(dx*dx+ty*ty)
        tanR[:,jj]=[dx/nt,ty/nt]
        lrMax=min(lrMax,(0.5*dx)/(tanR[0,jj]+tanR[0,jj-1]))
        jj+=1    
    tanS[:,-1]=[1.0,0.0];tanR[:,-1]=[1.0,0.0]
    vvPolygon['lambdaSCmax']=min(lsMax,(0.5*dx)/(tanS[0,jj]+tanS[0,jj-1]))    
    vvPolygon['lambdaRCmax']=min(lrMax,(0.5*dx)/(tanR[0,jj]+tanR[0,jj-1]))
    vvPolygon['pointsSC']=pointsS;vvPolygon['tanSC']=tanS
    vvPolygon['pointsRC']=pointsR;vvPolygon['tanRC']=tanR
    vvPolygon['sCb']=data['sCb'];vvPolygon['rCb']=data['rCb']    
    return vvPolygon


'''
function:'protoTheta'
--Calculates the values of the Bezier parameters 't' corresponding to a list of theta-angles  
--Input:
----data=precalculated data returned by function 'uaxLambda'
--Output:
--'thetaTList'=the correspondence on the gammaT curve (with top and bottom)
--'thetaCList'=the correspondence on the gammaC curve (with top and bottom)
''' 
#  
def protoTheta(data,nSections=10):
    ##print("protoTheta: using nSections = {}".format(nSections))
    vTheta=np.linspace(0.0,45.0,nSections)
    thetaData=data['pointsST'][0,:];nTheta=thetaData.shape[0]
    thetaTList=[]
    jj=0
    LshapeS=data['shapeUAX']*data['lambdaSTmax']
    LshapeR=data['shapeUAX']*data['lambdaRTmax']
    while(True):
        thetaA=thetaData[jj];thetaB=thetaData[jj+1]
        if(thetaB>45.1):break
        dd={'idxTop':(jj,jj+1),'vThetaTop':[],'vSttTop':[],'vRttTop':[],
            'idxBottom':(nTheta-jj-2,nTheta-jj-1),'vThetaBottom':[],'vSttBottom':[],'vRttBottom':[]}
        vv=vTheta[np.argwhere((vTheta>=thetaA)*(vTheta<thetaB))]
        dd['vThetaTop']=vv.reshape(vv.shape[0])
        dd['vSttTop']=curveSegTheta(dd['vThetaTop'],thetaA,data['tanST'][0,jj],thetaB,data['tanST'][0,jj+1],LshapeS)
        dd['vRttTop']=curveSegTheta(dd['vThetaTop'],thetaA,data['tanRT'][0,jj],thetaB,data['tanRT'][0,jj+1],LshapeR)
        dd['vThetaBottom']=90-dd['vThetaTop']
        dd['vSttBottom']=curveSegTheta(dd['vThetaBottom'],90-thetaB,data['tanST'][0,nTheta-jj-2],90-thetaA,data['tanST'][0,nTheta-jj-1],LshapeS)
        dd['vRttBottom']=curveSegTheta(dd['vThetaBottom'],90-thetaB,data['tanRT'][0,nTheta-jj-2],90-thetaA,data['tanRT'][0,nTheta-jj-1],LshapeR)
        thetaTList.append(dd)
        jj+=1
    thetaCList=[]
    if(not data['assym']):
        return thetaTList,thetaCList
    thetaData=data['pointsSC'][0,:];nTheta=thetaData.shape[0]
    jj=0
    LshapeS=data['shapeUAX']*data['lambdaSCmax']
    LshapeR=data['shapeUAX']*data['lambdaRCmax']
    while(True):
        thetaA=thetaData[jj];thetaB=thetaData[jj+1]
        if(thetaB>45.1):break
        dd={'idxTop':(jj,jj+1),'vThetaTop':[],'vSttTop':[],'vRttTop':[],
            'idxBottom':(nTheta-jj-2,nTheta-jj-1),'vThetaBottom':[],'vSttBottom':[],'vRttBottom':[]}
        vv=vTheta[np.argwhere((vTheta>=thetaA)*(vTheta<thetaB))]
        dd['vThetaTop']=vv.reshape(vv.shape[0])
        dd['vSttTop']=curveSegTheta(dd['vThetaTop'],thetaA,data['tanSC'][0,jj],thetaB,data['tanSC'][0,jj+1],LshapeS)
        dd['vRttTop']=curveSegTheta(dd['vThetaTop'],thetaA,data['tanRC'][0,jj],thetaB,data['tanRC'][0,jj+1],LshapeR)
        dd['vThetaBottom']=90-dd['vThetaTop']
        dd['vSttBottom']=curveSegTheta(dd['vThetaBottom'],90-thetaB,data['tanSC'][0,nTheta-jj-2],90-thetaA,data['tanSC'][0,nTheta-jj-1],LshapeS)
        dd['vRttBottom']=curveSegTheta(dd['vThetaBottom'],90-thetaB,data['tanRC'][0,nTheta-jj-2],90-thetaA,data['tanRC'][0,nTheta-jj-1],LshapeR)
        thetaCList.append(dd)
        jj+=1        
    return thetaTList,thetaCList    


def lambdaMax(bOne,tOne,bTwo,tTwo):
    bb=bTwo-bOne
    mdot=np.sum(tOne*tTwo)
    ddet=1.0-mdot*mdot
    t1=0.5*np.sum(bb*(tOne-mdot*tTwo))/ddet
    t2=0.5*np.sum(bb*(tTwo-mdot*tOne))/ddet
    #return min(t1,t2)
    return t1,t2

def lambdaMaxDebug(bOne,tOne,bTwo,tTwo):
    bb=bTwo-bOne
    mdot=np.sum(tOne*tTwo)
    ddet=1.0-mdot*mdot
    t1=0.5*np.sum(bb*(tOne-mdot*tTwo))/ddet
    t2=0.5*np.sum(bb*(tTwo-mdot*tOne))/ddet
    print("bOne,bTwo = \n",bOne,bTwo)
    print("TOne,tTwo = \n",tOne,tTwo)
    print("mdot = ",mdot)
    print("t1,t2 and min = ",t1,t2,min(t1,t2))
    return    

def vectorProd(a,b):
    vv=np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
    norm=np.sqrt(np.sum(vv**2))
    return vv/norm

'''
function:'protoData'
--Interpolates plane sections and calculates the  maximum shape parameter lambda for the yield surface
--Input:
----data=precalculated data returned by function 'uaxLambda'
--Output:
--the max lambda value (shape parameter)
--lists of precalculated values of points and tangents along the uniaxial 3D curves (tension/compression) 
  and at the two balanced-biaxial yielding points  
'''    
def protoData(uaxData,nSections=10):
    ###Generate the sequence of points and tangents for Bez5YS proto-model
    ##print("protoData: nSections = {}".format(nSections))
    thetaTList,thetaCList=protoTheta(uaxData,nSections)
    ###print("uaxData");print(uaxData)
    ##print("thetaTList");print(thetaTList)
    ##print("thetaCList");print(thetaCList)
    rad=np.pi/180.0
    LshapeST=uaxData['shapeUAX']*uaxData['lambdaSTmax']
    LshapeRT=uaxData['shapeUAX']*uaxData['lambdaRTmax']
    LshapeSC=uaxData['shapeUAX']*uaxData['lambdaSTmax']
    LshapeRC=uaxData['shapeUAX']*uaxData['lambdaRTmax']
    if(uaxData['assym']):
        LshapeSC=uaxData['shapeUAX']*uaxData['lambdaSCmax']
        LshapeRC=uaxData['shapeUAX']*uaxData['lambdaRCmax']
    vPatch=[]
    BS=np.zeros((2,1));TS=np.zeros((2,1));BE=np.zeros((2,1));TE=np.zeros((2,1))
    dGamma=np.zeros((3,))
    for dd in thetaTList:##for each segment of the interpolated directional properties 
        vThetaTop=rad*dd['vThetaTop'] ##the list of angles in the segment (in radians)
        vV=np.array([[-np.sin(2*t),np.sin(2*t),2*np.cos(2*t)] for t in vThetaTop]) ##vectors of PI-plane normals
        vTT=dd['vSttTop'] 
        BS[:]=uaxData['pointsST'][:,[dd['idxTop'][0]]]
        TS[:]=uaxData['tanST'][:,[dd['idxTop'][0]]]
        BE[:]=uaxData['pointsST'][:,[dd['idxTop'][1]]]
        TE[:]=uaxData['tanST'][:,[dd['idxTop'][1]]]
        ##print(BS,TS,BE,TE)
        vSTtop,vDSTtop=curveSegFDF(BS,TS,BE,TE,LshapeST,vTT)
        vDSTtop=vDSTtop[1,:]/(rad*vDSTtop[0,:])
        ##print("vSTtop=\n",vSTtop);print("vDSTtop=\n",vDSTtop)
        vTT=dd['vRttTop'] 
        BS[:]=uaxData['pointsRT'][:,[dd['idxTop'][0]]]
        TS[:]=uaxData['tanRT'][:,[dd['idxTop'][0]]]
        BE[:]=uaxData['pointsRT'][:,[dd['idxTop'][1]]]
        TE[:]=uaxData['tanRT'][:,[dd['idxTop'][1]]]
        vRTtop=curveSegF(BS,TS,BE,TE,LshapeRT,vTT)
        vWtop=np.array([[vRTtop[1,j]+(np.sin(t))**2,vRTtop[1,j]+(np.cos(t))**2,-np.cos(t)*np.sin(t)] for j,t in enumerate(vThetaTop)])
        vGammaTtop=np.array([[(np.cos(t))**2,(np.sin(t))**2,np.cos(t)*np.sin(t)] for t in vThetaTop])
        vDGammaTtop=np.array([[-np.sin(2*t),np.sin(2*t),np.cos(2*t)] for t in vThetaTop])
        vTT=dd['vSttBottom'] 
        BS[:]=uaxData['pointsST'][:,[dd['idxBottom'][0]]]
        TS[:]=uaxData['tanST'][:,[dd['idxBottom'][0]]]
        BE[:]=uaxData['pointsST'][:,[dd['idxBottom'][1]]]
        TE[:]=uaxData['tanST'][:,[dd['idxBottom'][1]]]
        vSTbottom,vDSTbottom=curveSegFDF(BS,TS,BE,TE,LshapeST,vTT)
        vDSTbottom=vDSTbottom[1,:]/(rad*vDSTbottom[0,:])
        vTT=dd['vRttBottom'] 
        BS[:]=uaxData['pointsRT'][:,[dd['idxBottom'][0]]]
        TS[:]=uaxData['tanRT'][:,[dd['idxBottom'][0]]]
        BE[:]=uaxData['pointsRT'][:,[dd['idxBottom'][1]]]
        TE[:]=uaxData['tanRT'][:,[dd['idxBottom'][1]]]
        vRTbottom=curveSegF(BS,TS,BE,TE,LshapeRT,vTT)
        vWbottom=np.array([[vRTbottom[1,j]+(np.sin(t))**2,vRTbottom[1,j]+(np.cos(t))**2,np.cos(t)*np.sin(t)] for j,t in enumerate(0.5*np.pi-vThetaTop)])
        vGammaTbottom=np.array([[(np.sin(t))**2,(np.cos(t))**2,-np.cos(t)*np.sin(t)] for t in vThetaTop])
        vDGammaTbottom=np.array([[np.sin(2*t),-np.sin(2*t),-np.cos(2*t)] for t in vThetaTop])
        for kk in range(len(vThetaTop)):
            vPoints=np.zeros((3,6));vTangents=np.zeros((3,6))
            vPoints[:,0]=vSTtop[1,kk]*vGammaTtop[kk,:]
            dGamma[:]=vDSTtop[kk]*vGammaTtop[kk]+vSTtop[1,kk]*vDGammaTtop[kk]
            vTangents[:,0]=vectorProd(vV[kk],vectorProd(vWtop[kk],dGamma))
            vPoints[:,1]=[uaxData['sTb'],uaxData['sTb'],0.0]
            vTangents[:,1]=vectorProd(vV[kk],np.array([1.0,uaxData['rTb'],0])/np.sqrt(1.0+uaxData['rTb']**2))
            vPoints[:,2]=vSTbottom[1,kk]*vGammaTbottom[kk,:]
            dGamma[:]=-vDSTbottom[kk]*vGammaTbottom[kk]+vSTbottom[1,kk]*vDGammaTbottom[kk]
            vTangents[:,2]=vectorProd(vV[kk],vectorProd(vWbottom[kk],dGamma))
            vPatch.append({'Points':vPoints,'Tangents':vTangents})
    #lbdMax=10.0**3
    vLbdMax=1.0e+3*np.ones(6)
    if(not uaxData['assym']):
        for j in range(len(vPatch)):
            vPatch[j]['Points'][:,3]=-vPatch[j]['Points'][:,0]
            vPatch[j]['Points'][:,4]=-vPatch[j]['Points'][:,1]
            vPatch[j]['Points'][:,5]=-vPatch[j]['Points'][:,2]
            vPatch[j]['Tangents'][:,3]=-vPatch[j]['Tangents'][:,0]
            vPatch[j]['Tangents'][:,4]=-vPatch[j]['Tangents'][:,1]
            vPatch[j]['Tangents'][:,5]=-vPatch[j]['Tangents'][:,2]
            for kk in range(5):
                #lbd=lambdaMax(vPatch[j]['Points'][:,kk],vPatch[j]['Tangents'][:,kk],vPatch[j]['Points'][:,kk+1],vPatch[j]['Tangents'][:,kk+1])
                lbd1,lbd2=lambdaMax(vPatch[j]['Points'][:,kk],vPatch[j]['Tangents'][:,kk],vPatch[j]['Points'][:,kk+1],vPatch[j]['Tangents'][:,kk+1])
                if(lbd1>0 and lbd2>0):
                    #lbdMax=min(lbdMax,lbd)
                    vLbdMax[kk]=min(vLbdMax[kk],lbd1)
                    vLbdMax[kk+1]=min(vLbdMax[kk+1],lbd2)
                else:
                    print("WARNING from 'protoData':")
                    print("--Negative plane section shape parameter lambdaMax = ",lbd)
                    print("--Data is not consistent with a convex model")
                    print("Calculations aborted");exit()
            lbd1,lbd2=lambdaMax(vPatch[j]['Points'][:,5],vPatch[j]['Tangents'][:,5],vPatch[j]['Points'][:,0],vPatch[j]['Tangents'][:,0])
            if(lbd1>0 and lbd2>0):
                vLbdMax[5]=min(vLbdMax[5],lbd1)
                vLbdMax[0]=min(vLbdMax[0],lbd2)
            else:
                print("WARNING from 'protoData':")
                print("--Negative plane section shape parameter lambdaMax = ",lbd)
                print("--Data is not consistent with a convex model")
                print("Calculations aborted");exit()        
        ##return lbdMax,vPatch
        vLbdMax=np.array([vLbdMax[0],vLbdMax[1],vLbdMax[0],vLbdMax[0],vLbdMax[1],vLbdMax[0]])
        if(uaxData['shapeBAX']==0):
            vLbdMax=np.min(vLbdMax)*np.ones(6)
        return vLbdMax,vPatch
    jj=0    ##;print("thetaTList\n",thetaTList);print("thetaCList\n",thetaCList)
    for dd in thetaCList:##for each segment of the interpolated directional properties 
        vThetaTop=rad*dd['vThetaTop'] ##the list of angles in the segment (in radians)
        vV=np.array([[-np.sin(2*t),np.sin(2*t),2*np.cos(2*t)] for t in vThetaTop]) ##vectors of PI-plane normals
        vTT=dd['vSttTop'] 
        BS[:]=uaxData['pointsSC'][:,[dd['idxTop'][0]]]
        TS[:]=uaxData['tanSC'][:,[dd['idxTop'][0]]]
        BE[:]=uaxData['pointsSC'][:,[dd['idxTop'][1]]]
        TE[:]=uaxData['tanSC'][:,[dd['idxTop'][1]]]
        vSTtop,vDSTtop=curveSegFDF(BS,TS,BE,TE,LshapeSC,vTT)
        vDSTtop=vDSTtop[1,:]/(rad*vDSTtop[0,:])
        vTT=dd['vRttTop'] 
        BS[:]=uaxData['pointsRC'][:,[dd['idxTop'][0]]]
        TS[:]=uaxData['tanRC'][:,[dd['idxTop'][0]]]
        BE[:]=uaxData['pointsRC'][:,[dd['idxTop'][1]]]
        TE[:]=uaxData['tanRC'][:,[dd['idxTop'][1]]]
        vRTtop=curveSegF(BS,TS,BE,TE,LshapeRC,vTT)
        vWtop=np.array([[vRTtop[1,j]+(np.sin(t))**2,vRTtop[1,j]+(np.cos(t))**2,-np.cos(t)*np.sin(t)] for j,t in enumerate(vThetaTop)])###surface normal points downwards
        vGammaTtop=np.array([[-(np.cos(t))**2,-(np.sin(t))**2,-np.cos(t)*np.sin(t)] for t in vThetaTop])##This is actually the bottom part 
        vDGammaTtop=np.array([[np.sin(2*t),-np.sin(2*t),-np.cos(2*t)] for t in vThetaTop])
        vTT=dd['vSttBottom'] 
        BS[:]=uaxData['pointsSC'][:,[dd['idxBottom'][0]]]
        TS[:]=uaxData['tanSC'][:,[dd['idxBottom'][0]]]
        BE[:]=uaxData['pointsSC'][:,[dd['idxBottom'][1]]]
        TE[:]=uaxData['tanSC'][:,[dd['idxBottom'][1]]]
        vSTbottom,vDSTbottom=curveSegFDF(BS,TS,BE,TE,LshapeSC,vTT)
        vDSTbottom=vDSTbottom[1,:]/(rad*vDSTbottom[0,:])
        vTT=dd['vRttBottom'] 
        BS[:]=uaxData['pointsRC'][:,[dd['idxBottom'][0]]]
        TS[:]=uaxData['tanRC'][:,[dd['idxBottom'][0]]]
        BE[:]=uaxData['pointsRC'][:,[dd['idxBottom'][1]]]
        TE[:]=uaxData['tanRC'][:,[dd['idxBottom'][1]]]
        vRTbottom=curveSegF(BS,TS,BE,TE,LshapeRC,vTT)
        vWbottom=np.array([[vRTbottom[1,j]+(np.sin(t))**2,vRTbottom[1,j]+(np.cos(t))**2,np.cos(t)*np.sin(t)] for j,t in enumerate(0.5*np.pi-vThetaTop)])
        vGammaTbottom=np.array([[-(np.sin(t))**2,-(np.cos(t))**2,np.cos(t)*np.sin(t)] for t in vThetaTop])##This is actually the top part 
        vDGammaTbottom=np.array([[-np.sin(2*t),np.sin(2*t),np.cos(2*t)] for t in vThetaTop])
        for kk in range(len(vThetaTop)):
            vPatch[jj]['Points'][:,3]=vSTtop[1,kk]*vGammaTtop[kk,:]
            dGamma[:]=vDSTtop[kk]*vGammaTtop[kk]+vSTtop[1,kk]*vDGammaTtop[kk]
            vPatch[jj]['Tangents'][:,3]=vectorProd(vV[kk],vectorProd(vWtop[kk],dGamma))
            vPatch[jj]['Points'][:,4]=[-uaxData['sCb'],-uaxData['sCb'],0.0]
            vPatch[jj]['Tangents'][:,4]=vectorProd(vV[kk],np.array([-1.0,-uaxData['rCb'],0])/np.sqrt(1.0+uaxData['rCb']**2))
            vPatch[jj]['Points'][:,5]=vSTbottom[1,kk]*vGammaTbottom[kk,:]
            dGamma[:]=-vDSTbottom[kk]*vGammaTbottom[kk]+vSTbottom[1,kk]*vDGammaTbottom[kk]
            vPatch[jj]['Tangents'][:,5]=vectorProd(vV[kk],vectorProd(vWbottom[kk],dGamma))
            for qq in range(5):
                lbd1,lbd2=lambdaMax(vPatch[jj]['Points'][:,qq],vPatch[jj]['Tangents'][:,qq],vPatch[jj]['Points'][:,qq+1],vPatch[jj]['Tangents'][:,qq+1])
                if(lbd1>0 and lbd2>0):
                    vLbdMax[qq]=min(vLbdMax[qq],lbd1)
                    vLbdMax[qq+1]=min(vLbdMax[qq+1],lbd2)    
                else:
                    print("WARNING from 'protoData':")
                    print("--Negative plane section shape parameter lambdaMax = ",lbd)
                    print("--Data is not consistent with a convex model")
                    print("Calculations aborted");exit()
            lbd1,lbd2=lambdaMax(vPatch[jj]['Points'][:,5],vPatch[jj]['Tangents'][:,5],vPatch[jj]['Points'][:,0],vPatch[jj]['Tangents'][:,0])
            if(lbd1>0 and lbd2>0):
                vLbdMax[5]=min(vLbdMax[5],lbd1)
                vLbdMax[0]=min(vLbdMax[0],lbd2)    
            else:
                print("WARNING from 'protoData':")
                print("--Negative plane section shape parameter lambdaMax = ",lbd)
                print("--Data is not consistent with a convex model")
                print("Calculations aborted");exit()
            jj+=1
    aa=min(vLbdMax[0],vLbdMax[2]);bb=min(vLbdMax[3],vLbdMax[5])        
    #vLbdMax=np.array([vLbdMax[0],vLbdMax[1],vLbdMax[0],vLbdMax[2],vLbdMax[3],vLbdMax[2]])
    vLbdMax=np.array([aa,vLbdMax[1],aa,bb,vLbdMax[4],bb])     
    if(uaxData['shapeBAX']==0):
            vLbdMax=np.min(vLbdMax)*np.ones(6)
    return vLbdMax,vPatch


'''
function:'protoDataPoints'
--Calculates a list of points on the proto-model of the yield surface 
'''
def protoDataPoints(lbd,shapeVBAX,uaxData,nSections,nPointsSegment):
    ##print("Using: nSections={}, nPoints/SectionSegment={}".format(nSections,nPointsSegment))
    lbd2,data=protoData(uaxData,nSections) 
    lbd=np.min((lbd2,lbd),axis=0)*shapeVBAX
    #lbd=lbd*shapeBAX
    vYSpoints=np.zeros(((nSections-1)*6*nPointsSegment,3))
    jj=0
    for dd in data:
        for kk in range(5):
            vYSpoints[jj:jj+nPointsSegment]=curveSeg3(
            dd['Points'][:,kk].reshape(3,1),dd['Tangents'][:,kk].reshape(3,1),dd['Points'][:,kk+1].reshape(3,1),dd['Tangents'][:,kk+1].reshape(3,1),lbd[kk],lbd[kk+1],False,nPointsSegment).T
            jj+=nPointsSegment
        vYSpoints[jj:jj+nPointsSegment]=curveSeg3(
        dd['Points'][:,5].reshape(3,1),dd['Tangents'][:,5].reshape(3,1),dd['Points'][:,0].reshape(3,1),dd['Tangents'][:,0].reshape(3,1),lbd[5],lbd[0],False,nPointsSegment).T
        jj+=nPointsSegment    
    return vYSpoints



'''
function:'genConstraintsPoints2D'
--Generates a list of points and corresponding (unit) tangent vectors (for plane stress) 
  on the unit sphere of the (u1,u2,u3) space 
--Returns a (nConstraints,15)-array with each row storing a point and four corresponding tangent vectors 
--Note: obsolete (replaced by optimized version 'genConstraintsPoints2DOpt') 
'''
def genConstraintsPoints2D(nPoints):
    ##delta=np.sqrt(4*np.pi/nPoints)
    ##delta=2.0*np.pi/nPoints
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1);vs1=np.sin(vt1)
    vP=[]
    vv=np.zeros(3)
    g1=np.zeros(3);g2=np.zeros(3);g3=np.zeros(3);g4=np.zeros(3)
    for kk in range(len(vt1)):
        ct1=vc1[kk];st1=vs1[kk]
        vt2=np.linspace(0,2*np.pi,int(nPoints*st1)+1)
        vc2=np.cos(vt2);vs2=np.sin(vt2)
        for jj in range(len(vt2)):
            ct2=vc2[jj];st2=vs2[jj]
            vv[:]=np.array([st1*ct2,st1*st2,ct1])
            g1[:]=np.array([ct1*ct2,ct1*st2,-st1])
            g2[:]=np.array([-st2,ct2,0.0])
            g3[:]=g1+g2;g3[:]/=np.sqrt(np.sum(g3*g3))
            g4[:]=g1-g2;g4[:]/=np.sqrt(np.sum(g4*g4))
            vP.append(np.array([vv,g1,g2,g3,g4]).reshape((15,)))
    return np.array(vP)
    
def genConstraintsPoints2DOpt(nPoints,number=False):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)
    if(0):
        tt=[(1,1),(-1,1),(1,0.5),(0.5,1),(-1,0.5),(-0.5,1),(1,0.25),(-1,0.25),
            (1,0.75),(-1,0.75),(0.75,1),(0.25,1),(-0.25,1),(-0.75,1),
            (1,0.1),(-1,0.1),(0.1,1),(-0.1,1),(1,0.9),(-1,0.9),(0.9,1),(-0.9,1),
            (0.4,1),(-0.4,1),(1,0.4),(-1,0.4),(0.6,1),(-0.6,1),(1,0.6),(-1,0.6)]
            ##(0.3,1),(-0.3,1),(1,0.3),(-1,0.3),(0.35,1),(-0.35,1),(1,0.35),(-1,0.35)] 
    ##Use the next denser grid of tangents when convexity is expected but not achieved
    ##(because of limited numerical precision)    
    tt=[(np.cos(t),np.sin(t)) for t in np.linspace(0,np.pi,51)[1:-1]]        
    nVec=len(tt)
    if(number): return (np.sum(vN2)+1,3*(nVec+3))
    vP=np.zeros((np.sum(vN2)+1,3*(nVec+3)))
    vN2[:]+=1
    vP[0,2]=1.0
    vP[0,3]=1.0  #g1[0]
    vP[0,7]=1.0  ##g2[1]
    i=9
    for t in tt:
        gg=t[0]*vP[0,3:6]+t[1]*vP[0,6:9]
        gg[:]/=np.sqrt(np.sum(gg*gg,axis=0))
        vP[0,i]=gg[0];vP[0,i+1]=gg[1];vP[0,i+2]=gg[2]
        i+=3
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,2*np.pi,N2)
        vc2=np.cos(vt2[0:N2-1]);vs2=np.sin(vt2[0:N2-1])
        ##print("kk={}, N2={}, jN2={}, jj={}".format(kk,N2,jN2,jj))
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        vP[jj:jN2,3]=-vs2  ##g1[0]
        vP[jj:jN2,4]=vc2
        vP[jj:jN2,6]=ct1*vc2  #g2[0]
        vP[jj:jN2,7]=ct1*vs2
        vP[jj:jN2,8]=-st1
        i=9
        for t in tt:
            gg=t[0]*vP[jj:jN2,3:6]+t[1]*vP[jj:jN2,6:9]
            vMod=np.sqrt(np.sum(gg*gg,axis=1))
            gg[:,0]/=vMod;gg[:,1]/=vMod;gg[:,2]/=vMod
            vP[jj:jN2,i]=gg[:,0];vP[jj:jN2,i+1]=gg[:,1];vP[jj:jN2,i+2]=gg[:,2]
            i+=3
        jj=jN2
    return vP        

'''
function: genConstraintsPoints2DOptPoints
--the same as 'genConstraintsPoints2DOpt' but returns only points (no tangents)
--used for testing/debugging
'''
def genConstraintsPoints2DOptPoints(nPoints):
    N1=int(0.25*nPoints)
    vt1=np.linspace(0,np.pi/2,N1+1)
    vc1=np.cos(vt1[1:]);vs1=np.sin(vt1[1:])
    vN2=np.int_(nPoints*vs1)     
    vP=np.zeros((np.sum(vN2)+1,3))
    vN2[:]+=1
    vP[0,2]=1.0
    jj=1
    for kk in range(0,len(vt1)-1):
        ct1=vc1[kk];st1=vs1[kk]
        N2=vN2[kk];jN2=jj+N2-1
        vt2=np.linspace(0,2*np.pi,N2)
        vc2=np.cos(vt2[0:N2-1]);vs2=np.sin(vt2[0:N2-1])
        #print("kk={}, N2={}, jN2={}, jj={}".format(kk,N2,jN2,jj))
        vP[jj:jN2,0]=st1*vc2
        vP[jj:jN2,1]=st1*vs2
        vP[jj:jN2,2]=ct1
        jj=jN2
    return vP        

'''
function:'genConstraints2D'
--Calculates the constraints based on points (locations) and tangents 
--Returns: matrix of constraints and vector of bounds 
'''
def genConstraints2D(ddMon,degQ,nP,nQ,cP0,cP1,cQ0,cQ1,nEquator,epsilon=0.01):
    print('generating constraints with nEquator = {} and epsilon = {} ...'.format(nEquator,epsilon))
    degPm1=degQ-2;degQm1=degQ-1;nCol=nQ+nP
    vP=ddMon['vP']
    vDD11P=ddMon['vH11P'];vDD12P=ddMon['vH12P'];vDD13P=ddMon['vH13P']
    vDD22P=ddMon['vH22P'];vDD23P=ddMon['vH23P'];vDD33P=ddMon['vH33P']
    vCD11P=ddMon['vCH11P'];vCD12P=ddMon['vCH12P'];vCD13P=ddMon['vCH13P']
    vCD22P=ddMon['vCH22P'];vCD23P=ddMon['vCH23P'];vCD33P=ddMon['vCH33P']
    vQ=ddMon['vQ']
    vDD11Q=ddMon['vH11Q'];vDD12Q=ddMon['vH12Q'];vDD13Q=ddMon['vH13Q']
    vDD22Q=ddMon['vH22Q'];vDD23Q=ddMon['vH23Q'];vDD33Q=ddMon['vH33Q']
    vCD11Q=ddMon['vCH11Q'];vCD12Q=ddMon['vCH12Q'];vCD13Q=ddMon['vCH13Q']
    vCD22Q=ddMon['vCH22Q'];vCD23Q=ddMon['vCH23Q'];vCD33Q=ddMon['vCH33Q']
    ##vUG=genConstraintsPoints2D()
    vUG=genConstraintsPoints2DOpt(nEquator)
    nPoints=vUG.shape[0]
    nVec=(vUG.shape[1]-3)//3##number of tangent vectors 
    nCons=nPoints*nVec##number of constraints
    GG=np.zeros((nPoints,6*nVec))##'6' is specific to 2D
    vPoints=np.zeros((nPoints,3*nVec))##'3' is specific to 2D
    jj=3;jG=0
    for kk in range(nVec):
        GG[:,jG]=vUG[:,jj]**2;jG+=1
        GG[:,jG]=vUG[:,jj+1]**2;jG+=1
        GG[:,jG]=vUG[:,jj+2]**2;jG+=1
        GG[:,jG]=2*vUG[:,jj]*vUG[:,jj+1];jG+=1
        GG[:,jG]=2*vUG[:,jj]*vUG[:,jj+2];jG+=1
        GG[:,jG]=2*vUG[:,jj+1]*vUG[:,jj+2];jG+=1
        vPoints[:,jj-3:jj]=vUG[:,0:3]
        jj+=3
    GG=GG.reshape((nCons,6))    
    vPoints=vPoints.reshape((nCons,3))
    vPowers1=np.zeros((nCons,nP))
    vPowers2=np.zeros((nCons,nP))
    vPowers3=np.zeros((nCons,nP))
    vv=np.zeros((nCons,nCol))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vP[k][0]
        vPowers2[:,k]=vPoints[:,1]**vP[k][1]
        vPowers3[:,k]=vPoints[:,2]**vP[k][2]
    vv[:,0:nP]=degPm1*vPowers1*vPowers2*vPowers3
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD11P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD11P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD11P[k][2]
    vv[:,0:nP]-=vCD11P[:]*vPowers1*vPowers2*vPowers3*GG[:,0].reshape((nCons,1))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD22P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD22P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD22P[k][2]
    vv[:,0:nP]-=vCD22P[:]*vPowers1*vPowers2*vPowers3*GG[:,1].reshape((nCons,1))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD33P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD33P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD33P[k][2]
    vv[:,0:nP]-=vCD33P[:]*vPowers1*vPowers2*vPowers3*GG[:,2].reshape((nCons,1))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD12P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD12P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD12P[k][2]
    vv[:,0:nP]-=vCD12P[:]*vPowers1*vPowers2*vPowers3*GG[:,3].reshape((nCons,1))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD13P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD13P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD13P[k][2]
    vv[:,0:nP]-=vCD13P[:]*vPowers1*vPowers2*vPowers3*GG[:,4].reshape((nCons,1))
    for k in range(nP):
        vPowers1[:,k]=vPoints[:,0]**vDD23P[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD23P[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD23P[k][2]
    vv[:,0:nP]-=vCD23P[:]*vPowers1*vPowers2*vPowers3*GG[:,5].reshape((nCons,1))
    vPowers1=np.zeros((nCons,nQ))
    vPowers2=np.zeros((nCons,nQ))
    vPowers3=np.zeros((nCons,nQ))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vQ[k][0]
        vPowers2[:,k]=vPoints[:,1]**vQ[k][1]
        vPowers3[:,k]=vPoints[:,2]**vQ[k][2]
    vv[:,nP:]=degQm1*vPowers1*vPowers2*vPowers3
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD11Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD11Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD11Q[k][2]
    vv[:,nP:]-=vCD11Q[:]*vPowers1*vPowers2*vPowers3*GG[:,0].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD22Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD22Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD22Q[k][2]
    vv[:,nP:]-=vCD22Q[:]*vPowers1*vPowers2*vPowers3*GG[:,1].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD33Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD33Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD33Q[k][2]
    vv[:,nP:]-=vCD33Q[:]*vPowers1*vPowers2*vPowers3*GG[:,2].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD12Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD12Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD12Q[k][2]
    vv[:,nP:]-=vCD12Q[:]*vPowers1*vPowers2*vPowers3*GG[:,3].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD13Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD13Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD13Q[k][2]
    vv[:,nP:]-=vCD13Q[:]*vPowers1*vPowers2*vPowers3*GG[:,4].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD23Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD23Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD23Q[k][2]
    vv[:,nP:]-=vCD23Q[:]*vPowers1*vPowers2*vPowers3*GG[:,5].reshape((nCons,1))
    MCC=np.zeros((nCons,nCol-4))
    MUB=np.zeros((nCons,1))
    MCC[:,0:nP-2]=vv[:,2:nP];MCC[:,nP-2:]=vv[:,nP+2:]
    MUB[:,0]=1.0-epsilon-(cP0*vv[:,0]+cP1*vv[:,1]+cQ0*vv[:,nP]+cQ1*vv[:,nP+1])
    return MCC,MUB    


def genConstraints2DSymm(ddMon,degQ,nQ,cQ1,nEquator,epsilon=0.01):
    print('generating constraints with nEquator = {} and epsilon = {} ...'.format(nEquator,epsilon))
    degPm1=degQ-2;degQm1=degQ-1;nCol=nQ
    vQ=ddMon['vQ']
    vDD11Q=ddMon['vH11Q'];vDD12Q=ddMon['vH12Q'];vDD13Q=ddMon['vH13Q']
    vDD22Q=ddMon['vH22Q'];vDD23Q=ddMon['vH23Q'];vDD33Q=ddMon['vH33Q']
    vCD11Q=ddMon['vCH11Q'];vCD12Q=ddMon['vCH12Q'];vCD13Q=ddMon['vCH13Q']
    vCD22Q=ddMon['vCH22Q'];vCD23Q=ddMon['vCH23Q'];vCD33Q=ddMon['vCH33Q']  
    vUG=genConstraintsPoints2DOpt(nEquator)    
    nPoints=vUG.shape[0]
    nVec=(vUG.shape[1]-3)//3##number of tangent vectors 
    nCons=nPoints*nVec##number of constraints
    GG=np.zeros((nPoints,6*nVec))##'6' is specific to 2D
    vPoints=np.zeros((nPoints,3*nVec))##'3' is specific to 2D
    jj=3;jG=0
    for kk in range(nVec):
        GG[:,jG]=vUG[:,jj]**2;jG+=1
        GG[:,jG]=vUG[:,jj+1]**2;jG+=1
        GG[:,jG]=vUG[:,jj+2]**2;jG+=1
        GG[:,jG]=2*vUG[:,jj]*vUG[:,jj+1];jG+=1
        GG[:,jG]=2*vUG[:,jj]*vUG[:,jj+2];jG+=1
        GG[:,jG]=2*vUG[:,jj+1]*vUG[:,jj+2];jG+=1
        vPoints[:,jj-3:jj]=vUG[:,0:3]
        jj+=3
    GG=GG.reshape((nCons,6))    
    vPoints=vPoints.reshape((nCons,3))
    vPowers1=np.zeros((nCons,nQ))
    vPowers2=np.zeros((nCons,nQ))
    vPowers3=np.zeros((nCons,nQ))
    vv=np.zeros((nCons,nCol))
    nP=0
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vQ[k][0]
        vPowers2[:,k]=vPoints[:,1]**vQ[k][1]
        vPowers3[:,k]=vPoints[:,2]**vQ[k][2]
    vv[:,nP:]=degQm1*vPowers1*vPowers2*vPowers3
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD11Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD11Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD11Q[k][2]
    vv[:,nP:]-=vCD11Q[:]*vPowers1*vPowers2*vPowers3*GG[:,0].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD22Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD22Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD22Q[k][2]
    vv[:,nP:]-=vCD22Q[:]*vPowers1*vPowers2*vPowers3*GG[:,1].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD33Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD33Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD33Q[k][2]
    vv[:,nP:]-=vCD33Q[:]*vPowers1*vPowers2*vPowers3*GG[:,2].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD12Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD12Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD12Q[k][2]
    vv[:,nP:]-=vCD12Q[:]*vPowers1*vPowers2*vPowers3*GG[:,3].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD13Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD13Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD13Q[k][2]
    vv[:,nP:]-=vCD13Q[:]*vPowers1*vPowers2*vPowers3*GG[:,4].reshape((nCons,1))
    for k in range(nQ):
        vPowers1[:,k]=vPoints[:,0]**vDD23Q[k][0]
        vPowers2[:,k]=vPoints[:,1]**vDD23Q[k][1]
        vPowers3[:,k]=vPoints[:,2]**vDD23Q[k][2]
    vv[:,nP:]-=vCD23Q[:]*vPowers1*vPowers2*vPowers3*GG[:,5].reshape((nCons,1))
    MCC=np.zeros((nCons,nCol-2))
    MUB=np.zeros((nCons,1))
    MCC[:,nP:]=vv[:,nP+2:]
    MUB[:,0]=1.0-epsilon-(cQ1*vv[:,nP+1])
    return MCC,MUB        


'''
function:'dataFitSHYqp'
--Calculates the SHY(Q)(P) coefficients by minimizing the weighted distance to the proto-model(Bezier5YS)
--Input: 
----data=the overall data structure of the material
----
'''
def dataFitSHYqp(data,uaxData,lbd,qpSolver='cvxopt',nEquator=200,epsilon=0.01,nSections=19):
    degQ=data['DEG'];degQm1=degQ-1;degQm2=degQ-2;degPm1=degQm2;degPm2=degPm1-1
    nQ,nP=nMonoms(degQ)
    nPm2=nP-2;nP1=nP+1;nP2=nP+2
    sq3=np.sqrt(3.0);sq32=np.sqrt(1.5);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    ddMon=vPoly(degQ)
    vP=ddMon['vP']
    vD1P=ddMon['vD1P'];vD2P=ddMon['vD2P'];vD3P=ddMon['vD3P']
    vC1P=ddMon['vC1P'];vC2P=ddMon['vC2P'];vC3P=ddMon['vC3P']
    vQ=ddMon['vQ']
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']
    cQ0=0.5*(1.0/data['sC']['0']-1.0);cP0=-cQ0
    rt=(1.0-data['rT']['0'])/(1.0+data['rT']['0'])
    rc=(1.0-data['rC']['0'])/(data['sC']['0']*(1.0+data['rC']['0']))
    cP1=0.5*(rt-rc)/sq3;cQ1=0.5*(rt+rc)/sq3
    vYSpoints=protoDataPoints(lbd,data['shapeVBAX'],uaxData,nSections=11,nPointsSegment=10)
    nYSpoints=vYSpoints.shape[0]
    nZZpoints=data['fileData'].shape[0]  ##number of additional data points read from 'fileData'
    nSeg=23
    nsT=(len(data['sT'])-1)*(nSeg-1);nsC=(len(data['sC'])-1)*(nSeg-1)
    nData=2*(nsT+nsC-2)+data['wsTb']*(1+data['wrTb'])+data['wsCb']*(1+data['wrCb'])
    nRow=nData+nYSpoints+nZZpoints
    nCol=nQ+nP
    vv=np.zeros(nCol)##stores monomial values on a (u1,u2,u3) tuple 
    AA=np.zeros((nRow,nCol-4));BB=np.zeros((nRow,1))
    vu=np.zeros((1,3))
    powerP=np.zeros((nP,3));powerP1=np.zeros((nP,3));powerP2=np.zeros((nP,3));powerP3=np.zeros((nP,3))
    powerQ=np.zeros((nQ,3));powerQ1=np.zeros((nQ,3));powerQ2=np.zeros((nQ,3));powerQ3=np.zeros((nQ,3))
    vTT=np.linspace(0.0,1.0,nSeg);vTTr=np.zeros(nSeg)
    vTheta=np.zeros(nsT);sTData=np.zeros(nsT);rTData=np.zeros(nsT)
    lbdST=data['shapeUAX']*uaxData['lambdaSTmax']
    lbdRT=data['shapeUAX']*uaxData['lambdaRTmax']
    hh=0
    for k in range(uaxData['pointsST'].shape[1]-1):
        BS=uaxData['pointsST'][:,k].reshape(2,1)
        BE=uaxData['pointsST'][:,k+1].reshape(2,1)
        TS=uaxData['tanST'][:,k].reshape(2,1)
        TE=uaxData['tanST'][:,k+1].reshape(2,1)
        vData=curveSegF(BS,TS,BE,TE,lbdST,vTT)
        vTheta[hh:hh+nSeg-1]=vData[0,1:]
        sTData[hh:hh+nSeg-1]=vData[1,1:]
        BS=uaxData['pointsRT'][:,k].reshape(2,1)
        BE=uaxData['pointsRT'][:,k+1].reshape(2,1)
        TS=uaxData['tanRT'][:,k].reshape(2,1)
        TE=uaxData['tanRT'][:,k+1].reshape(2,1)
        vTTr[:]=curveSegTheta(vData[0,:],BS[0,0],TS[0,0],BE[0,0],TE[0,0],lbdRT)
        vData=curveSegF(BS,TS,BE,TE,lbdRT,vTTr)
        rTData[hh:hh+nSeg-1]=vData[1,1:]
        hh+=nSeg-1
    vTheta[:]*=(np.pi/180.0)
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    ##print("sTData shape = {}\n".format(sTData.shape),sTData);print("rTData shape = {}\n".format(rTData.shape),rTData)
    vecDiag=np.ones(nRow)
    wData=data['weight']
    wDataZZ=0.0
    wwZZ=0.95
    if(nZZpoints):
        wData=wwZZ*data['weight']
        wDataZZ=(1-wwZZ)*data['weight']/nZZpoints    
    wDataS=0.8*wData/(0.5*nData);wDataR=0.2*wData/(0.5*nData)
    kRow=0
    for kk in range(1,nsT):
        vu[:]=cTheta[kk]**2-0.5*sTheta[kk]**2,0.5*sq3*sTheta[kk]**2,cTheta[kk]*sTheta[kk]*sq3
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/sTData[kk]-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataS
        kRow+=1
        alpha=rTData[kk]+sTheta[kk]**2
        beta=rTData[kk]+cTheta[kk]**2
        omega=cTheta[kk]*sTheta[kk]
        b2a=alpha-0.5*beta
        cc=-b2a*vu[0,0]-0.5*sq3*beta*vu[0,1]+omega*sq3*vu[0,2]
        powerP1[:]=vu**vD1P;powerP2[:]=vu**vD2P;powerP3[:]=vu**vD3P
        vv[0:nP]=cc*degPm1*vv[0:nP]+b2a*vC1P[:]*powerP1[:,0]*powerP1[:,1]*powerP1[:,2]+\
                 0.5*beta*sq3*vC2P[:]*powerP2[:,0]*powerP2[:,1]*powerP2[:,2]-\
                 omega*sq3*vC3P[:]*powerP3[:,0]*powerP3[:,1]*powerP3[:,2]
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q;powerQ3[:]=vu**vD3Q
        vv[nP:]=cc*degQm1*vv[nP:]+b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
                 0.5*beta*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]-\
                 omega*sq3*vC3Q[:]*powerQ3[:,0]*powerQ3[:,1]*powerQ3[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=cc-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataR
        kRow+=1
    vTheta=np.zeros(nsC);sCData=np.zeros(nsC);rCData=np.zeros(nsC)
    lbdSC=data['shapeUAX']*uaxData['lambdaSCmax']
    lbdRC=data['shapeUAX']*uaxData['lambdaRCmax']
    hh=0
    for k in range(uaxData['pointsSC'].shape[1]-1):
        BS=uaxData['pointsSC'][:,k].reshape(2,1)
        BE=uaxData['pointsSC'][:,k+1].reshape(2,1)
        TS=uaxData['tanSC'][:,k].reshape(2,1)
        TE=uaxData['tanSC'][:,k+1].reshape(2,1)
        vData=curveSegF(BS,TS,BE,TE,lbdSC,vTT)
        vTheta[hh:hh+nSeg-1]=vData[0,1:]
        sCData[hh:hh+nSeg-1]=vData[1,1:]
        BS=uaxData['pointsRC'][:,k].reshape(2,1)
        BE=uaxData['pointsRC'][:,k+1].reshape(2,1)
        TS=uaxData['tanRC'][:,k].reshape(2,1)
        TE=uaxData['tanRC'][:,k+1].reshape(2,1)
        vTTr[:]=curveSegTheta(vData[0,:],BS[0,0],TS[0,0],BE[0,0],TE[0,0],lbdRC)
        vData=curveSegF(BS,TS,BE,TE,lbdRC,vTTr)
        rCData[hh:hh+nSeg-1]=vData[1,1:]
        hh+=nSeg-1
    vTheta[:]*=(np.pi/180.0)
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    for kk in range(1,nsC):
        vu[:]=-cTheta[kk]**2+0.5*sTheta[kk]**2,-0.5*sq3*sTheta[kk]**2,-cTheta[kk]*sTheta[kk]*sq3
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/sCData[kk]-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataS
        kRow+=1
        alpha=rCData[kk]+sTheta[kk]**2
        beta=rCData[kk]+cTheta[kk]**2
        omega=cTheta[kk]*sTheta[kk]
        b2a=alpha-0.5*beta
        cc=-b2a*vu[0,0]-0.5*sq3*beta*vu[0,1]+omega*sq3*vu[0,2]
        powerP1[:]=vu**vD1P;powerP2[:]=vu**vD2P;powerP3[:]=vu**vD3P
        vv[0:nP]=cc*degPm1*vv[0:nP]+b2a*vC1P[:]*powerP1[:,0]*powerP1[:,1]*powerP1[:,2]+\
                 0.5*beta*sq3*vC2P[:]*powerP2[:,0]*powerP2[:,1]*powerP2[:,2]-\
                 omega*sq3*vC3P[:]*powerP3[:,0]*powerP3[:,1]*powerP3[:,2]
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q;powerQ3[:]=vu**vD3Q
        vv[nP:]=cc*degQm1*vv[nP:]+b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
                 0.5*beta*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]-\
                 omega*sq3*vC3Q[:]*powerQ3[:,0]*powerQ3[:,1]*powerQ3[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=cc-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataR
        kRow+=1    
    if(data['wsTb']):
        vu[:]=0.5,0.5*sq3,0.0   ##sTb
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/data['sTb']-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataS
        kRow+=1
    if(data['wsTb']*data['wrTb']):  ##rTb
        b2a=data['rTb']+0.5
        cc=0.5*(data['rTb']-1.0)
        powerP1[:]=vu**vD1P;powerP2[:]=vu**vD2P
        vv[0:nP]=cc*degPm1*vv[0:nP]-b2a*vC1P[:]*powerP1[:,0]*powerP1[:,1]*powerP1[:,2]+\
                 0.5*sq3*vC2P[:]*powerP2[:,0]*powerP2[:,1]*powerP2[:,2]
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q
        vv[nP:]=cc*degQm1*vv[nP:]-b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
                 0.5*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=cc-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataR
        kRow+=1
    if(data['wsCb']):
        vu[:]=-0.5,-0.5*sq3,0.0   ##sCb
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/data['sCb']-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataS
        kRow+=1
    if(data['wsCb']*data['wrCb']):  ##rCb
        b2a=data['rCb']+0.5
        cc=0.5*(-data['rCb']+1.0)
        powerP1[:]=vu**vD1P;powerP2[:]=vu**vD2P
        vv[0:nP]=cc*degPm1*vv[0:nP]-b2a*vC1P[:]*powerP1[:,0]*powerP1[:,1]*powerP1[:,2]+\
                 0.5*sq3*vC2P[:]*powerP2[:,0]*powerP2[:,1]*powerP2[:,2]
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q
        vv[nP:]=cc*degQm1*vv[nP:]-b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
                 0.5*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=cc-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataR
        kRow+=1
    wDataYS=(1.0-data['weight'])/nYSpoints    
    for kk in range(nYSpoints):
        vu[0,:]=(2*vYSpoints[kk,0]-vYSpoints[kk,1])/sq6,vYSpoints[kk,1]/sq2,vYSpoints[kk,2]*sq2
        modS=np.sqrt(vu[0,0]**2+vu[0,1]**2+vu[0,2]**2)
        vu[0,:]/=modS
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/(sq32*modS)-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataYS
        kRow+=1 
    vYSpoints=data['fileData']    
    for kk in range(nZZpoints):
        vu[0,:]=(2*vYSpoints[kk,0]-vYSpoints[kk,1])/sq6,vYSpoints[kk,1]/sq2,vYSpoints[kk,2]*sq2
        modS=np.sqrt(vu[0,0]**2+vu[0,1]**2+vu[0,2]**2)
        vu[0,:]/=modS
        powerP[:]=vu**vP
        vv[0:nP]=powerP[:,0]*powerP[:,1]*powerP[:,2]
        powerQ[:]=vu**vQ
        vv[nP:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,0:nPm2]=vv[2:nP];AA[kRow,nPm2:]=vv[nP2:]
        BB[kRow,0]=1.0/(sq32*modS)-1.0-(cP0*vv[0]+cP1*vv[1]+cQ0*vv[nP]+cQ1*vv[nP1])
        vecDiag[kRow]=wDataZZ
        kRow+=1            
    vDiag=np.diag(vecDiag)
    MAA=np.dot(AA.T,np.dot(vDiag,AA));MAA=0.5*(MAA.T+MAA)+1.0e-12*np.eye(MAA.shape[0])
    MBB=np.dot(AA.T,np.dot(vDiag,BB))
    print("Generating constraints....")
    MCC,MUB=genConstraints2D(ddMon,degQ,nP,nQ,cP0,cP1,cQ0,cQ1,nEquator,epsilon)
    print("{}: Calculating SHYqp parameters....".format(qpSolver))    
    #meq=0
    #vsq=qpg.solve_qp(MAA,MBB.reshape((MBB.shape[0],)),-MCC.T,-MUB.reshape((MUB.shape[0],)), meq)[0]
    if(qpSolver=='quadprog'):##use quadprog
        meq=0    
        vsq=qpg.solve_qp(MAA,MBB.reshape((MBB.shape[0],)),-MCC.T,-MUB.reshape((MUB.shape[0],)),meq)[0]
    elif(qpSolver=='cvxopt'):##use cvxopt    
        args = [cvxopt.matrix(MAA), cvxopt.matrix(-MBB.reshape((MBB.shape[0],))),cvxopt.matrix(MCC),cvxopt.matrix(MUB.reshape((MUB.shape[0],)))]
        vsq=cvxopt.solvers.qp(*args)['x']
        vsq=np.array(vsq).reshape(len(vsq))
    else:##unknown solver
        print('sqSolver = {}: unknown solver\nCalculations aborted'.format(sqSolver))
        exit()       
    vCoeff=np.zeros(nCol)
    vCoeff[0]=cP0;vCoeff[1]=cP1;vCoeff[2:nP]=vsq[0:nP-2];vCoeff[nP]=cQ0;vCoeff[nP+1]=cQ1;vCoeff[nP+2:]=vsq[nP-2:]
    return vCoeff,ddMon,nQ,nP


'''
function:'dataFitSHYqpSymm'
--Calculates the SHY(Q) coefficients by minimizing the weighted distance to the proto-model(Bezier5YS)
'''
def dataFitSHYqpSymm(data,uaxData,lbd,qpSolver='cvxopt',nEquator=200,epsilon=0.01,nSections=15):
    degQ=data['DEG'];degQm1=degQ-1;degQm2=degQ-2;degPm1=degQm2;degPm2=degPm1-1
    nQ,nP=nMonoms(degQ)
    nPm2=nP-2;nP1=nP+1;nP2=nP+2
    sq3=np.sqrt(3.0);sq32=np.sqrt(1.5);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    ddMon=vPoly(degQ)
    vQ=ddMon['vQ']
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']
    cQ0=0.0
    rt=(1.0-data['rT']['0'])/(1.0+data['rT']['0'])
    cQ1=rt/sq3
    vYSpoints=protoDataPoints(lbd,data['shapeVBAX'],uaxData,nSections,nPointsSegment=10)
    nYSpoints=vYSpoints.shape[0]
    ###nsT=len(data['sT']);nsC=len(data['sC'])
    nSeg=21
    nsT=(len(data['sT'])-1)*(nSeg-1)
    nData=2*(nsT-1)+data['wsTb']*(1+data['wrTb'])
    nZZpoints=data['fileData'].shape[0]  ##number of additional data points read from 'fileData'
    nRow=nData+nYSpoints+nZZpoints
    nCol=nQ
    vv=np.zeros(nCol)##stores monomial values on a (u1,u2,u3) tuple 
    AA=np.zeros((nRow,nCol-2));BB=np.zeros((nRow,1))
    vu=np.zeros((1,3))
    powerQ=np.zeros((nQ,3));powerQ1=np.zeros((nQ,3));powerQ2=np.zeros((nQ,3));powerQ3=np.zeros((nQ,3))
    vTT=np.linspace(0.0,1.0,nSeg);vTTr=np.zeros(nSeg)
    vTheta=np.zeros(nsT);sTData=np.zeros(nsT);rTData=np.zeros(nsT)
    lbdST=data['shapeUAX']*uaxData['lambdaSTmax']
    lbdRT=data['shapeUAX']*uaxData['lambdaRTmax']
    hh=0
    for k in range(uaxData['pointsST'].shape[1]-1):
        BS=uaxData['pointsST'][:,k].reshape(2,1)
        BE=uaxData['pointsST'][:,k+1].reshape(2,1)
        TS=uaxData['tanST'][:,k].reshape(2,1)
        TE=uaxData['tanST'][:,k+1].reshape(2,1)
        vData=curveSegF(BS,TS,BE,TE,lbdST,vTT)
        vTheta[hh:hh+nSeg-1]=vData[0,1:]
        sTData[hh:hh+nSeg-1]=vData[1,1:]
        BS=uaxData['pointsRT'][:,k].reshape(2,1)
        BE=uaxData['pointsRT'][:,k+1].reshape(2,1)
        TS=uaxData['tanRT'][:,k].reshape(2,1)
        TE=uaxData['tanRT'][:,k+1].reshape(2,1)
        vTTr[:]=curveSegTheta(vData[0,:],BS[0,0],TS[0,0],BE[0,0],TE[0,0],lbdRT)
        vData=curveSegF(BS,TS,BE,TE,lbdRT,vTTr)
        rTData[hh:hh+nSeg-1]=vData[1,1:]
        hh+=nSeg-1
    vTheta[:]*=(np.pi/180.0)
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    ##print("sTData shape = {}\n".format(sTData.shape),sTData);print("rTData shape = {}\n".format(rTData.shape),rTData)
    vecDiag=np.ones(nRow)
    wData=data['weight']
    wDataZZ=0.0
    wwZZ=0.95
    if(nZZpoints):
        wData=wwZZ*data['weight']
        wDataZZ=(1-wwZZ)*data['weight']/nZZpoints    
    wDataS=0.8*wData/(0.5*nData);wDataR=0.2*wData/(0.5*nData)
    kRow=0
    #vecDiag=np.ones(nRow)
    #wData=data['weight']/nData;
    #wDataS=0.8*data['weight']/(0.5*nData);wDataR=0.2*data['weight']/(0.5*nData)
    #kRow=0
    for kk in range(1,nsT):
        vu[:]=cTheta[kk]**2-0.5*sTheta[kk]**2,0.5*sq3*sTheta[kk]**2,cTheta[kk]*sTheta[kk]*sq3
        powerQ[:]=vu**vQ
        vv[:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=1.0/sTData[kk]-1.0-(cQ1*vv[1])
        vecDiag[kRow]=wDataS
        kRow+=1
        alpha=rTData[kk]+sTheta[kk]**2
        beta=rTData[kk]+cTheta[kk]**2
        omega=cTheta[kk]*sTheta[kk]
        b2a=alpha-0.5*beta
        cc=-b2a*vu[0,0]-0.5*sq3*beta*vu[0,1]+omega*sq3*vu[0,2]
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q;powerQ3[:]=vu**vD3Q
        vv[:]=cc*degQm1*vv[:]+b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
                 0.5*beta*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]-\
                 omega*sq3*vC3Q[:]*powerQ3[:,0]*powerQ3[:,1]*powerQ3[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=cc-(cQ1*vv[1])
        vecDiag[kRow]=wDataR
        kRow+=1
    if(data['wsTb']):
        vu[:]=0.5,0.5*sq3,0.0   ##sTb
        powerQ[:]=vu**vQ
        vv[:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=1.0/data['sTb']-1.0-(cQ1*vv[1])
        vecDiag[kRow]=wDataS
        kRow+=1
    if(data['wsTb']*data['wrTb']):  ##rTb
        b2a=data['rTb']+0.5
        cc=0.5*(data['rTb']-1.0)
        powerQ1[:]=vu**vD1Q;powerQ2[:]=vu**vD2Q
        vv[:]=cc*degQm1*vv[:]-b2a*vC1Q[:]*powerQ1[:,0]*powerQ1[:,1]*powerQ1[:,2]+\
              0.5*sq3*vC2Q[:]*powerQ2[:,0]*powerQ2[:,1]*powerQ2[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=cc-(cQ1*vv[1])
        vecDiag[kRow]=wDataR
        kRow+=1
    wDataYS=(1.0-data['weight'])/nYSpoints    
    for kk in range(nYSpoints):
        vu[:]=(2*vYSpoints[kk,0]-vYSpoints[kk,1])/sq6,vYSpoints[kk,1]/sq2,vYSpoints[kk,2]*sq2
        modS=np.sqrt(vu[0,0]**2+vu[0,1]**2+vu[0,2]**2)
        vu[:]/=modS
        powerQ[:]=vu**vQ
        vv[:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=1.0/(sq32*modS)-1.0-(cQ1*vv[1])
        vecDiag[kRow]=wDataYS
        kRow+=1
    vYSpoints=data['fileData']    
    for kk in range(nZZpoints):
        vu[:]=(2*vYSpoints[kk,0]-vYSpoints[kk,1])/sq6,vYSpoints[kk,1]/sq2,vYSpoints[kk,2]*sq2
        modS=np.sqrt(vu[0,0]**2+vu[0,1]**2+vu[0,2]**2)
        vu[:]/=modS
        powerQ[:]=vu**vQ
        vv[:]=powerQ[:,0]*powerQ[:,1]*powerQ[:,2]
        AA[kRow,:]=vv[2:]
        BB[kRow,0]=1.0/(sq32*modS)-1.0-(cQ1*vv[1])
        vecDiag[kRow]=wDataZZ
        kRow+=1          
    vDiag=np.diag(vecDiag)
    MAA=np.dot(AA.T,np.dot(vDiag,AA));MAA=0.5*(MAA.T+MAA)+1.0e-12*np.eye(MAA.shape[0])
    MBB=np.dot(AA.T,np.dot(vDiag,BB))
    print("Generating constraints....")
    MCC,MUB=genConstraints2DSymm(ddMon,degQ,nQ,cQ1,nEquator,epsilon)
    print("{}: Calculating SHYq parameters....".format(qpSolver))    
    if(qpSolver=='quadprog'):##use quadprog
        meq=0    
        vsq=qpg.solve_qp(MAA,MBB.reshape((MBB.shape[0],)),-MCC.T,-MUB.reshape((MUB.shape[0],)),meq)[0]
    elif(qpSolver=='cvxopt'):##use cvxopt    
        args = [cvxopt.matrix(MAA), cvxopt.matrix(-MBB.reshape((MBB.shape[0],))),cvxopt.matrix(MCC),cvxopt.matrix(MUB.reshape((MUB.shape[0],)))]
        vsq=cvxopt.solvers.qp(*args)['x']
        vsq=np.array(vsq).reshape(len(vsq))
    else:##unknown solver
        print('sqSolver = {}: unknown solver\nCalculations aborted'.format(sqSolver))
        exit()        
    vCoeff=np.zeros(nP+nCol)
    vCoeff[nP]=cQ0;vCoeff[nP+1]=cQ1;vCoeff[nP+2:]=vsq[:]
    return vCoeff,ddMon,nQ,nP


def SHYqp_uax_Plot(uaxData,vCoeff,ddMon,nQ,nP,savePng=False):
    degQ=ddMon['nQ']
    degQm1=degQ-1;degPm1=degQ-2
    sq3=np.sqrt(3.0);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)
    vD1P=ddMon['vD1P'];vD2P=ddMon['vD2P'];vD3P=ddMon['vD3P']
    vC1P=ddMon['vC1P'];vC2P=ddMon['vC2P'];vC3P=ddMon['vC3P']
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']
    vTheta=np.linspace(0,np.pi/2,100)
    nTheta=vTheta.shape[0]
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    vu=np.zeros((nTheta,3))
    vu[:,0]=cTheta**2-0.5*sTheta**2
    vu[:,1]=0.5*sq3*sTheta**2
    vu[:,2]=sq3*cTheta*sTheta
    vu3=vu[:,2]**2
    vMon=np.zeros((nTheta,nP))
    vMonD1=np.zeros((nTheta,nP));vMonD2=np.zeros((nTheta,nP))
    vTerms=np.zeros((nTheta,nPidx))
    vTermsD1=np.zeros((nTheta,nPidx));vTermsD2=np.zeros((nTheta,nPidx));vTermsD3=np.zeros((nTheta,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]    
    vvP=vTerms[:,-1];vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    for k in range(nPidx-2,-1,-1):
        vvP[:]=vTerms[:,k]+vu3*vvP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2 
    for k in range(nPidx-2,0,-1):        
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3
    vPD3[:]*=vu[:,2]    
    vMon=np.zeros((nTheta,nQ))
    vMonD1=np.zeros((nTheta,nQ));vMonD2=np.zeros((nTheta,nQ))
    vTerms=np.zeros((nTheta,nQidx))
    vTermsD1=np.zeros((nTheta,nQidx));vTermsD2=np.zeros((nTheta,nQidx));vTermsD3=np.zeros((nTheta,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]        
    vvQ=vTerms[:,-1];vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    for k in range(nQidx-2,-1,-1):
        vvQ[:]=vTerms[:,k]+vu3*vvQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
    for k in range(nQidx-2,0,-1):        
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3
    vQD3[:]*=vu[:,2]      
    sigThetaT=1.0/(1.0+vvP+vvQ)
    vPQ=1.0-degPm1*vvP-degQm1*vvQ
    D1=vu[:,0]*vPQ+vPD1+vQD1
    D2=vu[:,1]*vPQ+vPD2+vQD2
    DZ=(vu[:,2]*vPQ+vPD3+vQD3)/sq2
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rThetaT=(2*DZ*sTheta*cTheta-DX*sTheta**2-DY*cTheta**2)/(DX+DY)
    vu[:,0]=-vu[:,0]
    vu[:,1]=-vu[:,1]
    vu[:,2]=-vu[:,2]
    vMon=np.zeros((nTheta,nP))
    vMonD1=np.zeros((nTheta,nP));vMonD2=np.zeros((nTheta,nP))
    vTerms=np.zeros((nTheta,nPidx))
    vTermsD1=np.zeros((nTheta,nPidx));vTermsD2=np.zeros((nTheta,nPidx));vTermsD3=np.zeros((nTheta,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]    
    vvP=vTerms[:,-1];vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    for k in range(nPidx-2,-1,-1):
        vvP[:]=vTerms[:,k]+vu3*vvP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2 
    for k in range(nPidx-2,0,-1):        
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3
    vPD3[:]*=vu[:,2]    
    vMon=np.zeros((nTheta,nQ))
    vMonD1=np.zeros((nTheta,nQ));vMonD2=np.zeros((nTheta,nQ))
    vTerms=np.zeros((nTheta,nQidx))
    vTermsD1=np.zeros((nTheta,nQidx));vTermsD2=np.zeros((nTheta,nQidx));vTermsD3=np.zeros((nTheta,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]        
    vvQ=vTerms[:,-1];vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    for k in range(nQidx-2,-1,-1):
        vvQ[:]=vTerms[:,k]+vu3*vvQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
    for k in range(nQidx-2,0,-1):        
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3
    vQD3[:]*=vu[:,2]      
    sigThetaC=1.0/(1.0+vvP+vvQ)
    vPQ=1.0-degPm1*vvP-degQm1*vvQ
    D1=vu[:,0]*vPQ+vPD1+vQD1
    D2=vu[:,1]*vPQ+vPD2+vQD2
    DZ=(vu[:,2]*vPQ+vPD3+vQD3)/sq2
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rThetaC=(2*DZ*sTheta*cTheta-DX*sTheta**2-DY*cTheta**2)/(DX+DY)
    vTheta*=180.0/np.pi
    fgS=plt.figure();axS=fgS.add_subplot(1,1,1)
    axS.plot(vTheta,sigThetaT,'k-',linewidth=lineWidthUax)
    axS.plot(uaxData['pointsST'][0,:],uaxData['pointsST'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    axS.plot(vTheta,sigThetaC,'k--',linewidth=lineWidthUax)
    axS.plot(uaxData['pointsSC'][0,:],uaxData['pointsSC'][1,:],'bo',markerfacecolor='b',markersize=6.5,label='Exp data')
    ##axS.text(0,0.95*max(uaxData['pointsST'][1,:]),r'$\overline{\sigma}_{\theta}$',fontsize=14) 
    axS.text(0.022,0.93,r'$\overline{\sigma}_{\theta}$',fontsize=16,
    horizontalalignment='center',verticalalignment='center',transform = axS.transAxes)
    ##axS.set_ylim([0.9*np.min(uaxData['pointsST'][1,:]),1.05*np.max(uaxData['pointsST'][1,:])])  
    maxS=np.max([np.max(uaxData['pointsST'][1,:]),np.max(uaxData['pointsSC'][1,:]),np.max(sigThetaT),np.max(sigThetaC)])
    minS=np.min([np.min(uaxData['pointsST'][1,:]),np.min(uaxData['pointsSC'][1,:]),np.min(sigThetaT),np.min(sigThetaC)])
    axS.set_ylim([0.97*minS,1.03*maxS])   
    axS.set_xticks([0,15,30,45,60,75,90])
    axS.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=12) 
    axS.tick_params(axis="y", labelsize=12)    
    axS.grid()
    fgR=plt.figure();axR=fgR.add_subplot(1,1,1)
    axR.plot(vTheta,rThetaT,'k-',linewidth=lineWidthUax)
    axR.plot(uaxData['pointsRT'][0,:],uaxData['pointsRT'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    axR.plot(vTheta,rThetaC,'k--',linewidth=lineWidthUax)
    axR.plot(uaxData['pointsRC'][0,:],uaxData['pointsRC'][1,:],'bo',markerfacecolor='b',markersize=6.5,label='Exp data') 
    ##axR.text(0,0.95*max(uaxData['pointsRT'][1,:]),r'$r_{\theta}$',fontsize=14) 
    ##axR.set_ylim([0.9*np.min(uaxData['pointsRT'][1,:]),1.05*np.max(uaxData['pointsRT'][1,:])])     
    axR.text(0.022,0.92,r'$r_{\theta}$',fontsize=16,
    horizontalalignment='center',verticalalignment='center',transform = axR.transAxes)
    maxR=np.max([np.max(uaxData['pointsRT'][1,:]),np.max(uaxData['pointsRC'][1,:]),np.max(rThetaT),np.max(rThetaC)])
    minR=np.min([np.min(uaxData['pointsRT'][1,:]),np.min(uaxData['pointsRC'][1,:]),np.min(rThetaT),np.min(rThetaC)])
    if(minR<0.5):
        axR.set_ylim([-0.1,1.05*maxR])
    else:
        axR.set_ylim([0.9*minR,1.05*maxR])    
    axR.set_xticks([0,15,30,45,60,75,90])
    axR.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=12)
    axR.tick_params(axis="y", labelsize=12)     
    axR.grid()
    if(savePng):
        fgS.savefig('{name}.{ext}'.format(name=figDir+uaxData['name']+'_SHYqp_deg'+str(degQ)+'_UaxS',ext='png'),dpi=300,bbox_inches='tight') 
        fgR.savefig('{name}.{ext}'.format(name=figDir+uaxData['name']+'_SHYqp_deg'+str(degQ)+'_UaxR',ext='png'),dpi=300,bbox_inches='tight')    
    return


def SHYqp_bax_Plot(vCoeff,ddMon,nQ,nP,assym,name,savePng=False):
    sq3=np.sqrt(3.0);sq32=np.sqrt(1.5);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)##;print("vPidx\n",vPidx)
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)##;print("vQidx\n",vQidx)
    nx=400
    x=np.linspace(-1.5,1.5,nx)
    y=np.linspace(-1.5,1.5,nx)
    nSamples=nx*nx    
    vx,vy=np.meshgrid(x,y)
    vx=vx.reshape(nSamples)
    vy=vy.reshape(nSamples)
    zz=np.zeros(nSamples)
    vx=(2.0*vx-vy)/sq6;vy=vy/sq2;vx2y2=vx**2+vy**2
    if(assym):
        vsxy=np.linspace(0,1/(sq3*(1.0+vCoeff[-1])),6)  ##;print("c=",vCoeff[-1])
    else:  
        vsxy=np.linspace(0,0.965/(sq3*(1.0+vCoeff[-1])),7) 
        ##vsxy=np.linspace(0,0.9999/(sq3*(1.0+vCoeff[-1])),11)        
    vMod=np.zeros(nSamples)
    vux=np.zeros(nSamples);vuy=np.zeros(nSamples);vuz=np.zeros(nSamples)
    vMonP=np.zeros((nSamples,nP));vTermsP=np.zeros((nSamples,nPidx))
    vMonQ=np.zeros((nSamples,nQ));vTermsQ=np.zeros((nSamples,nQidx))
    print("Calculating biaxial sections for plots...")
    fg=plt.figure();ax=fg.add_subplot(1,1,1)
    for sxy in vsxy:
        print("--section sigma_xy = {}".format(sxy))
        vMod[:]=np.sqrt(vx2y2+2*sxy**2)
        vux[:]=vx/vMod;vuy[:]=vy/vMod;vuz[:]=(2*sxy**2)/vMod**2
        for k in range(nP):
            vMonP[:,k]=vux**vP[k][0]*vuy**vP[k][1]
        for k in range(nPidx):
            vTermsP[:,k]=np.dot(vMonP[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vsP=vTermsP[:,-1]
        for k in range(nPidx-2,-1,-1):
            ##print("vsP[{}]".format(k),vsP)
            vsP[:]=vTermsP[:,k]+vuz*vsP
        for k in range(nQ):
            vMonQ[:,k]=vux**vQ[k][0]*vuy**vQ[k][1]
        for k in range(nQidx):
            vTermsQ[:,k]=np.dot(vMonQ[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        ##vTermsQ[:,nQidx-1]=vCoeff[-1]        
        ##print("check vCoeff:",vCoeff[nP+vQidx[2][1][0]:nP+vQidx[2][1][-1]]);print(nP+vQidx[2][1][0],nP+vQidx[2][1][-1]);print(vMonQ[:,-1]) ;exit()           
        vsQ=vTermsQ[:,-1]
        for k in range(nQidx-2,-1,-1):
            ##print("vsQ[{}]".format(k),vsQ)
            vsQ[:]=vTermsQ[:,k]+vuz*vsQ
        zz[:]=sq32*vMod*(1.0+vsP+vsQ)
        ax.contour(x,y,zz.reshape(nx,nx),levels=[1])
    ax.set_aspect('equal')
    ax.grid()
    ax.text(0.06,0.925,r'$\overline{\sigma}_{yy}$',fontsize=16,
    horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    ax.text(0.94,0.05,r'$\overline{\sigma}_{xx}$',fontsize=16,
    horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    if(savePng):
        fg.savefig('{name}.{ext}'.format(name=figDir+name+'_SHYqp_deg'+str(ddMon['nQ'])+'_Biax',ext='png'),dpi=300,bbox_inches='tight')
    return

def SHYqp_Predictions(uaxData,vCoeff,ddMon,nQ,nP,qpSolver,cvxCheck,fileCoeff=None):
    degQ=ddMon['nQ']
    degQm1=degQ-1;degPm1=degQ-2
    sq3=np.sqrt(3.0);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)
    vD1P=ddMon['vD1P'];vD2P=ddMon['vD2P'];vD3P=ddMon['vD3P']
    vC1P=ddMon['vC1P'];vC2P=ddMon['vC2P'];vC3P=ddMon['vC3P']
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']
    if(fileCoeff):
        try:
            ff=open(fileCoeff,'r')
        except IOError as err:
            print(err)
            print('SHYqp_Predictions: early exit (with no report)')
        vCoeff=np.zeros(nP+nQ)
        k=0        
        for line in ff:
            vCoeff[k]=line
            k+=1
        ff.close()            
    vTheta=(np.pi/180.0)*np.array(uaxData['pointsST'][0,:])
    nTheta=vTheta.shape[0]
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    vu=np.zeros((nTheta,3))
    vu[:,0]=cTheta**2-0.5*sTheta**2
    vu[:,1]=0.5*sq3*sTheta**2
    vu[:,2]=sq3*cTheta*sTheta
    vu3=vu[:,2]**2
    vMon=np.zeros((nTheta,nP))
    vMonD1=np.zeros((nTheta,nP));vMonD2=np.zeros((nTheta,nP))
    vTerms=np.zeros((nTheta,nPidx))
    vTermsD1=np.zeros((nTheta,nPidx));vTermsD2=np.zeros((nTheta,nPidx));vTermsD3=np.zeros((nTheta,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]    
    vvP=vTerms[:,-1];vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    for k in range(nPidx-2,-1,-1):
        vvP[:]=vTerms[:,k]+vu3*vvP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2 
    for k in range(nPidx-2,0,-1):        
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3
    vPD3[:]*=vu[:,2]    
    vMon=np.zeros((nTheta,nQ))
    vMonD1=np.zeros((nTheta,nQ));vMonD2=np.zeros((nTheta,nQ))
    vTerms=np.zeros((nTheta,nQidx))
    vTermsD1=np.zeros((nTheta,nQidx));vTermsD2=np.zeros((nTheta,nQidx));vTermsD3=np.zeros((nTheta,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]        
    vvQ=vTerms[:,-1];vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    for k in range(nQidx-2,-1,-1):
        vvQ[:]=vTerms[:,k]+vu3*vvQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
    for k in range(nQidx-2,0,-1):        
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3
    vQD3[:]*=vu[:,2]      
    sigThetaT=1.0/(1.0+vvP+vvQ)
    vPQ=1.0-degPm1*vvP-degQm1*vvQ
    D1=vu[:,0]*vPQ+vPD1+vQD1
    D2=vu[:,1]*vPQ+vPD2+vQD2
    DZ=(vu[:,2]*vPQ+vPD3+vQD3)/sq2
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rThetaT=(2*DZ*sTheta*cTheta-DX*sTheta**2-DY*cTheta**2)/(DX+DY)
    vTheta=(np.pi/180.0)*np.array(uaxData['pointsSC'][0,:])
    nTheta=vTheta.shape[0]
    cTheta=np.cos(vTheta);sTheta=np.sin(vTheta)
    vu=np.zeros((nTheta,3))
    vu[:,0]=-cTheta**2+0.5*sTheta**2
    vu[:,1]=-0.5*sq3*sTheta**2
    vu[:,2]=-sq3*cTheta*sTheta
    vu3=vu[:,2]**2
    vMon=np.zeros((nTheta,nP))
    vMonD1=np.zeros((nTheta,nP));vMonD2=np.zeros((nTheta,nP))
    vTerms=np.zeros((nTheta,nPidx))
    vTermsD1=np.zeros((nTheta,nPidx));vTermsD2=np.zeros((nTheta,nPidx));vTermsD3=np.zeros((nTheta,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]    
    vvP=vTerms[:,-1];vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    for k in range(nPidx-2,-1,-1):
        vvP[:]=vTerms[:,k]+vu3*vvP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2 
    for k in range(nPidx-2,0,-1):        
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3
    vPD3[:]*=vu[:,2]    
    vMon=np.zeros((nTheta,nQ))
    vMonD1=np.zeros((nTheta,nQ));vMonD2=np.zeros((nTheta,nQ))
    vTerms=np.zeros((nTheta,nQidx))
    vTermsD1=np.zeros((nTheta,nQidx));vTermsD2=np.zeros((nTheta,nQidx));vTermsD3=np.zeros((nTheta,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]        
    vvQ=vTerms[:,-1];vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    for k in range(nQidx-2,-1,-1):
        vvQ[:]=vTerms[:,k]+vu3*vvQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
    for k in range(nQidx-2,0,-1):        
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3
    vQD3[:]*=vu[:,2]      
    sigThetaC=1.0/(1.0+vvP+vvQ)
    vPQ=1.0-degPm1*vvP-degQm1*vvQ
    D1=vu[:,0]*vPQ+vPD1+vQD1
    D2=vu[:,1]*vPQ+vPD2+vQD2
    DZ=(vu[:,2]*vPQ+vPD3+vQD3)/sq2
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rThetaC=(2*DZ*sTheta*cTheta-DX*sTheta**2-DY*cTheta**2)/(DX+DY)
    vv=np.array([0.5,0.5*sq3,0.0])
    PP=0.0;k=0
    for m in vP[0:degQ]:
        PP+=vCoeff[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    QQ=0.0;k=nP
    for m in vQ[0:degQ+1]:
        QQ+=vCoeff[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1   
    sTb=1.0/(1.0+PP+QQ)
    vPQ=1.0-degPm1*PP-degQm1*QQ
    D1P=0.0;k=0
    for m in vD1P[0:degQ]:
        D1P+=vCoeff[k]*vC1P[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    D1Q=0.0;k=0
    for m in vD1Q[0:degQ+1]:
        D1Q+=vCoeff[nP+k]*vC1Q[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1 
    D2P=0.0;k=0
    for m in vD2P[0:degQ]:
        D2P+=vCoeff[k]*vC2P[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    D2Q=0.0;k=0
    for m in vD2Q[0:degQ+1]:
        D2Q+=vCoeff[nP+k]*vC2Q[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1     
    D1=vv[0]*vPQ+D1P+D1Q
    D2=vv[1]*vPQ+D2P+D2Q
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rTb=DY/DX    
    ##print('predicted balanced-biaxial sTb = ',sTb)
    vv[:]=-vv
    PP=0.0;k=0
    for m in vP[0:degQ]:
        PP+=vCoeff[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    QQ=0.0;k=nP
    for m in vQ[0:degQ+1]:
        QQ+=vCoeff[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1   
    sCb=1.0/(1.0+PP+QQ)        
    vPQ=1.0-degPm1*PP-degQm1*QQ
    D1P=0.0;k=0
    for m in vD1P[0:degQ]:
        D1P+=vCoeff[k]*vC1P[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    D1Q=0.0;k=0
    for m in vD1Q[0:degQ+1]:
        D1Q+=vCoeff[nP+k]*vC1Q[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1 
    D2P=0.0;k=0
    for m in vD2P[0:degQ]:
        D2P+=vCoeff[k]*vC2P[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1
    D2Q=0.0;k=0
    for m in vD2Q[0:degQ+1]:
        D2Q+=vCoeff[nP+k]*vC2Q[k]*(vv[0]**m[0])*(vv[1]**m[1])
        k+=1     
    D1=vv[0]*vPQ+D1P+D1Q
    D2=vv[1]*vPQ+D2P+D2Q
    DX=np.sqrt(2.0/3.0)*D1
    DY=D2/sq2-D1/sq6
    rCb=DY/DX    
    ##print('predicted balanced-biaxial sCb = ',sCb)
    ff=open(figDir+uaxData['name']+'_SHYqp_deg'+str(degQ)+'_Err_and_Coeff.txt','w')
    ff.write(uaxData['name']+'\nSHYqp degree: '+str(degQ)+'\n')
    ff.write('Solver: '+qpSolver+'\n')
    ff.write('---Convexity check with Hessian leading principal minors: min(det_1), min(det_2), min(det_3)\n')
    ff.write('{}, {}, {}\n'.format(cvxCheck[0],cvxCheck[1],cvxCheck[2]))
    ff.write('---Convexity check with Gaussian curvature: minimum found  = {}\n'.format(cvxCheck[3]))
    ff.write('---All stresses square root error: {:.4f}\n'.format(
    np.sqrt((sTb-uaxData['sTb'])**2+(sCb-uaxData['sCb'])**2+np.sum((sigThetaT-uaxData['pointsST'][1,:])**2)+np.sum((sigThetaC-uaxData['pointsSC'][1,:])**2))))
    ff.write('---All r-Values square root error: {:.4f}\n'.format(
    np.sqrt((rTb-uaxData['rTb'])**2+(rCb-uaxData['rCb'])**2+np.sum((rThetaT-uaxData['pointsRT'][1,:])**2)+np.sum((rThetaC-uaxData['pointsRC'][1,:])**2))))
    ff.write('---Tension/Directional Stress: Theta, Data, Predicted\n')
    rg=range(uaxData['pointsST'].shape[1])
    for jj in rg:
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(uaxData['pointsST'][0,jj],uaxData['pointsST'][1,jj],sigThetaT[jj]))
    ff.write('---Tension/Directional r-Values: Theta, Data, Predicted\n')
    for jj in rg:
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(uaxData['pointsST'][0,jj],uaxData['pointsRT'][1,jj],rThetaT[jj]))
    ff.write('---Tension Balanced-Biaxial Stress: Data, Predicted\n')
    ff.write('{:.4f}, {:.4f}\n'.format(uaxData['sTb'],sTb))
    ff.write('---Tension Balanced-Biaxial r-Value: Data, Predicted\n')
    ff.write('{:.4f}, {:.4f}\n'.format(uaxData['rTb'],rTb))    
    ff.write('---Compression/Directional stress: Theta, Data, Predicted\n')
    rg=range(uaxData['pointsSC'].shape[1])
    for jj in rg:
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(uaxData['pointsSC'][0,jj],uaxData['pointsSC'][1,jj],sigThetaC[jj])) 
    ff.write('---Compression/Directional r-Values: Theta, Data, Predicted\n')
    for jj in rg:
        ff.write('{:4.1f}, {:.4f}, {:.4f}\n'.format(uaxData['pointsSC'][0,jj],uaxData['pointsRC'][1,jj],rThetaC[jj]))
    ff.write('---Compression Balanced-Biaxial Stress: Data, Predicted\n')
    ff.write('{:.4f}, {:.4f}\n'.format(uaxData['sCb'],sCb))
    ff.write('---Compression Balanced-Biaxial r-Value: Data, Predicted\n')
    ff.write('{:.4f}, {:.4f}\n'.format(uaxData['rCb'],rCb)) 
    ff.write('--------SHYqp-Coefficients\n')
    ff.write('---P-Coeffs\n')
    for jj in range(nP):
        ff.write('{}\n'.format(vCoeff[jj]))
    ff.write('---Q-Coeffs---------------------------------------------\n')    
    for jj in range(nP,nP+nQ):
        ff.write('{}\n'.format(vCoeff[jj]))
    ff.close()
    return
        

def SHYqp_HessGaussCheck(vCoeff,ddMon,nQ,nP):
    degQ=ddMon['nQ']
    degQm1=degQ-1;degPm1=degQ-2
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)
    vD1P=ddMon['vD1P'];vD2P=ddMon['vD2P'];vD3P=ddMon['vD3P']
    vC1P=ddMon['vC1P'];vC2P=ddMon['vC2P'];vC3P=ddMon['vC3P']
    vDD11P=ddMon['vH11P'];vDD12P=ddMon['vH12P'];vDD13P=ddMon['vH13P']
    vDD22P=ddMon['vH22P'];vDD23P=ddMon['vH23P'];vDD33P=ddMon['vH33P']
    vCD11P=ddMon['vCH11P'];vCD12P=ddMon['vCH12P'];vCD13P=ddMon['vCH13P']
    vCD22P=ddMon['vCH22P'];vCD23P=ddMon['vCH23P'];vCD33P=ddMon['vCH33P']
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)
    vD1Q=ddMon['vD1Q'];vD2Q=ddMon['vD2Q'];vD3Q=ddMon['vD3Q']
    vC1Q=ddMon['vC1Q'];vC2Q=ddMon['vC2Q'];vC3Q=ddMon['vC3Q']    
    vDD11Q=ddMon['vH11Q'];vDD12Q=ddMon['vH12Q'];vDD13Q=ddMon['vH13Q']
    vDD22Q=ddMon['vH22Q'];vDD23Q=ddMon['vH23Q'];vDD33Q=ddMon['vH33Q']
    vCD11Q=ddMon['vCH11Q'];vCD12Q=ddMon['vCH12Q'];vCD13Q=ddMon['vCH13Q']
    vCD22Q=ddMon['vCH22Q'];vCD23Q=ddMon['vCH23Q'];vCD33Q=ddMon['vCH33Q']
    if(0):
        nPoints=3500
        np.random.seed(99)
        vu=np.random.normal(0,1,(nPoints,3)) ##nPoints on 2D unit sphere of 3D 
        vNorm=np.sqrt(np.sum(vu**2,axis=1))
        vu[:,0]/=vNorm;vu[:,1]/=vNorm;vu[:,2]/=vNorm
    if(0):
        vu=genConstraintsPoints2DOptPoints(nPoints=200)
        nPoints=vu.shape[0]
    np.random.seed(99)     
    vuRandom=np.random.normal(0,1,(3500,3)) ##nPoints on 2D unit sphere of 3D
    vNorm=np.sqrt(np.sum(vuRandom**2,axis=1))
    vuRandom[:,0]/=vNorm;vuRandom[:,1]/=vNorm;vuRandom[:,2]/=vNorm
    vu=genConstraintsPoints2DOptPoints(nPoints=200)
    vu=np.concatenate((vu[:,0:3],vuRandom),axis=0)
    nPoints=vu.shape[0]    
    vu3=vu[:,2]**2
    ##vMonB=np.zeros((nPoints,nP));vMonD1B=np.zeros((nPoints,nP));vMonD11B=np.zeros((nPoints,nP))
    vMon=np.zeros((nPoints,nP))
    vMonD1=np.zeros((nPoints,nP));vMonD2=np.zeros((nPoints,nP))
    vMonD11=np.zeros((nPoints,nP));vMonD22=np.zeros((nPoints,nP))#;vMonD33=np.zeros((nPoints,nP))
    vMonD12=np.zeros((nPoints,nP));vMonD13=np.zeros((nPoints,nP));vMonD23=np.zeros((nPoints,nP))
    vTerms=np.zeros((nPoints,nPidx))
    vTermsD1=np.zeros((nPoints,nPidx));vTermsD2=np.zeros((nPoints,nPidx));vTermsD3=np.zeros((nPoints,nPidx))
    vTermsD11=np.zeros((nPoints,nPidx));vTermsD22=np.zeros((nPoints,nPidx));vTermsD33=np.zeros((nPoints,nPidx))
    vTermsD12=np.zeros((nPoints,nPidx));vTermsD13=np.zeros((nPoints,nPidx));vTermsD23=np.zeros((nPoints,nPidx))
    for k in range(nP):
        ##vMonB[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]*vu[:,2]**vP[k][2]
        ##vMonD1B[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]*vu[:,2]**vD1P[k][2]
        ##vMonD11B[:,k]=vCD11P[k]*vu[:,0]**vDD11P[k][0]*vu[:,1]**vDD11P[k][1]*vu[:,2]**vDD11P[k][2]
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
        vMonD1[:,k]=vC1P[k]*vu[:,0]**vD1P[k][0]*vu[:,1]**vD1P[k][1]
        vMonD2[:,k]=vC2P[k]*vu[:,0]**vD2P[k][0]*vu[:,1]**vD2P[k][1]
        vMonD11[:,k]=vCD11P[k]*vu[:,0]**vDD11P[k][0]*vu[:,1]**vDD11P[k][1]
        vMonD22[:,k]=vCD22P[k]*vu[:,0]**vDD22P[k][0]*vu[:,1]**vDD22P[k][1]
        ##vMonD33[:,k]=vCD33P[k]*vu[:,0]**vDD33P[k][0]*vu[:,1]**vDD33P[k][1]
        vMonD12[:,k]=vCD12P[k]*vu[:,0]**vDD12P[k][0]*vu[:,1]**vDD12P[k][1]
        vMonD13[:,k]=vCD13P[k]*vu[:,0]**vDD13P[k][0]*vu[:,1]**vDD13P[k][1]
        vMonD23[:,k]=vCD23P[k]*vu[:,0]**vDD23P[k][0]*vu[:,1]**vDD23P[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3P[vPidx[k][1][0]]*vTerms[:,k]
        vTermsD11[:,k]=np.dot(vMonD11[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD22[:,k]=np.dot(vMonD22[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD12[:,k]=np.dot(vMonD12[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD13[:,k]=np.dot(vMonD13[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD23[:,k]=np.dot(vMonD23[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
        vTermsD33[:,k]=vCD33P[vPidx[k][1][0]]*vTerms[:,k]
    ##vPB=np.dot(vMonB,vCoeff[0:nP]);vPD1B=np.dot(vMonD1B,vCoeff[0:nP]);vPD11B=np.dot(vMonD11B,vCoeff[0:nP])        
    vP=vTerms[:,-1]
    vPD1=vTermsD1[:,-1];vPD2=vTermsD2[:,-1];vPD3=vTermsD3[:,-1]
    vPD11=vTermsD11[:,-1];vPD22=vTermsD22[:,-1];vPD33=vTermsD33[:,-1]
    vPD12=vTermsD12[:,-1];vPD13=vTermsD13[:,-1];vPD23=vTermsD23[:,-1]
    for k in range(nPidx-2,-1,-1):
        vP[:]=vTerms[:,k]+vu3*vP
        vPD1[:]=vTermsD1[:,k]+vu3*vPD1
        vPD2[:]=vTermsD2[:,k]+vu3*vPD2
        vPD11[:]=vTermsD11[:,k]+vu3*vPD11
        vPD12[:]=vTermsD12[:,k]+vu3*vPD12
        vPD22[:]=vTermsD22[:,k]+vu3*vPD22 
    for k in range(nPidx-2,0,-1):
        vPD3[:]=vTermsD3[:,k]+vu3*vPD3    
        vPD33[:]=vTermsD33[:,k]+vu3*vPD33
        vPD13[:]=vTermsD13[:,k]+vu3*vPD13
        vPD23[:]=vTermsD23[:,k]+vu3*vPD23
    vPD3[:]*=vu[:,2];vPD13[:]*=vu[:,2];vPD23[:]*=vu[:,2] 
    ##vMonB=np.zeros((nPoints,nQ));vMonD1B=np.zeros((nPoints,nQ));vMonD11B=np.zeros((nPoints,nQ))
    vMon=np.zeros((nPoints,nQ))
    vMonD1=np.zeros((nPoints,nQ));vMonD2=np.zeros((nPoints,nQ))
    vMonD11=np.zeros((nPoints,nQ));vMonD22=np.zeros((nPoints,nQ))
    vMonD12=np.zeros((nPoints,nQ));vMonD13=np.zeros((nPoints,nQ));vMonD23=np.zeros((nPoints,nQ))
    vTerms=np.zeros((nPoints,nQidx))
    vTermsD1=np.zeros((nPoints,nQidx));vTermsD2=np.zeros((nPoints,nQidx));vTermsD3=np.zeros((nPoints,nQidx))
    vTermsD11=np.zeros((nPoints,nQidx));vTermsD22=np.zeros((nPoints,nQidx));vTermsD33=np.zeros((nPoints,nQidx))
    vTermsD12=np.zeros((nPoints,nQidx));vTermsD13=np.zeros((nPoints,nQidx));vTermsD23=np.zeros((nPoints,nQidx))
    for k in range(nQ):
        ##vMonB[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]*vu[:,2]**vQ[k][2]
        ##vMonD1B[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]*vu[:,2]**vD1Q[k][2]
        ##vMonD11B[:,k]=vCD11Q[k]*vu[:,0]**vDD11Q[k][0]*vu[:,1]**vDD11Q[k][1]*vu[:,2]**vDD11Q[k][2]
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
        vMonD1[:,k]=vC1Q[k]*vu[:,0]**vD1Q[k][0]*vu[:,1]**vD1Q[k][1]
        vMonD2[:,k]=vC2Q[k]*vu[:,0]**vD2Q[k][0]*vu[:,1]**vD2Q[k][1]
        vMonD11[:,k]=vCD11Q[k]*vu[:,0]**vDD11Q[k][0]*vu[:,1]**vDD11Q[k][1]
        vMonD22[:,k]=vCD22Q[k]*vu[:,0]**vDD22Q[k][0]*vu[:,1]**vDD22Q[k][1]
        vMonD12[:,k]=vCD12Q[k]*vu[:,0]**vDD12Q[k][0]*vu[:,1]**vDD12Q[k][1]
        vMonD13[:,k]=vCD13Q[k]*vu[:,0]**vDD13Q[k][0]*vu[:,1]**vDD13Q[k][1]
        vMonD23[:,k]=vCD23Q[k]*vu[:,0]**vDD23Q[k][0]*vu[:,1]**vDD23Q[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD1[:,k]=np.dot(vMonD1[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD2[:,k]=np.dot(vMonD2[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD3[:,k]=vC3Q[vQidx[k][1][0]]*vTerms[:,k]
        vTermsD11[:,k]=np.dot(vMonD11[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD22[:,k]=np.dot(vMonD22[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD12[:,k]=np.dot(vMonD12[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD13[:,k]=np.dot(vMonD13[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD23[:,k]=np.dot(vMonD23[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
        vTermsD33[:,k]=vCD33Q[vQidx[k][1][0]]*vTerms[:,k]
    ##vQB=np.dot(vMonB,vCoeff[nP:nP+nQ]);vQD1B=np.dot(vMonD1B,vCoeff[nP:nP+nQ]);vQD11B=np.dot(vMonD11B,vCoeff[nP:nP+nQ])   
    ##print("vTerms\n",vTerms[500,:])    
    vQ=vTerms[:,-1]
    vQD1=vTermsD1[:,-1];vQD2=vTermsD2[:,-1];vQD3=vTermsD3[:,-1]
    vQD11=vTermsD11[:,-1];vQD22=vTermsD22[:,-1];vQD33=vTermsD33[:,-1]
    vQD12=vTermsD12[:,-1];vQD13=vTermsD13[:,-1];vQD23=vTermsD23[:,-1]
    for k in range(nQidx-2,-1,-1):
        vQ[:]=vTerms[:,k]+vu3*vQ
        vQD1[:]=vTermsD1[:,k]+vu3*vQD1
        vQD2[:]=vTermsD2[:,k]+vu3*vQD2
        vQD11[:]=vTermsD11[:,k]+vu3*vQD11
        vQD12[:]=vTermsD12[:,k]+vu3*vQD12
        vQD22[:]=vTermsD22[:,k]+vu3*vQD22 
    for k in range(nQidx-2,0,-1):
        vQD3[:]=vTermsD3[:,k]+vu3*vQD3    
        vQD33[:]=vTermsD33[:,k]+vu3*vQD33
        vQD13[:]=vTermsD13[:,k]+vu3*vQD13
        vQD23[:]=vTermsD23[:,k]+vu3*vQD23
    vQD3[:]*=vu[:,2];vQD13[:]*=vu[:,2];vQD23[:]*=vu[:,2] 
    ##vPhiB=1.0-(degPm1*vPB+degQm1*vQB)
    ##vPhi2B=(degPm1**2-1)*vPB+(degQm1**2-1)*vQB-1.0
    ##H11B=vPhiB+vPhi2B*vu[:,0]**2-(degPm1*(2*vu[:,0]*vPD1B)+degQm1*(2*vu[:,0]*vQD1B))+vPD11B+vQD11B
    vPhi=1.0-(degPm1*vP+degQm1*vQ)
    vPhi2=(degQm1**2-1)*vP+(degQ**2-1)*vQ-1.0
    H11a=vPhi+vPhi2*vu[:,0]**2-(degPm1*(2*vu[:,0]*vPD1)+degQm1*(2*vu[:,0]*vQD1))+vPD11+vQD11
    H22a=vPhi+vPhi2*vu[:,1]**2-(degPm1*(2*vu[:,1]*vPD2)+degQm1*(2*vu[:,1]*vQD2))+vPD22+vQD22
    H33a=vPhi+vPhi2*vu3-(degPm1*(2*vu[:,2]*vPD3)+degQm1*(2*vu[:,2]*vQD3))+vPD33+vQD33
    H12a=vPhi2*vu[:,0]*vu[:,1]-(degPm1*(vu[:,0]*vPD2+vPD1*vu[:,1])+degQm1*(vu[:,0]*vQD2+vQD1*vu[:,1]))+vPD12+vQD12
    H13a=vPhi2*vu[:,0]*vu[:,2]-(degPm1*(vu[:,0]*vPD3+vPD1*vu[:,2])+degQm1*(vu[:,0]*vQD3+vQD1*vu[:,2]))+vPD13+vQD13
    H23a=vPhi2*vu[:,1]*vu[:,2]-(degPm1*(vu[:,1]*vPD3+vPD2*vu[:,2])+degQm1*(vu[:,1]*vQD3+vQD2*vu[:,2]))+vPD23+vQD23
    print("Hessian leading principal minors (min_value over the unit sphere):")
    k=np.argmin(H11a);m1=H11a[k];print("min det_1 = ",m1)
    ddet=H11a*H22a-H12a*H12a
    k=np.argmin(ddet);m2=ddet[k];print("min det_2 = ",m2)
    ddet=H33a*ddet+H13a*(H12a*H23a-H22a*H13a)+H23a*(H12a*H13a-H11a*H23a)
    k=np.argmin(ddet);m3=ddet[k];print("min det_3 = ",m3)
    sq3=np.sqrt(3.0);sq32=np.sqrt(1.5)
    H11=(2.0/3.0)*H11a
    H12=H12a/sq3-H11a/3
    H13=2*H13a/sq3 ## 2x tensor component  
    H22=H11a/6.0-H12a/sq3+H22a/2
    H23=2*0.5*(H23a-H13a/sq3) ##2x tensor component 
    H33=4*0.5*H33a ##4x tensor component
    H11s=H22*H33-H23*H23
    H12s=H23*H13-H12*H33
    H13s=H12*H23-H22*H13
    H22s=H11*H33-H13*H13
    H23s=H12*H13-H11*H23
    H33s=H11*H22-H12*H12
    G1a=sq32*(vu[:,0]*vPhi+vPD1+vQD1)
    G2a=sq32*(vu[:,1]*vPhi+vPD2+vQD2)
    G3=2*0.5*sq3*(vu[:,2]*vPhi+vPD3+vQD3) ##2x tensor component
    G1=np.sqrt(2/3.0)*G1a
    G2=G2a/np.sqrt(2)-G1a/np.sqrt(6)
    vNorm=np.sqrt(G1**2+G2**2+G3**2)
    G1[:]/=vNorm;G2[:]/=vNorm;G3[:]/=vNorm
    KG=H11s*G1*G1+H22s*G2*G2+H33s*G3*G3+2.0*(H12s*G1*G2+H13s*G1*G3+H23s*G2*G3)
    KG=sq32*KG/vNorm**2
    mm=np.min(KG);print("min(Gaussian curvature) = ",mm)
    if(mm<0):
        k=np.argmin(KG);print("KG[{}] = {}; location: u=[{},{},{}]".format(k,KG[k],vu[k,0],vu[k,1],vu[k,2]))
        ###test for numerical instability
        ##kg=SHYqp_KGCheckPoint(np.array([-vu[k,0],-vu[k,1],vu[k,2]]).reshape(1,3),vCoeff,ddMon,nQ,nP)
        ##print("kg = ",kg)
        ##A=np.array([[H11s[k],H12s[k],H13s[k]],[H12s[k],H22s[k],H23s[k]],[H13s[k],H23s[k],H33s[k]]]).reshape((3,3))
        ##print("eigenvalues = ",np.linalg.eigvalsh(A))
    return (m1,m2,m3,mm)
 

 


def SHYqp_surf_Plot(vCoeff,ddMon,nQ,nP,name,savePng=False):
    degQ=ddMon['nQ']
    degQm1=degQ-1;degPm1=degQ-2;sq3=np.sqrt(3.0)
    vP=ddMon['vP'];vPidx=ddMon['vPidx'];nPidx=len(vPidx)
    vQ=ddMon['vQ'];vQidx=ddMon['vQidx'];nQidx=len(vQidx)
    n2=50;n1=4*n2;nPoints=n1*n2
    t1=np.linspace(0,2*np.pi,n1)
    t2=np.linspace(0,0.5*np.pi,n2)
    vt1,vt2=np.meshgrid(t1,t2)
    vt1=vt1.reshape(nPoints)
    vt2=vt2.reshape(nPoints)
    ct1=np.cos(vt1);st1=np.sin(vt1)
    ct2=np.cos(vt2);st2=np.sin(vt2)
    ####print("min vtt = ",np.min(1.0+2*ct2*ct2-ct1*st1*st2*st2))
    vtt=np.sqrt(1.0+2*ct2*ct2-ct1*st1*st2*st2)
    vu=np.zeros((nPoints,3))
    vu[:,0]=(0.5*st2*(2*ct1-st1))/vtt
    vu[:,1]=(0.5*sq3*st2*st1)/vtt
    vu[:,2]=(sq3*ct2)/vtt
    vu3=vu[:,2]**2
    vMon=np.zeros((nPoints,nP));vTerms=np.zeros((nPoints,nPidx))
    for k in range(nP):
        vMon[:,k]=vu[:,0]**vP[k][0]*vu[:,1]**vP[k][1]
    for k in range(nPidx):
        vTerms[:,k]=np.dot(vMon[:,vPidx[k][1][0]:vPidx[k][1][-1]+1],vCoeff[vPidx[k][1][0]:vPidx[k][1][-1]+1])
    vP=vTerms[:,-1]
    for k in range(nPidx-2,-1,-1):
        vP[:]=vTerms[:,k]+vu3*vP
    vMon=np.zeros((nPoints,nQ));vTerms=np.zeros((nPoints,nQidx))
    for k in range(nQ):
        vMon[:,k]=vu[:,0]**vQ[k][0]*vu[:,1]**vQ[k][1]
    for k in range(nQidx):
        vTerms[:,k]=np.dot(vMon[:,vQidx[k][1][0]:vQidx[k][1][-1]+1],vCoeff[nP+vQidx[k][1][0]:nP+vQidx[k][1][-1]+1])
    vQ=vTerms[:,-1]
    for k in range(nQidx-2,-1,-1):
        vQ[:]=vTerms[:,k]+vu3*vQ        
    RR=1.0/(vtt*(1+vP+vQ))
    x=(RR*st2*ct1).reshape((n2,n1))
    y=(RR*st2*st1).reshape((n2,n1))
    z=(RR*ct2).reshape((n2,n1))
    fg=plt.figure();ax=fg.add_subplot(1,1,1,projection='3d')
    vColor=np.array([200/255,210/255,220/255])
    ax.plot_surface(x,y,z,color=vColor,alpha=0.45,linewidth=0.3, edgecolors=0.55*vColor)
    ax.plot_surface(x,y,-z,color=vColor,alpha=0.45,linewidth=0.3, edgecolors=0.55*vColor)
    ax.set_xlabel(r'$\overline{\sigma}_{xx}$',fontsize=14)
    ax.set_ylabel(r'$\overline{\sigma}_{yy}$',fontsize=14)
    if(savePng):
        fg.savefig('{name}.{ext}'.format(name=figDir+name+'_SHYqp_deg'+str(degQ)+'_Surf',ext='png'),dpi=300,bbox_inches='tight')
    return
    
'''
function: 'testUaxInterpPlot'
--Example plot of interpolated directional values 
--Used only for testing/illustration of the influence of the shape parameter
'''
def testUaxInterpPlot(data,savePng=False):
    uaxData=uaxLambda(data)
    nData=uaxData['pointsST'].shape[1]
    fg=plt.figure()
    ax1=fg.add_subplot(2,1,1);ax2=fg.add_subplot(2,1,2)
    vScale=[0.0,0.6,1.0]
    vLine=['--','-',':']
    for JJ in range(len(vScale)):    
        ##Lshape=data['shapeUAX']*uaxData['lambdaSTmax']
        LshapeS=vScale[JJ]*uaxData['lambdaSTmax']
        LshapeR=vScale[JJ]*uaxData['lambdaRTmax']
        for kk in range(1,nData):
            vv=curveSeg(uaxData['pointsST'][:,kk-1].reshape(2,1),uaxData['tanST'][:,kk-1].reshape(2,1),
                        uaxData['pointsST'][:,kk].reshape(2,1),uaxData['tanST'][:,kk].reshape(2,1),LshapeS)
            ax1.plot(vv[0,:],vv[1,:],color='k',linestyle=vLine[JJ])
            vv=curveSeg(uaxData['pointsRT'][:,kk-1].reshape(2,1),uaxData['tanRT'][:,kk-1].reshape(2,1),
                        uaxData['pointsRT'][:,kk].reshape(2,1),uaxData['tanRT'][:,kk].reshape(2,1),LshapeR)
            ax2.plot(vv[0,:],vv[1,:],color='k',linestyle=vLine[JJ])
    ax1.plot(uaxData['pointsST'][0,:],uaxData['pointsST'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    ax2.plot(uaxData['pointsRT'][0,:],uaxData['pointsRT'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    ax1.set_xticks([0,15,30,45,60,75,90])
    ax1.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=10) 
    ax1.grid()
    ##ax1.text(92,0.95*min(uaxData['pointsST'][1,:]),r'$\theta$',fontsize=14)
    ax1.text(-4,0.92*max(uaxData['pointsST'][1,:]),r'$\overline{\sigma}_{\theta}$',fontsize=14)
    ax2.set_xticks([0,15,30,45,60,75,90])
    ax2.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=10) 
    ax2.grid()
    ax2.text(-4,0.8*max(uaxData['pointsRT'][1,:]),r'$r_{\theta}$',fontsize=14)
    ax2.text(91.25,0.95*min(uaxData['pointsRT'][1,:]),r'$\theta$',fontsize=14)
    if(savePng):
        fg.savefig('{name}.{ext}'.format(name=figDir+'Bezier_uaxExample',ext='png'),dpi=300,bbox_inches='tight')
    return

'''
function:'protoBez5YS_uaxPlot'
--Plots the protoBez5YS interpolated directional properties 
--Utility for assesing the overall aspect (as determined by the shape parameters lambdaMax)
'''
def protoBez5YS_uaxPlot(uaxData,savePng=False):
    fg=plt.figure()
    ax1=fg.add_subplot(2,1,1);ax2=fg.add_subplot(2,1,2)
    nData=uaxData['pointsST'].shape[1]
    LshapeS=uaxData['shapeUAX']*uaxData['lambdaSTmax']
    LshapeR=uaxData['shapeUAX']*uaxData['lambdaRTmax']
    for kk in range(1,nData):
        vv=curveSeg(uaxData['pointsST'][:,kk-1].reshape(2,1),uaxData['tanST'][:,kk-1].reshape(2,1),
                    uaxData['pointsST'][:,kk].reshape(2,1),uaxData['tanST'][:,kk].reshape(2,1),LshapeS)
        ax1.plot(vv[0,:],vv[1,:],color='k',linestyle='-',linewidth=lineWidthUax)
        vv=curveSeg(uaxData['pointsRT'][:,kk-1].reshape(2,1),uaxData['tanRT'][:,kk-1].reshape(2,1),
                    uaxData['pointsRT'][:,kk].reshape(2,1),uaxData['tanRT'][:,kk].reshape(2,1),LshapeR)
        ax2.plot(vv[0,:],vv[1,:],color='k',linestyle='-',linewidth=lineWidthUax)
    ax1.plot(uaxData['pointsST'][0,:],uaxData['pointsST'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    ax2.plot(uaxData['pointsRT'][0,:],uaxData['pointsRT'][1,:],'bs',markerfacecolor='b',markersize=6,label='Exp data')
    ax1.set_xticks([0,15,30,45,60,75,90])
    ax1.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=10) 
    ax1.grid()
    ##ax1.text(92,0.95*min(uaxData['pointsST'][1,:]),r'$\theta$',fontsize=14)
    ##ax1.text(-4,0.88*max(uaxData['pointsST'][1,:]),r'$\overline{\sigma}_{\theta}$',fontsize=14)
    ax1.text(0.022,0.88,r'$\overline{\sigma}_{\theta}$',fontsize=14,
    horizontalalignment='center',verticalalignment='center',transform = ax1.transAxes)
    ax2.set_xticks([0,15,30,45,60,75,90])
    ax2.set_xticklabels(['$0^o$','$15^o$','$30^o$','$45^o$','$60^o$','$75^o$','$90^o$'],fontsize=10) 
    ax2.grid()
    ##ax2.text(-4,0.8*max(uaxData['pointsRT'][1,:]),r'$r_{\theta}$',fontsize=14)
    ax2.text(0.022,0.92,r'$r_{\theta}$',fontsize=14,
    horizontalalignment='center',verticalalignment='center',transform = ax2.transAxes)
    ##ax2.text(91.25,0.95*min(uaxData['pointsRT'][1,:]),r'$\theta$',fontsize=14)
    ##ax2.text(1,0,r'$\theta$',fontsize=14,
    ##horizontalalignment='center',verticalalignment='center',transform = ax2.transAxes)
    if(not uaxData['assym']):
        if(savePng):
            fg.savefig('{name}.{ext}'.format(name=figDir+uaxData['name']+'_Bezier5YS_Uax',ext='png'),dpi=300,bbox_inches='tight')
        return
    nData=uaxData['pointsSC'].shape[1]
    LshapeS=uaxData['shapeUAX']*uaxData['lambdaSCmax']
    LshapeR=uaxData['shapeUAX']*uaxData['lambdaRCmax']
    for kk in range(1,nData):
        vv=curveSeg(uaxData['pointsSC'][:,kk-1].reshape(2,1),uaxData['tanSC'][:,kk-1].reshape(2,1),
                    uaxData['pointsSC'][:,kk].reshape(2,1),uaxData['tanSC'][:,kk].reshape(2,1),LshapeS)
        ax1.plot(vv[0,:],vv[1,:],color='k',linestyle='--',linewidth=lineWidthUax)
        vv=curveSeg(uaxData['pointsRC'][:,kk-1].reshape(2,1),uaxData['tanRC'][:,kk-1].reshape(2,1),
                    uaxData['pointsRC'][:,kk].reshape(2,1),uaxData['tanRC'][:,kk].reshape(2,1),LshapeR)
        ax2.plot(vv[0,:],vv[1,:],color='k',linestyle='--',linewidth=lineWidthUax)        
    ax1.plot(uaxData['pointsSC'][0,:],uaxData['pointsSC'][1,:],'bo',markerfacecolor='b',markersize=6.5,label='Exp data')
    ax2.plot(uaxData['pointsRC'][0,:],uaxData['pointsRC'][1,:],'bo',markerfacecolor='b',markersize=6.5,label='Exp data')
    maxS=np.max([np.max(uaxData['pointsST'][1,:]),np.max(uaxData['pointsSC'][1,:])])
    minS=np.min([np.min(uaxData['pointsST'][1,:]),np.min(uaxData['pointsSC'][1,:])])
    ax1.set_ylim([0.97*minS,1.03*maxS])
    maxR=np.max([np.max(uaxData['pointsRT'][1,:]),np.max(uaxData['pointsRC'][1,:])])
    minR=np.min([np.min(uaxData['pointsRT'][1,:]),np.min(uaxData['pointsRC'][1,:])])
    ax2.set_ylim([0.05*minR,1.05*maxR])
    if(savePng):
        fg.savefig('{name}.{ext}'.format(name=figDir+uaxData['name']+'_Bezier5YS_Uax',ext='png'),dpi=300,bbox_inches='tight')
    return


def plotUAX3D(uaxData,ax):
    rad=np.pi/180
    lbd=uaxData['lambdaSTmax']*uaxData['shapeUAX']
    vsT=uaxData['pointsST']
    vsTtan=uaxData['tanST']
    for kk in range(vsT.shape[1]-1):
        vv=curveSeg(vsT[:,kk].reshape((2,1)),vsTtan[:,kk].reshape((2,1)),vsT[:,kk+1].reshape((2,1)),vsTtan[:,kk+1].reshape((2,1)),lbd)
        vv[0,:]=rad*vv[0,:]
        zz=[vv[1,:]*(np.cos(vv[0,:]))**2,vv[1,:]*(np.sin(vv[0,:]))**2,vv[1,:]*np.cos(vv[0,:])*np.sin(vv[0,:])]
        ax.plot(zz[0],zz[1],zz[2],color='b')
        ax.plot(zz[0],zz[1],-zz[2],color='b')
        if(not uaxData['assym']):
            ax.plot(-zz[0],-zz[1],-zz[2],color='b')
            ax.plot(-zz[0],-zz[1],zz[2],color='b')
    if(not uaxData['assym']):return
    lbd=uaxData['lambdaSCmax']*uaxData['shapeUAX']
    vsT=uaxData['pointsSC']
    vsTtan=uaxData['tanSC']
    for kk in range(vsT.shape[1]-1):
        vv=curveSeg(vsT[:,kk].reshape((2,1)),vsTtan[:,kk].reshape((2,1)),vsT[:,kk+1].reshape((2,1)),vsTtan[:,kk+1].reshape((2,1)),lbd)
        vv[0,:]=rad*vv[0,:]
        zz=[-vv[1,:]*(np.cos(vv[0,:]))**2,-vv[1,:]*(np.sin(vv[0,:]))**2,-vv[1,:]*np.cos(vv[0,:])*np.sin(vv[0,:])]
        ax.plot(zz[0],zz[1],zz[2],color='b')
        ax.plot(zz[0],zz[1],-zz[2],color='b')
    return   
    
'''
function:'protoBez5YS_Plot'
-Plots the yield surface (by plane sections) of the Bezier5YS proto-model
'''        
def protoBez5YS_Plot(maxLambda,shapeBAX,vPatch,uaxData,savePng=False):
    fg=plt.figure();ax=fg.add_subplot(1,1,1,projection='3d');linewidth=0.4
    lbd=maxLambda*shapeBAX
    #print('Bez5YS lbd:\n',lbd)
    pp=1
    #vColors=['k','g','r','b','y','m'] ##for visualizing each patch; below use this with vColors
    for dd in vPatch:
        for kk in range(5):
            vv=curveSeg3(dd['Points'][:,kk].reshape(3,1),dd['Tangents'][:,kk].reshape(3,1),dd['Points'][:,kk+1].reshape(3,1),dd['Tangents'][:,kk+1].reshape(3,1),lbd[kk],lbd[kk+1])
            ax.plot(vv[0,:],vv[1,:],vv[2,:],color='k',linewidth=pp*linewidth)
            ax.plot(vv[0,:],vv[1,:],-vv[2,:],color='k',linewidth=pp*linewidth)
            #ax.plot(vv[0,:],vv[1,:],vv[2,:],color=vColors[kk%6],linewidth=pp**linewidth)
            #ax.plot(vv[0,:],vv[1,:],-vv[2,:],color=vColors[kk%6],linewidth=pp*linewidth)
        vv=curveSeg3(dd['Points'][:,5].reshape(3,1),dd['Tangents'][:,5].reshape(3,1),dd['Points'][:,0].reshape(3,1),dd['Tangents'][:,0].reshape(3,1),lbd[5],lbd[0])
        ax.plot(vv[0,:],vv[1,:],vv[2,:],color='k',linewidth=pp*linewidth)
        ax.plot(vv[0,:],vv[1,:],-vv[2,:],color='k',linewidth=pp*linewidth)
        #ax.plot(vv[0,:],vv[1,:],vv[2,:],color=vColors[5],linewidth=pp*linewidth)
        #ax.plot(vv[0,:],vv[1,:],-vv[2,:],color=vColors[5],linewidth=pp*linewidth) 
    ax.set_xlabel(r'$\overline{\sigma}_{xx}$',fontsize=14)
    ax.set_ylabel(r'$\overline{\sigma}_{yy}$',fontsize=14)
    plotUAX3D(uaxData,ax)
    if(savePng):
        #ax.view_init(elev=90,azim=-90) ## view from top (z-axis) to check biaxial shape
        #ax.view_init(elev=15,azim=115) ## side view of discontinuity (in case of arbitrary six shape parameters)
        #shapeBAX=str(shapeBAX)
        shapeBAX=str(shapeBAX[0])
        shapeBAX='Lambda_'+shapeBAX[0]+'p'+shapeBAX[2:]
        fg.savefig('{name}.{ext}'.format(name=figDir+uaxData['name']+'_Bezier5YS_'+shapeBAX,ext='png'),dpi=300,bbox_inches='tight')    
    return
    


def testProtoYSBiax_Plot(lbd,data):
    fg=plt.figure();ax=fg.add_subplot(1,1,1)
    dd=data[0]##;print('dd=\n',dd)
    for kk in range(5):
        vv=curveSeg3(dd['Points'][:,kk].reshape(3,1),dd['Tangents'][:,kk].reshape(3,1),dd['Points'][:,kk+1].reshape(3,1),dd['Tangents'][:,kk+1].reshape(3,1),lbd)
        ax.plot(vv[0,:],vv[1,:],color='k')
    vv=curveSeg3(dd['Points'][:,5].reshape(3,1),dd['Tangents'][:,5].reshape(3,1),dd['Points'][:,0].reshape(3,1),dd['Tangents'][:,0].reshape(3,1),lbd)
    ax.plot(vv[0,:],vv[1,:],color='k')
    nn=101;nSamples=nn*nn
    vx=np.linspace(-1.2,1.2,nn)
    vy=np.linspace(-1.2,1.2,nn)
    x,y=np.meshgrid(vx,vy)
    zz=x**2+y**2-x*y###Mises
    ##x=x.reshape(nSamples)
    ##y=y.reshape(nSamples)
    ax.contour(y,x,zz,levels=[1],colors='r')
    ax.set_aspect('equal')
    ax.grid()    
    return
    
'''
function:'bxCurve2'
--Compares Bezier interpolated baxial curve to SHYqp-predicted curve 
'''
def bxCurve2(vCoeff,deg,nnP):
    nQ=deg;nQ1=nQ+1;nP=nQ-1
    vP=[(nP-k,k) for k in range(0,nQ)]
    vQ=[(nQ-k,k) for k in range(0,nQ1)]
    nx=400;nSamples=nx*nx
    x=np.linspace(-1.5,1.5,nx)
    y=np.linspace(-1.5,1.5,nx)    
    vx,vy=np.meshgrid(x,y)
    vx=vx.reshape(nSamples)
    vy=vy.reshape(nSamples)
    sq32=np.sqrt(1.5);sq6=np.sqrt(6.0);sq2=np.sqrt(2.0)
    vx=(2.0*vx-vy)/sq6;vy=vy/sq2
    vMod=np.sqrt(vx**2+vy**2)
    vx=vx/vMod;vy=vy/vMod
    u1u2=np.zeros((nSamples,2))
    u1u2[:,0]=vx;u1u2[:,1]=vy
    zz=np.zeros(nSamples)
    vPowersP1=np.zeros((nSamples,nQ))
    vPowersP2=np.zeros((nSamples,nQ))
    vPowersQ1=np.zeros((nSamples,nQ1))
    vPowersQ2=np.zeros((nSamples,nQ1))
    for k in range(nQ):
        vPowersP1[:,k]=u1u2[:,0]**vP[k][0]
        vPowersP2[:,k]=u1u2[:,1]**vP[k][1]
    for k in range(nQ1):
        vPowersQ1[:,k]=u1u2[:,0]**vQ[k][0]
        vPowersQ2[:,k]=u1u2[:,1]**vQ[k][1]
    ###print(vCoeff[0],vCoeff[1],vCoeff[nnP],vCoeff[nnP+1])    
    zz[:]=sq32*vMod[:]*(1.0+np.dot(vPowersP1*vPowersP2,vCoeff[0:nQ])+np.dot(vPowersQ1*vPowersQ2,vCoeff[nnP:nnP+nQ1]))    
    ##for k in range(nSamples):
    ##    powerP=u1u2[k]**vP
    ##    powerQ=u1u2[k]**vQ
    ##    zz[k]=1.0+np.dot(bxCoeff[0:nQ],powerP[:,0]*powerP[:,1])+np.dot(bxCoeff[nQ:],powerQ[:,0]*powerQ[:,1])
    ##zz*=sq32*vMod
    fg=plt.figure();ax=fg.add_subplot(1,1,1)
    ax.contour(x,y,zz.reshape(nx,nx),levels=[1])
    ax.set_aspect('equal')
    ax.grid()
    return    



### function 'getCoeff': utility function
### use it to extract back the SHYqp coefficients from a '*_Err_and_Coeff.txt' file    
def getCoeff(fName):
    ff=open(fName,'r')
    line=ff.readline()
    line=ff.readline()
    deg=int(line.strip().split(':')[1])
    print('degQ = ',deg)
    while('P-Coeffs' not in line):
        line=ff.readline()
    line=ff.readline()    
    vCoeff=[]    
    while('Q-Coeffs' not in line):
        vCoeff.append(float(line.strip()))
        line=ff.readline()
    print('P-coeff: Done')        
    while(True):
        line=ff.readline().strip()
        if(line):
            vCoeff.append(float(line))
        else:
            break
    print('Q-coeff: Done')        
    return deg,np.array(vCoeff)       



if(not (osp.exists(figDir) and osp.isdir(figDir))):
    print("The local folder for saving reports and figures was not found")
    print("Make sure a folder \'FIGS\' exists at the same location as this script")
    print("Calculations aborted");exit()
    

if __name__ == "__main__": 
    ###Read mechanical data and other global parameters from text file  
    data=readData('mat000File.txt')
    ### echo data 
    prtData(data)  
    ###Generate the sequence of points and tangents from directional data 
    uaxData=uaxLambda(data)
    ###################plot Bezier5YS
    savePngProtoModel=True
    protoBez5YS_uaxPlot(uaxData,savePngProtoModel)
    lbd,vPatch=protoData(uaxData,31)
    print("max sectional shape parameter = ",lbd)
    protoBez5YS_Plot(lbd,data['shapeBAX'],vPatch,uaxData,savePngProtoModel)
    if(1):###calculate SHYqp parameters and plots 
        ##qpSolver='quadprog'
        qpSolver='cvxopt'        
        if(data['assym']):
            vCoeff,ddMon,nQ,nP=dataFitSHYqp(data,uaxData,lbd,qpSolver,nEquator=200,epsilon=0.01)
        else:
            vCoeff,ddMon,nQ,nP=dataFitSHYqpSymm(data,uaxData,lbd,qpSolver,nEquator=200,epsilon=0.01)
        cvxCheck=SHYqp_HessGaussCheck(vCoeff,ddMon,nQ,nP)
        savePngSHYqp=True
        SHYqp_Predictions(uaxData,vCoeff,ddMon,nQ,nP,qpSolver,cvxCheck)
        SHYqp_uax_Plot(uaxData,vCoeff,ddMon,nQ,nP,savePngSHYqp)
        SHYqp_bax_Plot(vCoeff,ddMon,nQ,nP,data['assym'],uaxData['name'],savePngSHYqp)
        SHYqp_surf_Plot(vCoeff,ddMon,nQ,nP,uaxData['name'],savePngSHYqp)    
    plt.show()
