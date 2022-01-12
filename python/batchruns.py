import os
import tempfile
import numpy as np
import itertools
from subprocess import check_output

def changeFile(inFile,outFile,inString,outString):
    fin = open(inFile,"rt")
    
    lines=[]
    for line in fin:
        lines.append(line)
        
    fin.close()
    fout = open(outFile,"wt")

    for line in lines:
        fout.write(line.replace(inString,outString))
    
    fout.close()






inputFile='../namelist.in'
tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
tmpFile=tmpFileObj.name

runsToDo = [['hm_flag=.false.', 'hm_flag=.true.'],\
            ['mode2_flag=.false.', 'mode2_flag=.true.'], \
            ['break_flag=.false.', 'break_flag=.true.'] ]
# runsToDo = [['hm_flag=.false.', 'hm_flag=.true.'],\
#             ['mode2_flag=.false.', 'mode2_flag=.true.'], \
#             ['break_flag=.false.', 'break_flag=.true.'],\
#             ['aer=1','aer=2'] ]

columns1=['hm','mode2','break']
#columns1=['hm','mode2','break','aerosol']


(r,c)=np.shape(runsToDo)

loops = c**r
k=0
for elements in itertools.product(*runsToDo[:loops]):
    n=str(k)
    print(str(k) )
    print(elements)
    
    
    
    changeFile(inputFile,tmpFile,runsToDo[0][0],elements[0])
    changeFile(tmpFile,tmpFile,runsToDo[1][0],elements[1])
    changeFile(tmpFile,tmpFile,runsToDo[2][0],elements[2])
    #changeFile(tmpFile,tmpFile,runsToDo[3][0],elements[3])
    changeFile(tmpFile,tmpFile,\
                'outputfile = \'/tmp/output1.nc\'',\
                'outputfile = \'/tmp/output' + n.zfill(3) + '.nc\'')
                
                
    # run the model here
    str1='./main.exe ' + tmpFile       
    result = check_output(str1, shell=True,cwd='../').decode()
    print (str1)
    k += 1
    
tmpFileObj.close()
os.unlink(tmpFileObj.name)    
    
    