#!/usr/bin/env python
# coding: utf-8
#!/usr/bin/env python
# coding: utf-8
import sys
import os
inputfile=sys.argv[1]
outputfile=sys.argv[2]
subclass,subclass_unique=[],[]
if os.path.isfile(inputfile) :
    with open(inputfile,'r') as file:
        lines=file.readlines()
        with open (outputfile,"w")as out:
            for line in lines[:]:
                splitedline=line.split("\t")
                outline=splitedline.copy()
                ppt=splitedline[17] #PPT column number
                pbs=splitedline[22] #PBS column number
                ltr1=int(splitedline[6])-int(splitedline[5]) #left LTR end then start column number
                ltr2=int(splitedline[9])-int(splitedline[8]) #right LTR end then start column number
                ltr_lenght= int(ltr1)+int(ltr2)
                full_lenght=int(splitedline[4])  #LTR-RT full length
                ir_lenght=int(full_lenght-ltr_lenght)
                clade=splitedline[32] #family column
                subclass=splitedline[36] #domains column
                subclass.replace(" ","")
                subb_splited=[sub[:sub.index("|")] if "|" in subclass else 'NONE' for sub in subclass.split(" ")]
                subb_splited+=" " 
                non_existing=[elem for elem in ["GAG","PROT","RT","RH","INT"] if elem not in subb_splited]
                subb="|".join(subb_splited).strip()
                # get all the region that have pbs and ppt
                if  pbs !='' and ppt !='':
                    if len(subb_splited[:-1])== 0:
                        #get TRIM
                        if ltr2 in  range(100,251) and  ltr1 in  range(100,251) and ir_lenght in range(100,301):
                            outline=[ f'Nonautonomous:{clade}1']+outline  
                         #get LARD    
                        elif (ltr2 in  range(4000,4501) and  ltr1 in  range(4000,4501) and  ir_lenght in range(3000,3601)):  
                            outline=[ f'Nonautonomous:{clade}']+outline 
                        else:
                            #not trim nor a lard and no domaine exsist
                             outline=[ f'Nonautonomous:{clade}']+outline 
                    else:
                        if "GAG" in subb_splited :
                            gag_ing=subb_splited.index('GAG')
                             #get autonoumos gypsy 
                            if ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "RT"==subb_splited[gag_ing+2]
                                and "RH"==subb_splited[gag_ing+3] and "INT"==subb_splited[gag_ing+4]):
                                    #append to list
                                    outline=[ f'Autonomous:{clade}']+outline
                            #get autonoumos copia         
                            elif ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "INT"==subb_splited[gag_ing+2]
                                and "RT"==subb_splited[gag_ing+3] and "RH"==subb_splited[gag_ing+4]):
                                     #append to list
                                     outline=[ f'Autonomous:{clade}']+outline
                            #get  TR-GAG           
                            elif len(subb_splited[:-1])==1 and "GAG"==subb_splited[gag_ing] :
                                
                                outline=[ f'Nonautonomous:TR-GAG']+outline
                            elif len(non_existing)==0:
                                # all domaine exist but different order
                                outline=[ f'Nonautonomous:{clade}']+outline
                            else:
                                outline=[ f'Nonautonomous:{clade}']+outline
                        # gag is not in domaine and prot is
                        elif  "GAG" not in subb_splited  and 'PROT' in subb_splited: 
                            prot_ing=subb_splited.index('PROT')
                            #get bare-2
                            if  ("PROT"==subb_splited[prot_ing] and "INT"==subb_splited[prot_ing+1]
                                and "RT"==subb_splited[prot_ing+2] and "RH"==subb_splited[prot_ing+3]):
                                    outline=[ f'Nonautonomous:BARE-2']+outline
                            elif ("PROT" in subb_splited and "INT" in subb_splited
                                and "RT"in subb_splited and "RH" in subb_splited):
                                # gag in not in domaine and prot in domaine but in diff order
                                outline=[ f'Nonautonomous:BARE-2']+outline
                            else:
                                # gag inin not in domaine and prot in domaine but some domaine are missing
                                outline=[ f'Nonautonomous:{clade}']+outline
                        else:
                            # some domains are missing or  the order in not the same
                            outline=[ f'Nonautonomous:{clade}']+outline
                
                
                #second part of the code--------------------------------------------------------------
                elif   pbs =='' and ppt !='': 
                    
                    if len(subb_splited[:-1])== 0:
                            #not trim nor a lard and no domaine exsist
                             outline=[ f'Nonautonomous:Unknown']+outline 
                    else:
                        if "GAG" in subb_splited :
                            gag_ing=subb_splited.index('GAG')
                             #get autonoumos gypsy 
                            if ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "RT"==subb_splited[gag_ing+2]
                                and "RH"==subb_splited[gag_ing+3] and "INT"==subb_splited[gag_ing+4]):
                                    #append to list
                                    outline=[ f'Nonautonomous:{clade}']+outline
                            elif ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "INT"==subb_splited[gag_ing+2]
                                and "RT"==subb_splited[gag_ing+3] and "RH"==subb_splited[gag_ing+4]):
                                     #append to list
                                     outline=[ f'Nonautonomous:{clade}']+outline        
                            else:
                                # gag in domains but some domaine are missing
                                 outline=[ f'Nonautonomous:{clade}']+outline
                        # gag is not in domaine and prot is
                        else :
                                outline=[ f'Nonautonomous:{clade}']+outline
                
                
                #third part of the code------------------------------------------------------------------
                elif   pbs !='' and ppt =='':
                    
                    if len(subb_splited[:-1])== 0:
                            
                            #not trim nor a lard and no domaine exsist
                            outline=[ f'Nonautonomous:Unknown']+outline 
                    else:
                        
                        if "GAG" in subb_splited :
                            
                            gag_ing=subb_splited.index('GAG')
                             #get autonoumos gypsy 
                            if ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "RT"==subb_splited[gag_ing+2]
                                and "RH"==subb_splited[gag_ing+3] and "INT"==subb_splited[gag_ing+4]):
                                    #append to list
                                    
                                    outline=[ f'Nonautonomous:{clade}']+outline
                            #get autonoumos copia         
                            elif ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "INT"==subb_splited[gag_ing+2]
                                and "RT"==subb_splited[gag_ing+3] and "RH"==subb_splited[gag_ing+4]):
                                     #append to list
                                    
                                    outline=[ f'Nonautonomous:{clade}']+outline        
                            else:
                               
                                # gag in domains but some domaine are missing
                                outline=[ f'Nonautonomous:{clade}']+outline
                        # gag is not in domaine and prot is
                        else: 
                                outline=[ f'Nonautonomous:{clade}']+outline
                       
                #third part of the code------------------------------------------------------------------
                #pbs =='' and ppt =='': 
                else   :  
                    
                    if len(subb_splited[:-1])== 0:
                            #not trim nor a lard and no domaine exsist
                             outline=[ f'Nonautonomous:{clade}']+outline 
                    else:
                        if "GAG" in subb_splited :
                            gag_ing=subb_splited.index('GAG')
                             #get autonoumos gypsy 
                            if ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "RT"==subb_splited[gag_ing+2]
                                and "RH"==subb_splited[gag_ing+3] and "INT"==subb_splited[gag_ing+4]):
                                    #append to list
                                    outline=[ f'Nonautonomous:{clade}']+outline
                            #get autonoumos copia         
                            elif ("GAG"==subb_splited[gag_ing] and "PROT"==subb_splited[gag_ing+1] and "INT"==subb_splited[gag_ing+2]
                                and "RT"==subb_splited[gag_ing+3] and "RH"==subb_splited[gag_ing+4]):
                                    #append to list
                                    outline=[ f'Nonautonomous:{clade}']+outline        
                            else:
                                # gag in domains but some domaine are missing
                                 outline=[ f'Nonautonomous:{clade}']+outline
                        # gag is not in domaine and prot is
                        else :
                                outline=[ f'Nonautonomous:{clade}']+outline
               
                out.write("\t".join(outline))       
                                             
