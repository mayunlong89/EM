#!usr/bin/env python

# -*- coding: utf-8 -*-
"""
@author: Yunlong Ma
@E-mail:glb-biotech@zju.edu.cn

This script was desiged to combine all the expression files downloaded from the ArrayExpress database.

"""

#import packages
import os 
import linecache  

#set the work directory
#all downloaded files were stored in the test file
root = 'F:\\Desktop\\EM-RNA-data\\test' 
file_names = os.listdir(root) 
file_ob_list = [] 
for file_name in file_names:
    fileob = root + '\\' + file_name 
    file_ob_list.append(fileob)                              

  
ldata = []
data = []

line_num = 1 
total_line = len(open(file_ob_list[0]).readlines())
while line_num <= total_line:  
    for file_ob in file_ob_list:
        line = linecache.getline(file_ob,line_num) 
        line = line.strip() 
        if line is None or len(line) == 0:
            break
        fields =  line.split('\t') 
        prob = fields[1] 
              
        if file_ob == file_ob_list[0]: 
            data = [fields[0],prob] 
        else:                      
            data.append(prob)                            
    line_num = line_num +1 
    ldata.append(data) 
    data = []          

ldata.pop(0) 


#add header for data
file_names_new1 = []
for col_name in file_names:
    col_name = col_name.strip(".txt")
    file_names_new1.append(col_name)
file_names_new2 = ["Gene",] + file_names_new1
temp = [file_names_new2,] + ldata

#write data out
f = open("F:\\Desktop\EM-RNA-data\\test_combined.txt", "w+") 
for i, p in enumerate(temp):
    for j, q in enumerate(p): 
        f.write(q+"\t") 
    print (i) 
    f.write("\n") 

f.close() 


#End
