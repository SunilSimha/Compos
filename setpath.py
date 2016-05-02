import os
import shutil

compospath = os.getcwd() #get the path of compos packages

cambpath = '/home/yanziang/research/software/camb' #path of camb, set by user

con = open('./compos/const.py','r')

text = con.read()

begintext1 = text.find('\'path\' : ')
endtext1 = text.find(',#path of compos')
cop = text[begintext1:endtext1]

begintext2 = text.find('\'cambpath\' : ')
endtext2 = text.find('#path of camb')
cap = text[begintext2:endtext2]

text = text.replace(cop,'\'path\' : \'' + compospath + '\'')
if cambpath != '':
    text = text.replace(cap,'\'cambpath\' : \''+ cambpath + '\'')
con.close()

con = open('./compos/const.py','w')
con.write(text)
con.close

if cambpath != '\'\'':
    shutil.copyfile('paramsformps.ini', cambpath + 'paramsformps.ini')
