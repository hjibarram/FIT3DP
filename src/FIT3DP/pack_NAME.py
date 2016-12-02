#! /usr/bin/env python
import sys

def mycall(comand, alias_py='python'):
    alias_py=''
    from subprocess import call
    line=comand.split(" ")
    if alias_py != '' :
        fcomand=[alias_py]
    else:
        fcomand=[]#2.7']
    fcomand.extend(line)
    linp=''
    nx=len(fcomand)
    for i in range(0, nx):
        linp=linp+fcomand[i]+" "
    print linp+" | PROCESSING"
    call(fcomand)
    print "DONE"
    
sys.argv=filter(None,sys.argv)
if len(sys.argv) < 3:
    print "USE: pack_NAME.py NAME CONF"
    sys.exit(0)

name=sys.argv[1]
conf=sys.argv[2]
type=sys.argv[3]#"v1.5"
call="pack_results_name.py "+name+" "+conf+"/pack_CS.csv SSP"
mycall(call)
call="pack_results_name.py "+name+" "+conf+"/pack_elines_"+type+".csv ELINES"
mycall(call)
call="sum_mass_age.py "+name+"  Mass."+name+".csv"
mycall(call)
call="pack_results_name.py "+name+" "+conf+"/pack_SFH.csv SFH"
mycall(call)