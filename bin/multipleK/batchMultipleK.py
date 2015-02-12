import os
import sys

start_klength = 6
delta = 0
end_klength = 7

for i in range(start_klength, end_klength):
	cmd = "python multiplek.py --repeat 1 --klength " + str(i) + " --multiplek " + str(i+delta) + " > result/output"+ str(i) +"_"+ str(i+delta)
	os.system(cmd)	
