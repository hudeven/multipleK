import os
import sys

start_klength = 5
delta = 1
end_klength = 6

for i in range(start_klength, end_klength):
	cmd = "python multiplek.py --repeat 1 --klength " + str(i) + " --multiplek " + str(i+delta) + " > result/output"+ str(i) +"_"+ str(i+delta)
	os.system(cmd)	
