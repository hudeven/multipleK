import os
import sys

start_klength = 4
delta = 3
end_klength = 10

for i in range(start_klength, end_klength):
	cmd = "python multiplek.py --repeat 1 --klength " + str(i) + " --multiplek " + str(i+delta) + " > result/output"+ str(i) +"_"+ str(i+delta)
	os.system(cmd)	
