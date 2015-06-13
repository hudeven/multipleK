import re
#seq = 'TTTT' # pattern 0
#seq = 'ATGCTGCTGCTGCA'  # pattern 1
#seq = 'ATGCTGCTGCTGCT'  # pattern 2
#seq = 'CTGCTGCTGCTGCTGCTGCCGGCG'  # pattern 3 
#seq = 'AGTAATAGAATAGCCGTA' # pattern 4
def getPattern(seq):
  patternId = 5;
  # pattern 0:  AAAAA
  match = re.search(r"^([A,T,G,C]?)\1+$", seq)
  if match:
    #print "pattern 0: " + match.group(1)
    patternId = 4 
    return patternId
  # pattern 1 or 2 or 3:
  match = re.search(r"^[A,T,G,C]([A,T,G,C]+?)\1+([A,T,G,C]*)$", seq)
  if match:
    pattern = match.group(1)
    tail = match.group(2)
    #print "pattern = " + pattern
    #print "tail = " + tail

    if len(tail) < len(pattern):
      #print tail[0:-1] + " = " + pattern[0:len(tail)-1]
      if tail[0:-1] == pattern[0:len(tail)-1]:
        if len(tail) == 0 or tail[-1] == pattern[len(tail)-1]:
          if seq[0] == pattern[-1]:
            patternId = 3 
          else:
            patternId = 2 
        else:
            patternId = 1 
    
  #print "seq = "  + seq
  #print "pattern id = " + str(patternId)
  return patternId


