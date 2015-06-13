import random

def sequentialy(qset):
    rset = []
    for i in range(0, len(qset)):
	rset.append(qset[i])
    
    return rset


def randomly(_qset):
    qset = list(_qset)
    rset = []
    while qset:
        sel = random.choice(qset)
	rset.append(sel)
	qset.remove(sel)

    return rset

def fullBFS(qset):
    if len(qset) <= 2:
        return qset

    rset = []
    rset.append(qset[0])
    rset.append(qset[-1])
    end = len(qset) - 2

    queue = []
    cur = Node((1 + end + 1) / 2, 1 ,end)
    queue.append(cur)

    while len(queue) != 0 :
        cur = queue.pop(0)
        rset.append(qset[cur.val])

        if cur.left < cur.val:
            val = (cur.left + cur.val) / 2
            leftSon = Node(val, cur.left, cur.val - 1)
            queue.append(leftSon)

        if cur.val < cur.right:
            val = (cur.right + cur.val +1) / 2
            rightSon = Node(val, cur.val+1, cur.right)
            queue.append(rightSon)

    return rset


##################### below is abandoned #######################

def binary(qset, num):
    if num > len(qset) or num <= 0:
	print "subquery amount overflow! num = " + str(num)
	return qset

    rset = []
    rset.append(qset[0])
    if num == 1:
	return rset
	
    end = len(qset) - 1
    rset.append(qset[end])
    count = []
    count.append(num - 2) # minus 2: head and tail element
    bs_helper(qset, rset, 1, end-1, count)
    return rset
    
def bs_helper(qset, rset, left, right, count):
    if left > right or count[0] <= 0:
	return
    mid = (left + right)/2
#    print qset[mid]
    rset.append(qset[mid])
    count[0] -= 1
    bs_helper(qset, rset, left, mid - 1, count)
    bs_helper(qset, rset, mid + 1, right, count)

def skip(qset, num):
    if num > len(qset) or num <= 0:
	print "subquery amount overflow! num = " + str(num)
	return qset
    
    visited = [False] * len(qset)
    rset = []
    rset.append(qset[0])
    visited[0] = True
    if num == 1:
	return rset

    end = len(qset) - 1
    rset.append(qset[end])
    visited[end] = True
    num -= 2
    step = len(qset) / 2
    while num > 0 and step > 0:
	cur = 0
	while cur < end:
	    if not visited[cur]:
		rset.append(qset[cur])
		visited[cur] = True
	    cur += step
	
	step /= 2

    return rset

class Node(object):
    def __init__(self, val, left, right):
	self.val = val
	self.left = left
	self.right = right


def BTLT(qset, num):
    if num >= len(qset) or num <= 0:
#	print "subquery amount overflow! num = " + str(num)
	return qset
 
    rset = []
    rset.append(qset[0])
    if num == 1:
	return rset
	
    end = len(qset) - 1
    rset.append(qset[end])
    if num == 2:
	return rset

    visited = [False] * len(qset)
    visited[0] = True
    visited[1] = True
    num -= 2
    queue = []
    cur = Node(len(qset) / 2, 0 ,end)
    queue.append(cur)
    while len(queue) != 0 and num > 0:
	cur = queue.pop(0)
	if not visited[cur.val]:
            rset.append(qset[cur.val])
	    visited[cur.val] = True
	    num -= 1

	leftVal = (cur.left + cur.val) / 2 
	leftSon = Node(leftVal, cur.left, cur.val)
	rightVal = (cur.right + cur.val) / 2
	rightSon = Node(rightVal, cur.val, cur.right)
	queue.append(leftSon)
	queue.append(rightSon)

    return rset

