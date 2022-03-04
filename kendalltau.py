from itertools import combinations, cycle
import sys

def kendall_tau_distance(starting, ending):
    pairs = combinations(range(1, len(starting)+1), 2)
    distance = 0
    for x, y in pairs:
        a = starting.index(x) - starting.index(y)
        b = ending.index(x) - ending.index(y)
        if a * b < 0:
            distance += 1
    return distance

print(kendall_tau_distance([1,4,2,3,5,8,6,7,9,12,10,11,13,16,14,15], [2,3,1,4,6,7,5,8,10,11,9,12,14,15,13,16]))

