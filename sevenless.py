#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

for int in list(range(0,99)):
    if int % 7 != 0:
        print(int)
