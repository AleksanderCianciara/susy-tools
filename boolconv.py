

def bool(number):
	nmbrs = list(number)
	lst = list(map(int,nmbrs))
	count = 0
	for idx,val in enumerate(nmbrs):
		count += int(val)*2**idx
	print(count)

	# if len(nmbrs)==8:
	# 	x0,x1,x2,x3,x4,x5,x6,x7 = list(map(int,nmbrs))
	# 	print((x0)*2**0+(x1)*2**1+(x2)*2**2+(x3)*2**3+(x4)*2**4+(x5)*2**5+(x6)*2**6+(x7)*2**7)
	# if len(nmbrs)==4:
	# 	x0,x1,x2,x3 = list(map(int,nmbrs))
	# 	print((x0)*2**0+(x1)*2**1+(x2)*2**2+(x3)*2**3)	




bool("01010101")



# bool("10110101")

# 0011
# 0010

# bool("00110010")
# bool("00100011")