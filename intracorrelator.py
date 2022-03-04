from itertools import combinations, cycle
import sys
import numpy as np


class IntraCorrelators:
	def __init__(self,octet):
		if octet == 'RA':
			self.multiplet = 'Rana multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["12345678","21436587","34127856","43218765","56781234","65872143","78563412","87654321"]

		if octet == 'CC':
			self.multiplet = 'Chiral-Chiral multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14235867","23146758","32417685","41328576","58671423","67582314","76853241","85764132"]

		if octet == 'CT':
			self.multiplet = 'Chiral-Tensor multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14235786","23146875","32417568","41328657","58671342","67582431","76853124","85764213"]

		if octet == 'CV':
			self.multiplet = 'Chiral-Vector multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14236857","23145768","32418675","41327586","58672413","67581324","76854231","85763142"]


		if octet == 'TT':
			self.multiplet = 'Tensor-Tensor multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["13425786","24316875","31247568","42138657","57861342","68752431","75683124","86574213"]
		
		if octet == 'TV':
			self.multiplet = 'Tensor-Vector multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["13426857","24315768","31248675","42137586","57862413","68751324","75684231","86573142"]

		if octet == 'VV':
			self.multiplet = 'Vector-Vector multiplet'
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["24136857","13245768","42318675","31427586","68572413","57681324","86754231","75863142"]



		self.indices = self.index_combinations()

	def index_combinations(self):
		comb = combinations([1,2,3,4,5,6,7,8],2)
		combs = []
		for i in comb:
			combs.append(i)
		inds = [list(el) for el in combs]
		return inds

	def kendall_tau_distance(self, starting, ending):
		pairs = combinations(range(1, len(starting)+1), 2)
		distance = 0
		for x, y in pairs:
			a = starting.index(str(x)) - starting.index(str(y))
			b = ending.index(str(x)) - ending.index(str(y))
			if a * b < 0:
				distance += 1
		return distance

	def intra_correlator_calculation(self):
		p1_p1=p2_p2=p3_p3=p4_p4=p5_p5=p6_p6=p7_p7=p8_p8='0'

		p1_p2 = str(self.kendall_tau_distance(self.p1,self.p2))
		p1_p3 = str(self.kendall_tau_distance(self.p1,self.p3))
		p1_p4 = str(self.kendall_tau_distance(self.p1,self.p4))
		p1_p5 = str(self.kendall_tau_distance(self.p1,self.p5))
		p1_p6 = str(self.kendall_tau_distance(self.p1,self.p6))
		p1_p7 = str(self.kendall_tau_distance(self.p1,self.p7))
		p1_p8 = str(self.kendall_tau_distance(self.p1,self.p8))

		print('for the ' + self.multiplet + ' we find:')

		print('distances from p1 are (in order) ' + str([p1_p1,p1_p2,p1_p3,p1_p4,p1_p5,p1_p6,p1_p7,p1_p8])
		+ ' for a sum of ' + str(sum([int(p1_p2),int(p1_p3),int(p1_p4),int(p1_p5),int(p1_p6),int(p1_p7),int(p1_p8)])))

		p2_p1 = str(self.kendall_tau_distance(self.p2,self.p1))
		p2_p3 = str(self.kendall_tau_distance(self.p2,self.p3))
		p2_p4 = str(self.kendall_tau_distance(self.p2,self.p4))
		p2_p5 = str(self.kendall_tau_distance(self.p2,self.p5))
		p2_p6 = str(self.kendall_tau_distance(self.p2,self.p6))
		p2_p7 = str(self.kendall_tau_distance(self.p2,self.p7))
		p2_p8 = str(self.kendall_tau_distance(self.p2,self.p8))

		print('distances from p2 are (in order) ' + str([p2_p1,p2_p2,p2_p3,p2_p4,p2_p5,p2_p6,p2_p7,p2_p8])
		+ ' for a sum of ' + str(sum([int(p2_p1),int(p2_p3),int(p2_p4),int(p2_p5),int(p2_p6),int(p2_p7),int(p2_p8)])))

		p3_p1 = str(self.kendall_tau_distance(self.p3,self.p1))
		p3_p2 = str(self.kendall_tau_distance(self.p3,self.p2))
		p3_p4 = str(self.kendall_tau_distance(self.p3,self.p4))
		p3_p5 = str(self.kendall_tau_distance(self.p3,self.p5))
		p3_p6 = str(self.kendall_tau_distance(self.p3,self.p6))
		p3_p7 = str(self.kendall_tau_distance(self.p3,self.p7))
		p3_p8 = str(self.kendall_tau_distance(self.p3,self.p8))

		print('distances from p3 are (in order) ' + str([p3_p1,p3_p2,p3_p3,p3_p4,p3_p5,p3_p6,p3_p7,p3_p8])
		+ ' for a sum of ' + str(sum([int(p3_p1),int(p3_p2),int(p3_p4),int(p3_p5),int(p3_p6),int(p3_p7),int(p3_p8)])))
			
		p4_p1 = str(self.kendall_tau_distance(self.p4,self.p1))
		p4_p2 = str(self.kendall_tau_distance(self.p4,self.p2))
		p4_p3 = str(self.kendall_tau_distance(self.p4,self.p3))
		p4_p5 = str(self.kendall_tau_distance(self.p4,self.p5))
		p4_p6 = str(self.kendall_tau_distance(self.p4,self.p6))
		p4_p7 = str(self.kendall_tau_distance(self.p4,self.p7))
		p4_p8 = str(self.kendall_tau_distance(self.p4,self.p8))

		print('distances from p4 are (in order) ' + str([p4_p1,p4_p2,p4_p3,p4_p4,p4_p5,p4_p6,p4_p7,p4_p8])
		+ ' for a sum of ' + str(sum([int(p4_p1),int(p4_p2),int(p4_p3),int(p4_p5),int(p4_p6),int(p4_p7),int(p4_p8)])))

		p5_p1 = str(self.kendall_tau_distance(self.p5,self.p1))
		p5_p2 = str(self.kendall_tau_distance(self.p5,self.p2))
		p5_p3 = str(self.kendall_tau_distance(self.p5,self.p3))
		p5_p4 = str(self.kendall_tau_distance(self.p5,self.p4))
		p5_p6 = str(self.kendall_tau_distance(self.p5,self.p6))
		p5_p7 = str(self.kendall_tau_distance(self.p5,self.p7))
		p5_p8 = str(self.kendall_tau_distance(self.p5,self.p8))

		print('distances from p5 are (in order) ' + str([p5_p1,p5_p2,p5_p3,p5_p4,p5_p5,p5_p6,p5_p7,p5_p8])
		+ ' for a sum of ' + str(sum([int(p5_p1),int(p5_p2),int(p5_p3),int(p5_p4),int(p5_p6),int(p5_p7),int(p5_p8)])))

		p6_p1 = str(self.kendall_tau_distance(self.p6,self.p1))
		p6_p2 = str(self.kendall_tau_distance(self.p6,self.p2))
		p6_p3 = str(self.kendall_tau_distance(self.p6,self.p3))
		p6_p4 = str(self.kendall_tau_distance(self.p6,self.p4))
		p6_p5 = str(self.kendall_tau_distance(self.p6,self.p5))
		p6_p7 = str(self.kendall_tau_distance(self.p6,self.p7))
		p6_p8 = str(self.kendall_tau_distance(self.p6,self.p8))

		print('distances from p6 are (in order) ' + str([p6_p1,p6_p2,p6_p3,p6_p4,p6_p5,p6_p6,p6_p7,p6_p8])
		+ ' for a sum of ' + str(sum([int(p6_p1),int(p6_p2),int(p6_p3),int(p6_p4),int(p6_p5),int(p6_p7),int(p6_p8)])))

		p7_p1 = str(self.kendall_tau_distance(self.p7,self.p1))
		p7_p2 = str(self.kendall_tau_distance(self.p7,self.p2))
		p7_p3 = str(self.kendall_tau_distance(self.p7,self.p3))
		p7_p4 = str(self.kendall_tau_distance(self.p7,self.p4))
		p7_p5 = str(self.kendall_tau_distance(self.p7,self.p5))
		p7_p6 = str(self.kendall_tau_distance(self.p7,self.p6))
		p7_p8 = str(self.kendall_tau_distance(self.p7,self.p8))

		print('distances from p7 are (in order) ' + str([p7_p1,p7_p2,p7_p3,p7_p4,p7_p5,p7_p6,p7_p7,p7_p8])
		+ ' for a sum of ' + str(sum([int(p7_p1),int(p7_p2),int(p7_p3),int(p7_p4),int(p7_p5),int(p7_p6),int(p7_p8)])))

		p8_p1 = str(self.kendall_tau_distance(self.p8,self.p1))
		p8_p2 = str(self.kendall_tau_distance(self.p8,self.p2))
		p8_p3 = str(self.kendall_tau_distance(self.p8,self.p3))
		p8_p4 = str(self.kendall_tau_distance(self.p8,self.p4))
		p8_p5 = str(self.kendall_tau_distance(self.p8,self.p5))
		p8_p6 = str(self.kendall_tau_distance(self.p8,self.p6))
		p8_p7 = str(self.kendall_tau_distance(self.p8,self.p7))

		print('distances from p8 are (in order) ' + str([p8_p1,p8_p2,p8_p3,p8_p4,p8_p5,p8_p6,p8_p7,p8_p8])
		+ ' for a sum of ' + str(sum([int(p8_p1),int(p8_p2),int(p8_p3),int(p8_p4),int(p8_p5),int(p8_p6),int(p8_p7)])))


		arr = np.array([[p1_p1,p1_p2,p1_p3,p1_p4,p1_p5,p1_p6,p1_p7,p1_p8],
						[p2_p1,p2_p2,p2_p3,p2_p4,p2_p5,p2_p6,p2_p7,p2_p8],
						[p3_p1,p3_p2,p3_p3,p3_p4,p3_p5,p3_p6,p3_p7,p3_p8],
						[p4_p1,p4_p2,p4_p3,p4_p4,p4_p5,p4_p6,p4_p7,p4_p8],
						[p5_p1,p5_p2,p5_p3,p5_p4,p5_p5,p5_p6,p5_p7,p5_p8],
						[p6_p1,p6_p2,p6_p3,p6_p4,p6_p5,p6_p6,p6_p7,p6_p8],
						[p7_p1,p7_p2,p7_p3,p7_p4,p7_p5,p7_p6,p7_p7,p7_p8],
						[p8_p1,p8_p2,p8_p3,p8_p4,p8_p5,p8_p6,p8_p7,p8_p8]
						],dtype=int)
		print('\n')
		print('matrix of correlators is ')
		print(arr)
		print('\n')
		print('with eigenvalues ')
		print(np.round(np.real(np.linalg.eigvals(arr))))
		print('and sum of eigenvalues equal to ')
		print(np.round(sum(np.real(np.linalg.eigvals(arr)))))



if __name__=="__main__":
	cont = True
	while cont:
		raw = raw_input('Choose from Chiral-Chiral (CC), Rana basis (RA), Chiral-Vector (CV), Chiral-Tensor (CT), Vector-Vector (VV), Tensor-Vector (TV), or Tensor-Tensor (TT) multiplets: ')
		IC = IntraCorrelators(raw)
		IC.intra_correlator_calculation()
		raw = raw_input('Do you wish to calculate another (Y/N)? ')
		if raw.lower() == "y":
			pass
		if raw.lower() == "n":
			cont = False
	


