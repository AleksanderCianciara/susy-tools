from itertools import combinations, cycle
import sys
import numpy as np

class InterCorrelators:
	def __init__(self, octet1, octet2):
		self.multiplets = [octet1,octet2]
		if octet1 == 'RA':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["12345678","21436587","34127856","43218765","56781234","65872143","78563412","87654321"]

		elif octet1 == 'CC':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14235867","23146758","32417685","41328576","58671423","67582314","76853241","85764132"]

		elif octet1 == 'CT':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14235786","23146875","32417568","41328657","58671342","67582431","76853124","85764213"]

		elif octet1 == 'CV':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["14236857","23145768","32418675","41327586","58672413","67581324","76854231","85763142"]

		elif octet1 == 'TT':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["13425786","24316875","31247568","42138657","57861342","68752431","75683124","86574213"]
		
		elif octet1 == 'TV':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["13426857","24315768","31248675","42137586","57862413","68751324","75684231","86573142"]

		elif octet1 == 'VV':
			self.p1,self.p2,self.p3,self.p4,self.p5,self.p6,self.p7,self.p8 = ["24136857","13245768","42318675","31427586","68572413","57681324","86754231","75863142"]	

		if octet2 == 'RA':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["12345678","21436587","34127856","43218765","56781234","65872143","78563412","87654321"]

		elif octet2 == 'CC':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["14235867","23146758","32417685","41328576","58671423","67582314","76853241","85764132"]

		elif octet2 == 'CT':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["14235786","23146875","32417568","41328657","58671342","67582431","76853124","85764213"]

		elif octet2 == 'CV':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["14236857","23145768","32418675","41327586","58672413","67581324","76854231","85763142"]

		elif octet2 == 'TT':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["13425786","24316875","31247568","42138657","57861342","68752431","75683124","86574213"]
		
		elif octet2 == 'TV':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["13426857","24315768","31248675","42137586","57862413","68751324","75684231","86573142"]

		elif octet2 == 'VV':
			self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8 = ["24136857","13245768","42318675","31427586","68572413","57681324","86754231","75863142"]	



	def kendall_tau_distance(self, starting, ending):
		pairs = combinations(range(1, len(starting)+1), 2)
		distance = 0
		for x, y in pairs:
			a = starting.index(str(x)) - starting.index(str(y))
			b = ending.index(str(x)) - ending.index(str(y))
			if a * b < 0:
				distance += 1
		return distance

	def inter_correlator_calculation(self):

		p1_q1 = str(self.kendall_tau_distance(self.p1,self.q1))
		p1_q2 = str(self.kendall_tau_distance(self.p1,self.q2))
		p1_q3 = str(self.kendall_tau_distance(self.p1,self.q3))
		p1_q4 = str(self.kendall_tau_distance(self.p1,self.q4))
		p1_q5 = str(self.kendall_tau_distance(self.p1,self.q5))
		p1_q6 = str(self.kendall_tau_distance(self.p1,self.q6))
		p1_q7 = str(self.kendall_tau_distance(self.p1,self.q7))
		p1_q8 = str(self.kendall_tau_distance(self.p1,self.q8))


		p2_q1 = str(self.kendall_tau_distance(self.p2,self.q1))
		p2_q2 = str(self.kendall_tau_distance(self.p2,self.q2))
		p2_q3 = str(self.kendall_tau_distance(self.p2,self.q3))
		p2_q4 = str(self.kendall_tau_distance(self.p2,self.q4))
		p2_q5 = str(self.kendall_tau_distance(self.p2,self.q5))
		p2_q6 = str(self.kendall_tau_distance(self.p2,self.q6))
		p2_q7 = str(self.kendall_tau_distance(self.p2,self.q7))
		p2_q8 = str(self.kendall_tau_distance(self.p2,self.q8))


		p3_q1 = str(self.kendall_tau_distance(self.p3,self.q1))
		p3_q2 = str(self.kendall_tau_distance(self.p3,self.q2))
		p3_q3 = str(self.kendall_tau_distance(self.p3,self.q3))
		p3_q4 = str(self.kendall_tau_distance(self.p3,self.q4))
		p3_q5 = str(self.kendall_tau_distance(self.p3,self.q5))
		p3_q6 = str(self.kendall_tau_distance(self.p3,self.q6))
		p3_q7 = str(self.kendall_tau_distance(self.p3,self.q7))
		p3_q8 = str(self.kendall_tau_distance(self.p3,self.q8))


			
		p4_q1 = str(self.kendall_tau_distance(self.p4,self.q1))
		p4_q2 = str(self.kendall_tau_distance(self.p4,self.q2))
		p4_q3 = str(self.kendall_tau_distance(self.p4,self.q3))
		p4_q4 = str(self.kendall_tau_distance(self.p4,self.q4))
		p4_q5 = str(self.kendall_tau_distance(self.p4,self.q5))
		p4_q6 = str(self.kendall_tau_distance(self.p4,self.q6))
		p4_q7 = str(self.kendall_tau_distance(self.p4,self.q7))
		p4_q8 = str(self.kendall_tau_distance(self.p4,self.q8))

		p5_q1 = str(self.kendall_tau_distance(self.p5,self.q1))
		p5_q2 = str(self.kendall_tau_distance(self.p5,self.q2))
		p5_q3 = str(self.kendall_tau_distance(self.p5,self.q3))
		p5_q4 = str(self.kendall_tau_distance(self.p5,self.q4))
		p5_q5 = str(self.kendall_tau_distance(self.p5,self.q5))
		p5_q6 = str(self.kendall_tau_distance(self.p5,self.q6))
		p5_q7 = str(self.kendall_tau_distance(self.p5,self.q7))
		p5_q8 = str(self.kendall_tau_distance(self.p5,self.q8))


		p6_q1 = str(self.kendall_tau_distance(self.p6,self.q1))
		p6_q2 = str(self.kendall_tau_distance(self.p6,self.q2))
		p6_q3 = str(self.kendall_tau_distance(self.p6,self.q3))
		p6_q4 = str(self.kendall_tau_distance(self.p6,self.q4))
		p6_q5 = str(self.kendall_tau_distance(self.p6,self.q5))
		p6_q6 = str(self.kendall_tau_distance(self.p6,self.q6))
		p6_q7 = str(self.kendall_tau_distance(self.p6,self.q7))
		p6_q8 = str(self.kendall_tau_distance(self.p6,self.q8))


		p7_q1 = str(self.kendall_tau_distance(self.p7,self.q1))
		p7_q2 = str(self.kendall_tau_distance(self.p7,self.q2))
		p7_q3 = str(self.kendall_tau_distance(self.p7,self.q3))
		p7_q4 = str(self.kendall_tau_distance(self.p7,self.q4))
		p7_q5 = str(self.kendall_tau_distance(self.p7,self.q5))
		p7_q6 = str(self.kendall_tau_distance(self.p7,self.q6))
		p7_q7 = str(self.kendall_tau_distance(self.p7,self.q7))
		p7_q8 = str(self.kendall_tau_distance(self.p7,self.q8))


		p8_q1 = str(self.kendall_tau_distance(self.p8,self.q1))
		p8_q2 = str(self.kendall_tau_distance(self.p8,self.q2))
		p8_q3 = str(self.kendall_tau_distance(self.p8,self.q3))
		p8_q4 = str(self.kendall_tau_distance(self.p8,self.q4))
		p8_q5 = str(self.kendall_tau_distance(self.p8,self.q5))
		p8_q6 = str(self.kendall_tau_distance(self.p8,self.q6))
		p8_q7 = str(self.kendall_tau_distance(self.p8,self.q7))
		p8_q8 = str(self.kendall_tau_distance(self.p8,self.q8))

		print("for " + str(self.multiplets[0]) + "-" + str(self.multiplets[1]) + " we have the intercorrelator matrix ")

		arr = np.array([[p1_q1,p1_q2,p1_q3,p1_q4,p1_q5,p1_q6,p1_q7,p1_q8],
						[p2_q1,p2_q2,p2_q3,p2_q4,p2_q5,p2_q6,p2_q7,p2_q8],
						[p3_q1,p3_q2,p3_q3,p3_q4,p3_q5,p3_q6,p3_q7,p3_q8],
						[p4_q1,p4_q2,p4_q3,p4_q4,p4_q5,p4_q6,p4_q7,p4_q8],
						[p5_q1,p5_q2,p5_q3,p5_q4,p5_q5,p5_q6,p5_q7,p5_q8],
						[p6_q1,p6_q2,p6_q3,p6_q4,p6_q5,p6_q6,p6_q7,p6_q8],
						[p7_q1,p7_q2,p7_q3,p7_q4,p7_q5,p7_q6,p7_q7,p7_q8],
						[p8_q1,p8_q2,p8_q3,p8_q4,p8_q5,p8_q6,p8_q7,p8_q8]
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
	oct1 = raw_input('Choose starting Octet ')
	oct2 = raw_input('Choose ending Octet ')
	IC = InterCorrelators(octet1=oct1,octet2=oct2)
	IC.inter_correlator_calculation()

	


