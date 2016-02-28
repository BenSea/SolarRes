def Y2D(a,b,c):
	q = (a*b + a*c + b*c) / c
	r = (a*b + a*c + b*c) / a
	s = (a*b + a*c + b*c) / b
	return q r s
	
def D2Y(a,b,c):
	q = ((b*c)/(a+b+c))
	r = ((a*c)/(a+b+c))
	s = ((a*b)/(a+b+c))
	return q r s
		