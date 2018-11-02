s = '''
{  1,   0,   0,   0},
{  sqrt(3)/2,   0,   0,   0.5},
{  sqrt(3)/2,   0,   0,  -0.5},
{  0.5,   0,   0,   sqrt(3)/2},
{  0.5,  -0,  -0,  -sqrt(3)/2},
{  0,   1,   0,   0},
{  0,   sqrt(3)/2,   0.5,   0},
{  0,   sqrt(3)/2,  -0.5,   0},
{  0,   0.5,   sqrt(3)/2,   0},
{ -0,   0.5,  -sqrt(3)/2,  -0},
{  0,   0,   1,   0},
{  0,   0,   0,   1},
'''

def go():
	data = []
	lines = s.split('\n')
	lines = [e for e in lines if len(e) > 1]
	for line in lines:
		line = line.replace('{', '')
		line = line.replace('}', '')
		line = line.replace(' ', '')
		line = line.replace('\t', '')
		line = line.split(',')[:-1]
		data += [line]

	nmax = max([max([len(e) for e in row]) for row in data])
	print nmax

	output = []
	for row in data:
		print "{ " + ", ".join([e.rjust(nmax) for e in row]) + " },"
go()
