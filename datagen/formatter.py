s = '''
{                    0,  2./63 + 4*sqrt(3)/63,                    0 },
{    sqrt(3)/63 + 2./21, -2*sqrt(3)/63 - 1./63,                    0 },
{   -2./21 - sqrt(3)/63, -2*sqrt(3)/63 - 1./63,                    0 },
{   -2./21 - sqrt(3)/63,  1./21 + 2*sqrt(3)/21,                    0 },
{    sqrt(3)/63 + 2./21,  1./21 + 2*sqrt(3)/21,                    0 },
{  2*sqrt(3)/63 + 4./21,                    0,                    0 },
{    sqrt(3)/63 + 2./21, -2*sqrt(3)/21 - 1./21,                    0 },
{   -2./21 - sqrt(3)/63, -2*sqrt(3)/21 - 1./21,                    0 },
{ -4./21 - 2*sqrt(3)/63,                    0,                    0 },
{                    0,                    0,                    0 },
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
