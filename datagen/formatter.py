s = '''
{ 0,  0,  0},
{ 11/(2*(sqrt(3) + 3 + 6*sqrt(2))),  11/(2*(sqrt(3) + 3 + 6*sqrt(2))),  11/(2*(sqrt(3) + 3 + 6*sqrt(2)))},
{-11/(2*(sqrt(3) + 3 + 6*sqrt(2))),  11/(2*(sqrt(3) + 3 + 6*sqrt(2))), -11/(2*(sqrt(3) + 3 + 6*sqrt(2)))},
{ 11/(2*(sqrt(3) + 3 + 6*sqrt(2))), -11/(2*(sqrt(3) + 3 + 6*sqrt(2))), -11/(2*(sqrt(3) + 3 + 6*sqrt(2)))},
{-11/(2*(sqrt(3) + 3 + 6*sqrt(2))), -11/(2*(sqrt(3) + 3 + 6*sqrt(2))),  11/(2*(sqrt(3) + 3 + 6*sqrt(2)))},
{ 0,  0, -11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 0, -11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
{-11/(sqrt(3) + 3 + 6*sqrt(2)),  0,  0},
{ 0,  0,  11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 0,  11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
{ 11/(sqrt(3) + 3 + 6*sqrt(2)),  0,  0},
{ 0, -11/(sqrt(3) + 3 + 6*sqrt(2)), -11/(sqrt(3) + 3 + 6*sqrt(2))},
{-11/(sqrt(3) + 3 + 6*sqrt(2)),  0, -11/(sqrt(3) + 3 + 6*sqrt(2))},
{-11/(sqrt(3) + 3 + 6*sqrt(2)), -11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
{-11/(sqrt(3) + 3 + 6*sqrt(2)),  11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
{ 11/(sqrt(3) + 3 + 6*sqrt(2)), -11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
{-11/(sqrt(3) + 3 + 6*sqrt(2)),  0,  11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 11/(sqrt(3) + 3 + 6*sqrt(2)),  0, -11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 0, -11/(sqrt(3) + 3 + 6*sqrt(2)),  11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 0,  11/(sqrt(3) + 3 + 6*sqrt(2)), -11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 0,  11/(sqrt(3) + 3 + 6*sqrt(2)),  11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 11/(sqrt(3) + 3 + 6*sqrt(2)),  0,  11/(sqrt(3) + 3 + 6*sqrt(2))},
{ 11/(sqrt(3) + 3 + 6*sqrt(2)),  11/(sqrt(3) + 3 + 6*sqrt(2)),  0},
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
