s = '''
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{                                   0,                                   0,       4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,                                   0,                                   0 },

{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{                                   0,                                   0,      -4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,                                   0,                                   0 },

{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
{                                   0,                                   0,      -4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                     -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),                      4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
{                                   0,                                   0,                                   0 },
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
