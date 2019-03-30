#ifndef PTM_LINEAR_ASSIGNMENT
#define PTM_LINEAR_ASSIGNMENT

namespace ptm {

double linear_assignment_sum(	const std::vector < std::vector < double > >& cost,
				std::vector < int >& lmate, std::vector < int >& rmate);
}

#endif

