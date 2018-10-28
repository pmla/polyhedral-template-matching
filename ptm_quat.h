#ifndef PTM_QUAT_H
#define PTM_QUAT_H

namespace ptm {

int rotate_quaternion_into_cubic_fundamental_zone(double* q);
int rotate_quaternion_into_diamond_cubic_fundamental_zone(double* q);
int rotate_quaternion_into_icosahedral_fundamental_zone(double* q);
int rotate_quaternion_into_hcp_fundamental_zone(double* q);
int rotate_quaternion_into_hcp_crystalline_fundamental_zone(double* q);
int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double* q);

void quaternion_to_rotation_matrix(double* q, double* U);

}

#endif

