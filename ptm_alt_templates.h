#ifndef PTM_ALT_TEMPLATES_H
#define PTM_ALT_TEMPLATES_H

#include <cmath>


const double ptm_template_hcp_alt1[PTM_NUM_POINTS_HCP][3] = {
								{          0,          0,          0 },
								{          1,          0,          0 },
								{       -0.5, -sqrt(3)/2,          0 },
								{       -0.5, -sqrt(3)/6, -sqrt(6)/3 },
								{          0,  sqrt(3)/3, -sqrt(6)/3 },
								{        0.5, -sqrt(3)/6, -sqrt(6)/3 },
								{         -1,          0,          0 },
								{       -0.5,  sqrt(3)/2,          0 },
								{        0.5,  sqrt(3)/2,          0 },
								{        0.5, -sqrt(3)/2,          0 },
								{        0.5, -sqrt(3)/6,  sqrt(6)/3 },
								{          0,  sqrt(3)/3,  sqrt(6)/3 },
								{       -0.5, -sqrt(3)/6,  sqrt(6)/3 },
};

const double ptm_template_dcub_alt1[PTM_NUM_POINTS_DCUB][3] = {
								{                      0,                      0,                      0 },
								{  4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)) },
								{  4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)) },
								{ -4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)), -4/(sqrt(3)+6*sqrt(2)) },
								{ -4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)),  4/(sqrt(3)+6*sqrt(2)) },
								{  8/(sqrt(3)+6*sqrt(2)),                      0,  8/(sqrt(3)+6*sqrt(2)) },
								{                      0, -8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)) },
								{  8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)),                      0 },
								{                      0,  8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)) },
								{  8/(sqrt(3)+6*sqrt(2)),                      0, -8/(sqrt(3)+6*sqrt(2)) },
								{  8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)),                      0 },
								{ -8/(sqrt(3)+6*sqrt(2)),                      0, -8/(sqrt(3)+6*sqrt(2)) },
								{                      0, -8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)) },
								{ -8/(sqrt(3)+6*sqrt(2)), -8/(sqrt(3)+6*sqrt(2)),                      0 },
								{ -8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)),                      0 },
								{ -8/(sqrt(3)+6*sqrt(2)),                      0,  8/(sqrt(3)+6*sqrt(2)) },
								{                      0,  8/(sqrt(3)+6*sqrt(2)),  8/(sqrt(3)+6*sqrt(2)) },
};

const double ptm_template_dhex_alt1[PTM_NUM_POINTS_DHEX][3] = {
	{                                   0,                                   0,                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),   -4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{                                   0,                                   0,       4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
};

const double ptm_template_dhex_alt2[PTM_NUM_POINTS_DHEX][3] = {
	{                                   0,                                   0,                                   0 },
	{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{                                   0,                                   0,      -4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),   4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{                                   0,   -8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
};

const double ptm_template_dhex_alt3[PTM_NUM_POINTS_DHEX][3] = {
	{                                   0,                                   0,                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),    4*sqrt(3)/(3*sqrt(3)+18*sqrt(2)) },
	{                                   0,                                   0,      -4*sqrt(3)/(sqrt(3)+6*sqrt(2)) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{      -8*sqrt(2)/(sqrt(3)+6*sqrt(2)),                                   0,                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),      -4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)),  16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),       4*sqrt(6)/(sqrt(3)+6*sqrt(2)),                                   0 },
	{      -4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{                                   0,    8*sqrt(6)/(3*sqrt(3)+18*sqrt(2)), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
	{       4*sqrt(2)/(sqrt(3)+6*sqrt(2)),  -4*sqrt(6)/(3*(sqrt(3)+6*sqrt(2))), -16*sqrt(3)/(3*(sqrt(3)+6*sqrt(2))) },
};


const double ptm_template_graphene_alt1[PTM_NUM_POINTS_GRAPHENE][3] = {
									{                    0,                    0,                    0 },
									{   3*sqrt(3)/22-9./11,  -3./22+3*sqrt(3)/11,                    0 },
									{   9./11-3*sqrt(3)/22,  -3./22+3*sqrt(3)/11,                    0 },
									{                    0,  -6*sqrt(3)/11+3./11,                    0 },
									{ -18./11+3*sqrt(3)/11,                    0,                    0 },
									{   3*sqrt(3)/22-9./11,  -9./22+9*sqrt(3)/11,                    0 },
									{   9./11-3*sqrt(3)/22,  -9./22+9*sqrt(3)/11,                    0 },
									{ -3*sqrt(3)/11+18./11,                    0,                    0 },
									{   9./11-3*sqrt(3)/22,  -9*sqrt(3)/11+9./22,                    0 },
									{   3*sqrt(3)/22-9./11,  -9*sqrt(3)/11+9./22,                    0 },
};

#endif

