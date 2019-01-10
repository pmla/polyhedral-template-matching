#ifndef PTM_ARB_TREEFILTER_H
#define PTM_ARB_TREEFILTER_H


namespace ptm {

typedef struct
{
	const int num_rows;
	const int row_starts[3];
	const int row_lengths[3];
	const uint8_t (*data)[3];
} treefilter_t;

const uint8_t data_fluorite_ca[16][3] = 
{
	{0},

	{0, 1}, {0, 9},

	{0, 1, 2}, {0, 1, 4}, {0, 1, 8}, {0, 1, 9}, {0, 1, 12}, {0, 1, 15},
	 {0, 9, 1}, {0, 9, 3}, {0, 9, 7}, {0, 9, 10}, {0, 9, 12}, {0, 9, 13}, {0, 9, 15}
};

const treefilter_t filter_fluorite_ca = {3, {0, 1, 3}, {1, 2, 13}, data_fluorite_ca};


const uint8_t data_fluorite_f[27][3] = 
{
	{0},

	{0, 1}, {0, 5}, {0, 11},

	{0, 1, 2}, {0, 1, 5}, {0, 1, 8}, {0, 1, 11}, {0, 1, 14}, {0, 1, 20},
	 {0, 5, 1}, {0, 5, 2}, {0, 5, 6}, {0, 5, 8}, {0, 5, 11}, {0, 5, 13}, {0, 5, 16},
	 {0, 11, 1}, {0, 11, 2}, {0, 11, 3}, {0, 11, 5}, {0, 11, 7}, {0, 11, 8}, {0, 11, 12}, {0, 11, 14}, {0, 11, 18}, {0, 11, 20}
};

const treefilter_t filter_fluorite_f = {3, {0, 1, 4}, {1, 3, 23}, data_fluorite_f};

}

#endif

