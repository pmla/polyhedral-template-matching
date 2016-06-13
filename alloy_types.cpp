#include "index_ptm.h"

int32_t find_fcc_alloy_type(int8_t* mapping, int32_t* numbers)
{
	int len = 13;
	int8_t temp[13];

	int num_cu = 1, cu = numbers[0];
	for (int i=1;i<len;i++)
		if (numbers[i] == cu)
			num_cu++;
	if (num_cu == len)
		return PTM_ALLOY_PURE;

	int num_au = 0, au = -1;
	for (int i=1;i<len;i++)
	{
		if (numbers[i] != cu)
		{
			if (au == -1)
			{
				au = numbers[i];
				num_au = 1;
			}
			else if (numbers[i] == au)
			{
				num_au++;
			}
			else
			{
				return PTM_ALLOY_NONE;
			}
		}
	}

	if (num_au == len - 1)
		return PTM_ALLOY_L12_AU;

	for (int i=0;i<len-1;i++)
		temp[i] = numbers[mapping[i + 1]];

	if (num_au == 4)
	{
		for (int j = 0;j<3;j++)
		{
			int n = 0;
			for (int i=j*4;i<(j+1)*4;i++)
				if (temp[i] == au)
					n++;
			if (n == 4)
				return PTM_ALLOY_L12_CU;
		}
	}
	else if (num_au == 8)
	{
		for (int j = 0;j<3;j++)
		{
			int n = 0;
			for (int i=j*4;i<(j+1)*4;i++)
				if (temp[i] == cu)
					n++;
			if (n == 4)
				return PTM_ALLOY_L10;
		}
	}

	return PTM_ALLOY_NONE;
}

int32_t find_bcc_alloy_type(int8_t* mapping, int32_t* numbers)
{
	int len = 15;
	int8_t temp[15];

	int num_cu = 1, cu = numbers[0];
	for (int i=1;i<len;i++)
		if (numbers[i] == cu)
			num_cu++;
	if (num_cu == len)
		return PTM_ALLOY_PURE;

	for (int i=0;i<len;i++)
		temp[i] = numbers[mapping[i + 1]];

	int num_au_inner = 0, num_cu_outer = 0, au = -1;
	for (int i=0;i<8;i++)
	{
		if (temp[i] != cu)
		{
			if (au == -1)
			{
				au = temp[i];
				num_au_inner = 1;
			}
			else if (temp[i] == au)
			{
				num_au_inner++;
			}
		}
	}


	for (int i=8;i<14;i++)
		if (temp[i] == cu)
			num_cu_outer++;

	if (num_au_inner == 8 && num_cu_outer == 6)
		return PTM_ALLOY_B2;

	return PTM_ALLOY_NONE;
}

