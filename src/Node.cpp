#include "Node.h"
Node::Node(int* alphabet_sizes)
{
    // compute auxiliary table
    int tmp_byte = 0;
    int i, j;
    for(i = 0; i < DIM; i++)
    {
        DMBR_start_byte_lut[i] = tmp_byte;
        for(j = 0; j < alphabet_sizes[i]; j++)
        {
            DMBR_byte_lut[i][j] = tmp_byte + j / BITS_PER_BYTE;
            DMBR_bit_lut[i][j] = j % BITS_PER_BYTE;
        }
        tmp_byte += (alphabet_sizes[i] - 1) / BITS_PER_BYTE + 1;
        DMBR_end_byte_lut[i] = tmp_byte - 1;
    }
    for(i = 0; i < MAX_ALPHABET_SIZE; i++)
    {
        bitmap_byte_lut[i] = i / BITS_PER_BYTE; 
        bitmap_bit_lut[i] = i % BITS_PER_BYTE; 
    }
}

double Node::cal_normed_area(unsigned char* DMBR, int* alphabet_sizes)
{
    int i, j;
    int tmp_sum;
    double area = 1.0;
    for(i = 0; i < DIM; i++)
    {
        tmp_sum = 0;
        for(j = DMBR_start_byte_lut[i]; j <= DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[DMBR[j]];
        if(tmp_sum == 0)
            return 0.0;
        area *= static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    return area;
}


double Node::cal_normed_overlap(unsigned char* DMBR1, unsigned char* DMBR2, int* alphabet_sizes)
{
    int i, j;
    int tmp_sum;
    double overlap = 1.0;
    for(i = 0; i < DIM; i++)
    {
        tmp_sum = 0;
        for(j = DMBR_start_byte_lut[i]; j <= DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[DMBR1[j] & DMBR2[j]];
        if(tmp_sum == 0) // tmp_sum is zero
            return 0.0;
        else
            overlap *= static_cast<double>(tmp_sum) / alphabet_sizes[i];
    }
    return overlap;
}

bool Node::set_equal(unsigned char* set1, unsigned char* set2, int byte_size)
{
    int i;
    for(i = 0; i < byte_size; i++)
        if(set1[i] != set2[i])return false;
    return true;
}

void Node::set_intersect(unsigned char* set1, unsigned char* set2, int byte_size, unsigned char* intersect_set, bool& is_intersect)
{
    int i;
    is_intersect = false;
    for(i = 0; i < byte_size; i++)
    {
        intersect_set[i] = set1[i] & set2[i];
        if(intersect_set[i] > 0)is_intersect = true;
    }
}

void Node::bitmap_to_letters(unsigned char* bitmap_set, int byte_size, unsigned char* letters, int &letter_count)
{
    int i, j;
    int tmp_letter_inc;
    int tmp_letter;
    letter_count = 0;
    tmp_letter_inc = 0;
    for(i = 0; i < byte_size; i++)
    {
        for(j = 0; j < BITS_PER_BYTE; j++)
	    if(bitmap_set[i] & MASKS[j])
            {
                tmp_letter = tmp_letter_inc + j;
                letters[letter_count] = tmp_letter;
                letter_count++;            
            }
        tmp_letter_inc += BITS_PER_BYTE;
    }
}

int Node::set_size(unsigned char* bitmap_set, int byte_size)
{
    int i;
    int letter_count = 0;

    for(i = 0; i < byte_size; i++)
        letter_count += bit_1_count_lut[bitmap_set[i]];

    return letter_count;
}

void Node::log_DMBR(const unsigned char* const DMBR)
{
    int tmp_sum;
    unsigned char tmp_DMBR[DMBR_SIZE];

    for(int i=0;i<DMBR_SIZE;i++)
        tmp_DMBR[i] = DMBR[i];

    for(int i = 0; i < DIM; i++)
    {
        tmp_sum = 0;
        for(int j = DMBR_start_byte_lut[i]; j <= DMBR_end_byte_lut[i]; j++)
            tmp_sum += bit_1_count_lut[DMBR[j]];
        logO.log2File(tmp_sum);/*logO.log2File(" ");*/
    }
    logO.log2File("\n");

}


void Node::log_OneDirEntry_OnOneDscDim(  unsigned char*  DMBR, int dimNum,int* alphabet_sizes )
{


	unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
	int size_of_alphabet2 = 0;


	assert(dimNum<DIM);

	Node::bitmap_to_letters(&(DMBR[dimNum*BYTES_PER_DIM_IN_DMBR]), BYTES_PER_DIM_IN_DMBR, alphabet2, size_of_alphabet2);


	for(int j=0;j<size_of_alphabet2;j++)
	{
		logO.log2File((int)alphabet2[j]);logO.log2File(" ");
	}

}





