#ifndef NDTREE_HELPER_H
#define NDTREE_HELPER_H

// for fasta parser and zlib
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "Dir_entry.h"
#include "ND_tree.h"
#include "Leaf_node.h"
#include <vector>
#include "logClass.h"
#include <limits.h>
#include "config.h"
#include <time.h>
#include <iostream>
using namespace std;

class Item
{
public:
	int id;
	int count;
};


class ndtreeHelper {
public:
ofstream OutStream;

// GLOBAL VARIABLES - Added by Alok to make interface of this
// tree consistent with the rest of the trees.
//#define MAX_COPIES_OF_DATA_FILE 50

//double dir_min_util = 0.0000000000003;
const double dir_min_util = ((nodeSplitType==ORIGINAL)||(enforce_minUtil_for_exhaustSplit==1))?0.3:0.0000000000003;
const double leaf_min_util = 0.3;

int tmp_DMBR_byte_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding byte in DMBR
//int tmp_DMBR_byte_lut[MAX_K][MAX_ALPHABET_SIZE]; //changed to support multiple K 
int tmp_DMBR_bit_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding bit in DMBR
//int tmp_DMBR_bit_lut[MAX_K][MAX_ALPHABET_SIZE]; //changed to support mutiple K

//debug_boxQ_leaf_accessed=0;
//int debug_boxQ_leaf_hit=0;
//int debug_boxQ_leaf_hit_peak=0;
vector<int> debug_boxQ_leaf_hit_for_all;

int debug_height;

int duplicateDataPoints; //used in batchBuild(...), batchGrow(...)

ndtreeHelper();


int letter2num(char letter);
char num2letter(int num);



void LocalDMBRInfoCalculation();
//Dir_entry* makeRandomBoxQueryData(ifstream & query_file, int K, char queryK[]);
Dir_entry* makeRandomBoxQueryData(char *seq, int K, char queryK[]);
Dir_entry makeBoxQueryData(ifstream & query_file);
void batchGrow_with_duplicate(long skipSize, long sizeGrowTo);
//this one only reads data file up to the give number of data points 
void batchBuild_with_duplicate(long  size);
void clear_record();
//Insert kmers with a link to the record(readid,annotation)
void batchBuild_with_duplicate_record(bool newtree, long  size);
void batchRangeQuery();
//output the cancer types in the query result to file
void output_records(Leaf_entry* query_results, int query_results_size);
//output the cancer types in the query result to file
void output_records(Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE], int query_results_size[], char queryK[], int maxShift);
string output_records_fasta(Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE], int query_results_size, char query_id[],  char queryK[]);
//output the cancer types in the query result to file
void output_records_array(Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE], int query_results_size[], char queryK[], int maxShift);

void clear_result();
void batchRandomBoxQuery(int K);
void onlineBoxQuery(int K, char *queryStr);
void batchBoxQuery();
vector< vector<int> >  makeBoxQueryData_for_linearScan(ifstream & query_file);
//result of this LinearScanBoxQuery should be 
void LinearScanBoxQuery(int dataNUM);

void display_help();


};





ndtreeHelper::ndtreeHelper() {
	debug_boxQ_leaf_accessed=0;
 	debug_boxQ_leaf_hit=0;
 	debug_boxQ_leaf_hit_peak=0;
}

int ndtreeHelper::letter2num(char letter){
    int num = 0;
    switch(letter){
	case 'A':
	    num = 0;
	    break;
	case 'T':
	    num = 1;
	    break;
	case 'G':
	    num = 2;
	    break;
	case 'C':
	    num = 3;
	    break;
        case 'N':
	    num = 4;
	    break;
	case UNI_ELEM:
	    num = UNI_ELEM_NUM; // MUST be continous, 5
	    break;
	default:
	    num = 6;
    }

    return num;
}

char ndtreeHelper::num2letter(int num){
    char letter;
    switch(num){
	case 0:
	    letter = 'A';
	    break;
 	case 1:
	    letter = 'T';
	    break;
	case 2:
	    letter = 'G';
	    break;
	case 3:
	    letter = 'C';
	    break;
	case 4:
	    letter = 'N';
	    break;
	case UNI_ELEM_NUM:
	    letter = UNI_ELEM;
	    break;
	default:
	    letter = '?';
    }
    return letter;

}

void ndtreeHelper::LocalDMBRInfoCalculation()
{
    int tmp_byte = 0;
    int i, j;
    for(i = 0; i < DIM; i++)
    {
        //tmp_DMBR_start_byte_lut[i] = tmp_byte;

        for(j = 0; j < A[i]; j++)
        {
            tmp_DMBR_byte_lut[i][j] = tmp_byte + j / BITS_PER_BYTE;
            tmp_DMBR_bit_lut[i][j] = j % BITS_PER_BYTE; ////letters must starts from 0
        }


        /** how many bytes occupied by this DBMR[i][j] **/
        tmp_byte += (A[i] - 1) / BITS_PER_BYTE + 1;

        //tmp_DMBR_end_byte_lut[i] = tmp_byte - 1;
    }



}

Dir_entry* ndtreeHelper::makeRandomBoxQueryData(char *seqStr, int K, char queryK[])
{
    int queryK_length = 0;
    Dir_entry* dirEntry = new Dir_entry[MAX_K];
    for(int i=0; i < MAX_K; i++)
    	for(int j = 0; j < DMBR_SIZE; j++)
            dirEntry[i].DMBR[j]=0;
    
    string line[MAX_K]; // length is equal to K

    int byte_no,bit_no;
/*
    for(int i=0;i<K;i++)
    {
        getline(query_file, line[i]);
	istringstream instr(line[i]);
	char v;
	instr >> v;//consume the boxsize in the head of line
	instr >> v;//read the base 
	queryK[i] = v;
    }

    for (int i=0; i <= K - DIM; i++) {
	for(int j=0; j < DIM; j++){
            istringstream instr(line[j+i]);
            int v;
            int box_size;
            instr>>box_size;
		
	    for(int t=0;t<box_size;t++)
	    {
		instr>>v;
		byte_no = tmp_DMBR_byte_lut[j][v];
		bit_no = tmp_DMBR_bit_lut[j][v];
		dirEntry[i].DMBR[byte_no] |= MASKS[bit_no];
	    }
      }
*/

int cur = 0;
int j = 0;
int v = 0;
    
    while (j < DIM) {
	if(seqStr[cur] == '('){
	    while(seqStr[++cur] != ')'){
		v = letter2num(seqStr[cur]);
		byte_no = tmp_DMBR_byte_lut[j][v];
                bit_no = tmp_DMBR_bit_lut[j][v];
                dirEntry[0].DMBR[byte_no] |= MASKS[bit_no];
 	    }
	} else {
	    v = letter2num(seqStr[cur]); 
	    byte_no = tmp_DMBR_byte_lut[j][v];
            bit_no = tmp_DMBR_bit_lut[j][v];
            dirEntry[0].DMBR[byte_no] |= MASKS[bit_no];
	}
	cur++;
	j++;

    }




  // }



//Why to consume the rest lines. What are the float values? 
/* 	
    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);
    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);
*/
    return dirEntry;
}

Dir_entry ndtreeHelper::makeBoxQueryData(ifstream & query_file)
{
    Dir_entry dirEntry;

    for(int j = 0; j < DMBR_SIZE; j++)
        dirEntry.DMBR[j]=0;
    string line;

    for(int j=0;j<DIM;j++)
    {
        int byte_no,bit_no;
        getline(query_file, line);
        istringstream instr(line);
        int v;
        for(int t=0;t<boxSize;t++)
        {
            instr>>v;
            byte_no = tmp_DMBR_byte_lut[j][v];
            bit_no = tmp_DMBR_bit_lut[j][v];/** todo: here DMBR_bit_lut[j] should support enough letters **/
            dirEntry.DMBR[byte_no] |= MASKS[bit_no];
        }
    }
    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);
    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);
    return dirEntry;
}

void ndtreeHelper::batchGrow_with_duplicate(long skipSize, long sizeGrowTo)
{
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    ndt.read_existing_tree(globalIndexFilename);
    ifstream data_file;
    data_file.open(globalDataFilename.c_str());
    if (data_file.fail())
    {
        cout<<"cant open file "<<globalDataFilename<<endl;
        exit(1);

    }
    int number_of_io = 0;
    long distinctDataPoints=0;
    string line;
    getline(data_file, line);
    for(long i = 0; i < sizeGrowTo; i++)
    {
#ifdef LOG_VECTOR_INDEX
logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
        vector_index = i;
#endif
        if(i<skipSize)
        {
            getline(data_file, line);
            continue;
        }
        else
        {
            stringstream instr(line);
           for(int j = 0; j < DIM; j++)
           {
               int n;
               instr >> n;
               new_data.key[j] = n - '0';
           }
           new_data.record_count = 1; 
           result = ndt.insert_use_link(new_data, number_of_io);
           if(result == duplicate_error)
            {
                duplicateDataPoints++;
            }
            else
            {
                distinctDataPoints++;
            }
            getline(data_file, line);
        }
    }
    data_file.close();
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<sizeGrowTo <<endl;
    ndt.print_information( );
}

//this one only reads data file up to the give number of data points 
void ndtreeHelper::batchBuild_with_duplicate(long  size)
{
    duplicateDataPoints=0;
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    ndt.create_empty_tree(A, dir_min_util, leaf_min_util, globalIndexFilename);
    ndt.read_existing_tree(globalIndexFilename);
    long num_of_points=size;
    ifstream data_file;
    data_file.open(globalDataFilename.c_str());
    if (data_file.fail())
    {
        cout<<"cant open file "<<globalDataFilename<<endl;
        exit(1);

    }

    int number_of_io = 0;
    string line;
    getline(data_file, line);
    int distinctDataPoints = 0;
    char n;
    for(long  i = 0; i < num_of_points && !data_file.eof(); i++)
    {
#ifdef LOG_VECTOR_INDEX
        logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif
#ifdef LOG_SPLITTING_VECTOR_INDEX
        vector_index = i;
#endif
        istringstream instr(line);
        for(int j = 0; j < DIM; j++)
        {
            instr >> n;
            new_data.key[j] = n - '0';
        }
        new_data.record_count = 1; 
	//new_data.record = 0;
        result = ndt.insert_use_link(new_data, number_of_io);
        if(result == duplicate_error)
        {
            duplicateDataPoints++;
        }
        else
        {
            distinctDataPoints++;
        }
        getline(data_file, line);
    }
    data_file.close();
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<size <<endl;
    ndt.print_information( );
}


void ndtreeHelper::clear_record()
{
fstream typeid_file;
const char* typeid_filename = (globalRecordFilename+".typeid").c_str();
typeid_file.open(typeid_filename, ios::binary | ios::out | ios::trunc);
if(typeid_file.fail())
{
    cout<<"can't open record.typeid file "<< typeid_filename <<endl;
    exit(1);
}
typeid_file.clear();
typeid_file.close();
}




//Insert kmers with a link to the record(readid,annotation)
void ndtreeHelper::batchBuild_with_duplicate_record(bool newtree, long  size)
{
    clear_record();

    duplicateDataPoints=0;
    int distinctDataPoints = 0;
    ND_tree ndt;
    Leaf_entry new_data/*, query_data, db_data*/;
    Error_code result;
    if(newtree) {
	cout << "create new index tree"<< endl;
        ndt.create_empty_tree(A, dir_min_util, leaf_min_util, globalIndexFilename);
    }
    ndt.read_existing_tree(globalIndexFilename);
    long num_of_points=size;

    gzFile kmer_file;
    kseq_t *seq;
    int tmp;

    kmer_file = gzopen(globalDataFilename.c_str(), "r");
    seq = kseq_init(kmer_file);

    int count_records = 0;
    int total_records = 0;
    while((tmp = kseq_read(seq)) >= 0)
	total_records++; 

    gzclose(kmer_file);
    // reset the file pointer to the beginning
    kmer_file = gzopen(globalDataFilename.c_str(), "r");
    seq = kseq_init(kmer_file);

    int number_of_io = 0;

    time_t start, pause;
    time(&start);

    while((tmp = kseq_read(seq)) >= 0)
    {
	// show rate of progress, very useful for big data
	count_records++;
	if(count_records % 50 == 0) {
		double percent = count_records * 100.0 / total_records;
		time(&pause);	
		double diff = difftime(pause, start);
		int remaining =(int)(diff * (total_records - count_records)/ count_records);
		int min = remaining / 60;
		int sec = remaining % 60;
		double avgLinksPerInsert = total_links_num_global * 1.0 / valid_insert_num_global;
		printf("\r%.0lf%%\t%d min %d sec remaining\tAvg LPI: %.2lf", percent, min, sec, avgLinksPerInsert);

	}
	for (int j = 0; j < seq->seq.l; j++) 
	    new_data.key[j] = letter2num(seq->seq.s[j]);

	typeid_global = atoi(seq->name.s);

	new_data.record_count = 1; 
        result = ndt.insert_use_link(new_data, number_of_io);
        if(result == duplicate_error)
        {
            duplicateDataPoints++;
        }
        else
        {
            distinctDataPoints++;
        }
    }
    cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;
    cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints<<endl;
    cout<<"Total data points read:"<<size <<endl;
    ndt.print_information( );

    // close kmer file and destroyt seq
    kseq_destroy(seq);
    gzclose(kmer_file);
}


void ndtreeHelper::batchRangeQuery()
{

    Error_code result;

    string input=globalIndexFilename;
    ND_tree ndt;
    Leaf_entry query_data;

    ndt.read_existing_tree(input);


                Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
                int query_results_size;

    string query_fn = globalRQFilename;

    string str_num_of_points;
    int num_of_points;
    ////cout << endl << "How many query records:" << endl;
    ////getline(cin, str_num_of_points);
    ////num_of_points = string_to_int(str_num_of_points);

    num_of_points=TOTAL_RANGE_QUERY_NUM;

    for(int range=0;range<RANGE_STOP_BEFORE;range+=RANGE_SIZE_STEP)
    {
        //start = clock();

        ifstream query_file;
        query_file.open(query_fn.c_str(),ios::in);
        if (query_file.fail())
        {
            cout<<"cant open file "<<query_fn<<endl;
            exit(1);

        }                

        int number_of_io = 0;
        int total_number_of_io = 0;
        int total_results_size = 0;

        string line;
        int n;
        for(int i = 0; i < num_of_points; i++)
        {
            query_results_size = 0;

            //cout <<endl<<"range query for query point "<< i << endl;
            getline(query_file, line);

            if(query_file.bad()||query_file.fail() )

            {
                cout<<"reading query file error\n";
                exit(0);
            }

            istringstream instr(line);



            for(int j = 0; j < DIM; j++)
            {
                instr >> n;
                query_data.key[j] = n;
            }


            result = ndt.range_query_by_Hamming_dist(query_data, range, query_results, query_results_size, number_of_io);

            //cout << "  " << query_results_size << " " << number_of_io << endl;
            total_number_of_io += number_of_io;
            total_results_size += query_results_size;
        }
        query_file.close();

        //finish = clock();
        cout<<total_number_of_io<<" : "<<num_of_points<<endl;
        //elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
        //cout << "Time (in seconds): " << elapsed_time / num_of_points << endl;
        cout << "For range: "<<range<<" range query I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
        cout << "Result size: " << static_cast<double>(total_results_size) / num_of_points << endl;




    }//end of for(int range=1;range<RANGE_STOP_BEFORE;++range)


    if(RUNNING_ENVIRONMENT == WINDOWS)
        system("pause");
}

//output the cancer types in the query result to file
void ndtreeHelper::output_records(Leaf_entry* query_results, int query_results_size)
{
    fstream typeid_file, query_result_file;
    const char* typeid_filename = (globalRecordFilename+".typeid").c_str();
    const char* query_result_filename = (globalBQFilename+ ".result").c_str();
    typeid_file.open(typeid_filename, ios_base::binary | ios_base::in);
    query_result_file.open(query_result_filename, ios_base::app);

    if(typeid_file.fail())
    {
	cout<<"can't open file "<<typeid_filename<<endl;
    }
    if(query_result_file.fail())
    {
	cout<<"can't open file "<<query_result_filename<<endl;
    }

    int type_array[TYPE_ARRAY_SIZE]={0};
    int output[256]={0};
    output[0]=0;
    int type_num;
    int record_id;

for(int k=0; k<query_results_size; k++)
{
    record_id = query_results[k].record;

    typeid_file.seekg(record_id * sizeof(type_array), ios::beg);
    typeid_file.read((char*)type_array, sizeof(type_array));

    type_num = type_array[0];
    int i;

    for(i=1; i<=type_num; i++)
    {
	int p;
	for(p=1; p<=output[0]; p++)
	{
	    if(type_array[i]==output[p])
		break;
	}
	if(p>output[0])
	{
	    output[0]++;
	    output[output[0]]=type_array[i];
	}
    }
}
 
    query_result_file << output[0]<<" ";
    for(int i=1; i<=output[0]; i++)
	query_result_file << output[i]<<" ";
    query_result_file << endl;

    typeid_file.close();
    query_result_file.close();

}


//output the cancer types in the query result to file
void ndtreeHelper::output_records(Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE], int query_results_size[], char queryK[], int maxShift)
{
  string typeid_filename = (globalRecordFilename+".typeid").c_str();
  string query_result_filename = (globalBQFilename +".result").c_str();


const int item_size = 256;
int item_count = 0;
Item item[item_size];
for(int i=0; i<item_size; i++)
{
	item[i].id = -1;
	item[i].count = 0;
}

    fstream typeid_file, query_result_file;
  
    typeid_file.open(typeid_filename.c_str(), ios_base::binary | ios_base::in);
    query_result_file.open(query_result_filename.c_str(), fstream::in | fstream::out| fstream::app);

    if(typeid_file.fail()) {
	cout << "can't open file typeid " << typeid_filename << endl;
    }
    if(query_result_file.fail())
    {
	cout<<"can't open file query result"<<query_result_filename<<endl;
    }

    int type_array[TYPE_ARRAY_SIZE]={0};
    int output[256]={0};
    output[0]=0;
    int type_num;
    int record_id;

for(int m = 0; m <= maxShift; m++) {
for(int k=0; k<query_results_size[m]; k++){
    record_id = query_results[m][k].record;
    typeid_file.seekg(record_id * sizeof(type_array), ios::beg);
    typeid_file.read((char*)type_array, sizeof(type_array));
    type_num = type_array[0];
    int i;
if(type_num > TYPE_ARRAY_SIZE - 1) {
	cout << "type_num > " << TYPE_ARRAY_SIZE - 1 << " :"<<type_num<<endl;
//	cout << query_K <<endl;
	type_num = TYPE_ARRAY_SIZE - 1;
}
    for(i=1; i<=type_num; i++){
	if(m==0) {
	    if(item_count > item_size - 1) {
		cout << "item overflow!" <<endl;
		continue;
	    }
	    item[item_count].id = type_array[i];
    	    item[item_count++].count = 0;
	}else {
	    for(int d = 0; d < item_count; d++) {
		if(item[d].id == type_array[i]){
		   item[d].count++;
		   break;
		}
	    }
	}

    }
}
}

    int valid_read_count = 0;

    for(int i = 0; i < item_count; i++) {
	if(item[i].count == maxShift){
	    valid_read_count++;
	}
    }
    if(valid_read_count > 0){
	//query_result_file << valid_read_count<< " ";
        for(int i = 0; i < item_count; i++) {
	    if(item[i].count == maxShift){
	        query_result_file << item[i].id << " ";
	    }
        }
        query_result_file << endl;
 
        int ct = 0;
        while(true){
 	    if(queryK[ct] == '\0')
		break;

        query_result_file << queryK[ct++] << " ";
        }
        query_result_file << endl;

   }
/*
    query_result_file << output[0]<<" ";
    for(int i=1; i<=output[0]; i++)
	query_result_file << output[i]<<" ";
 
    query_result_file << endl;
*/
    typeid_file.close();
    query_result_file.close();

}


//output the cancer types in the query result to file
string ndtreeHelper::output_records_fasta(Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE], int query_results_size, char query_id[], char queryK[])
{
  string typeid_filename = (globalRecordFilename+".typeid").c_str();
  string query_result_filename = (globalBQFilename +".result").c_str();

    fstream typeid_file, query_result_file;
  
    typeid_file.open(typeid_filename.c_str(), ios_base::binary | ios_base::in);
    query_result_file.open(query_result_filename.c_str(), fstream::in | fstream::out| fstream::app);

    if(typeid_file.fail()) {
	cout << "can't open file typeid " << typeid_filename << endl;
    }
    if(query_result_file.fail())
    {
	cout<<"can't open file query result"<<query_result_filename<<endl;
    }

    int type_array[TYPE_ARRAY_SIZE]={0};
    int output[256]={0};
    output[0]=0;
//    int type_num;
    int record_id;

string retStr = "";

for(int k=0; k<query_results_size; k++){
    string kmer="";
    for(int p = 0; p < DIM; p++){
	kmer += num2letter(query_results[k].key[p]);
    }
    string tuple = ">" + string(query_id) + " "+string(queryK)+" "+kmer+"\n";

    record_id = query_results[k].record;
    typeid_file.seekg(record_id * sizeof(type_array), ios::beg);
    typeid_file.read((char*)type_array, sizeof(type_array));
    //type_num = type_array[0];
    //int i;
    /*if(type_num > TYPE_ARRAY_SIZE - 1) {
	cout << "type_num > " << TYPE_ARRAY_SIZE - 1 << " :"<<type_num<<endl;
	type_num = TYPE_ARRAY_SIZE - 1;
    }
    for(i=1; i<type_num; i++){
	tuple += to_string(type_array[i]) + ",";
    }
    
    tuple += to_string(type_array[type_num]);
    */
//    cout<<tuple<<endl;

    //replace
    int cur = 0;
    int i = 0;
    string idStr = "";
    while(type_array[i] != END_ARRAY){
	if(i == TYPE_ARRAY_SIZE -1){
	   cur = type_array[i];
	   typeid_file.seekg(cur * sizeof(type_array), ios::beg);
	   typeid_file.read((char*)type_array, sizeof(type_array));
	   i = 0;
	   continue;
	}
	idStr += to_string(type_array[i]) + ",";
	i++;
    }
    tuple += idStr.substr(0, idStr.size()-1);
    query_result_file << tuple << endl;
    retStr += tuple;
}

    typeid_file.close();
    query_result_file.close();
    return retStr;
}

/*
//output the cancer types in the query result to file
void ndtreeHelper::output_records_array(Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE], int query_results_size[], char queryK[], int maxShift)
{

const int item_size = 256;
int item_count = 0;
Item item[item_size];
for(int i=0; i<item_size; i++)
{
	item[i].id = -1;
	item[i].count = 0;
}

    fstream typeid_file, query_result_file;
    const char* query_result_filename = (globalBQFilename+ ".result").c_str();
    query_result_file.open(query_result_filename, fstream::in | fstream::out| fstream::app);

    if(query_result_file.fail())
    {
	cout<<"can't open file "<<query_result_filename<<endl;
    }

    int type_array[TYPE_ARRAY_SIZE]={0};
    int output[256]={0};
    output[0]=0;
    int type_num;
    int record_id;

for(int m = 0; m <= maxShift; m++) {
for(int k=0; k<query_results_size[m]; k++){
    record_id = query_results[m][k].record;
    type_num = record_type[record_id][0];
    int i;
if(type_num > TYPE_ARRAY_SIZE -1) {
	cout << "type_num > "<< TYPE_ARRAY_SIZE -1 <<" :"<<type_num<<endl;
//	cout << query_K <<endl;
	type_num = TYPE_ARRAY_SIZE - 1;
}
    for(i=1; i<=type_num; i++){
	if(m==0) {
	    if(item_count > item_size -1) {
		continue;
	    }
	    item[item_count].id = record_type[record_id][i];
    	    item[item_count++].count = 0;
//	    cout << "m==0 : " << record_type[record_id][i] << endl;
	}else {
	    for(int d = 0; d < item_count; d++) {
		if(item[d].id == record_type[record_id][i]){
		   item[d].count++;
		   break;
		}
	    }
//	    cout << "m==1 : " << record_type[record_id][i] << endl;
	}

    }
}
}

    int valid_read_count = 0;

    for(int i = 0; i < item_count; i++) {
	if(item[i].count == maxShift){
	    valid_read_count++;
	}
    }
    if(valid_read_count > 0){
	//query_result_file << valid_read_count<< " ";
        for(int i = 0; i < item_count; i++) {
	    if(item[i].count == maxShift){
	        query_result_file << item[i].id << " ";
	    }
        }
        query_result_file << endl;
 
        int ct = 0;
        while(true){
 	    if(queryK[ct] == '\0')
		break;

        query_result_file << queryK[ct++] << " ";
        }
        query_result_file << endl;

   }
   typeid_file.close();
    query_result_file.close();

}
*/

void ndtreeHelper::clear_result()
{
fstream  query_result_file;
string tmp = globalBQFilename;
const char* query_result_filename = (tmp.append(".result")).c_str();
query_result_file.open(query_result_filename, ios::out);
if(query_result_file.fail())
{
    cout<<"can't open result file "<<query_result_filename<<endl;
    exit(1);
}
query_result_file.clear();
query_result_file.close();
}



void ndtreeHelper::onlineBoxQuery(int K, char *queryStr)
{
    string input=globalIndexFilename;
    ND_tree ndt;
    ndt.read_existing_tree(input);
    LocalDMBRInfoCalculation();
    
    Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE];
    int query_results_size[MAX_K];
    Dir_entry *boxQueryData;
    char queryK[255]; 
boxQueryData = makeRandomBoxQueryData(queryStr, K, queryK);
int maxShift = K - DIM;
int number_of_io ;
for(int i=0; i <= maxShift; i++) {
	ndt.box_query(boxQueryData[i], query_results[i], query_results_size[i], number_of_io);
}

int total_record_num = 0;
//calculate total record in the result
for(int i=0; i <= maxShift; i++) {
    for(int k = 0; k < query_results_size[i]; k++)
	total_record_num += query_results[i][k].record_count;

}

if(query_results_size[0] == 0){
	cout << "Not found!" << endl;
} else {
// Return output string !!!!
	string output = "Found: \n";
        char tmp[1] = {'i'};
 	output += output_records_fasta(query_results[0], query_results_size[0],tmp, queryStr);
	cout << output << endl;
}

}








void ndtreeHelper::batchRandomBoxQuery(int K)
{
clear_result();
    string input=globalIndexFilename;
    ND_tree ndt;
    ndt.read_existing_tree(input);
    LocalDMBRInfoCalculation();
    
    Leaf_entry query_results[MAX_K][QUERY_RESULTS_BUFFER_SIZE];
    int query_results_size[MAX_K];
    Dir_entry *boxquerydata;

    int num_of_points=TOTAL_BOX_QUERY_NUM;

    gzFile query_file;
    kseq_t *seq;
    int tmp;
    query_file = gzopen(globalBQFilename.c_str(), "r");
    seq = kseq_init(query_file);


    debug_boxQ_leaf_hit_for_all.clear();
    debug_boxQ_leaf_accessed=0;

    int number_of_io ;
    int total_number_of_io = 0;
    int total_results_size=0;
    int total_record_num = 0;
   
    while((tmp = kseq_read(seq)) >= 0)
    {
        debug_boxQ_leaf_hit_peak=0;
	char queryK[MAX_K] = {'\0'};

        boxquerydata = makeRandomBoxQueryData(seq->seq.s, K, queryK);
	int maxshift = K - DIM;
	for(int i=0; i <= maxshift; i++) {
	        ndt.box_query(boxquerydata[i], query_results[i], query_results_size[i], number_of_io);
	}


	//calculate total record in the result
	for(int i=0; i <= maxshift; i++) {
            for(int k = 0; k < query_results_size[i]; k++)
		total_record_num += query_results[i][k].record_count;

	}

	// output_records(query_results, query_results_size);


        total_number_of_io += number_of_io;
        total_results_size += query_results_size[0];
        debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);

 	//output_records_array(query_results, query_results_size, queryk, maxshift);
 	output_records_fasta(query_results[0], query_results_size[0],seq->name.s, seq->seq.s);

    }

//    query_file.close();

kseq_destroy(seq);
gzclose(query_file);


    cout<<"boxsize= random, avg boxquery i/o: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
  cout<<" avg matched data point found="<< static_cast<double>(total_results_size)/num_of_points<< endl; 
    cout << " avg leaf node accessed: " << static_cast<double>(debug_boxQ_leaf_accessed) / num_of_points << endl;
  cout<<"total boxquery i/o="<<static_cast<double>(total_number_of_io)<<endl;
     cout<<"total matched data points="<< static_cast<double>(total_results_size)<< endl; 
     cout<<"total matched records = "<< static_cast<double>(total_record_num)<<endl;
   
    //assert(debug_boxq_leaf_hit_for_all.size()==num_of_points);
    int debug_tatol=0;
    for(unsigned int i=0;i<debug_boxQ_leaf_hit_for_all.size();i++)
    {
        
        debug_tatol+=debug_boxQ_leaf_hit_for_all.at(i);
    
    }

    cout << " avg leaf node hit peak: " << static_cast<double>(debug_tatol) / num_of_points << endl;
    
}


vector< vector<int> >  ndtreeHelper::makeBoxQueryData_for_linearScan(ifstream & query_file)
{

    //int maxdimandcontdim = 16;//01/09/2007


    string line;

    vector< vector<int> >  boxquerydata;
    boxquerydata.resize(DIM);

    for(int j=0;j<DIM;j++)
    {

        getline(query_file, line);
        istringstream instr(line);
        boxquerydata[j].clear();

        int v;
        for(int t=0;t<(boxSize*A[j])/10;t++)
        {
            instr>>v;
            boxquerydata[j].push_back(v);

        }
    }


    for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
        getline(query_file, line);


    //consumes the rest lines holding float values
    for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
        getline(query_file, line);

    assert(boxquerydata.size()==DIM);
    assert(boxquerydata[0].size()==boxSize);

    return boxquerydata;
}








//result of this linearscanboxquery should be 
void ndtreeHelper::LinearScanBoxQuery(int datanum)
{
//    leaf_entry query_results [query_results_buffer_size];


//    int query_results_size;


    int db_size=datanum;
    string query_fn="c:\\temp\\boxqueryall.txt";


    ifstream query_file;
    query_file.open(query_fn.c_str());
    if (query_file.fail())
    {
        cout<<"cant open file "<<query_fn<<endl;
        exit(1);

    }




    for(boxSize=1;boxSize<=9;boxSize+=1)
    {


        query_file.seekg(0, ios::beg);


        int num_of_points=200;

        int totalnumofmatches=0;

        for(int i = 0; i < num_of_points; i++) // for every query point
        {

            vector< vector<int> > matchedpoints;
            matchedpoints.clear();

            vector< vector<int> >  boxquerydata;
            boxquerydata = makeBoxQueryData_for_linearScan(query_file);



            vector<int> db_data;

            string db_fn=globalDataFilename;
            ifstream db_file;
            db_file.open(db_fn.c_str());

    if (db_file.fail())
    {
        cout<<"cant open file "<<db_fn<<endl;
        exit(1);

    }




            for(int j = 0; j < db_size; j++)
            {

                //cout<<"query data"<<i<<"db data"<<j<<endl;
                db_data.clear();
                string line;
                getline(db_file, line);
                istringstream dbstr(line);
                for(int k = 0; k < DIM; k++)
                {
                    int n;
                    dbstr >> n;
                    db_data.push_back(n);
                }


                bool matchonalldim = true;        
                for(unsigned int s=0;(s<boxquerydata.size())&&matchonalldim ;s++)
                {
                    bool matchononedim=false;
                    for(unsigned int x=0;x<boxquerydata[s].size();x++)    
                        if(boxquerydata[s][x]==db_data[s])
                            matchononedim=true;

                    if(!matchononedim)
                        matchonalldim=false;

                }



                if(matchonalldim)
                {
                    bool findduplicate=false;
                    for(unsigned int i=0;i<matchedpoints.size();i++)
                        if(matchedpoints.at(i)==db_data)
                            findduplicate=true;

                    if(!findduplicate)
                    {
                        matchedpoints.push_back(db_data);
                        //for(int t1=0;t1<boxquerydata.size();t1++)
                        //{
                        //    for(int t2=0;t2<boxquerydata.at(t1).size();t2++)
                        //        cout<<boxquerydata.at(t1).at(t2)<<" ";
                        //    cout<<endl;
 
                        //}    

                        //for(int t2=0;t2<db_data.size();t2++)
                        //    cout<<db_data.at(t2)<<": ";
                        //cout<<endl;


                    }
                }


            } 
            db_file.close();
            totalnumofmatches+=(int)matchedpoints.size();

        }// end of         for(int i = 0; i < num_of_points; i++) // for every query point




        cout<<"linear scan boxsize="<<boxSize<<endl;
        cout <<  "matched data # from linear scan:" << static_cast<double>(totalnumofmatches) / num_of_points << endl;

    }
    query_file.close();



}



void ndtreeHelper::display_help()
{
    cout<<"this implementation of bond tree was modified by alok"<<endl;
    cout<<"it only works with discrete dimensions and performs random sized box queries."<<endl;
    cout<<"please feel free to contact alok at watvealo@cse.msu.edu for any comments"<<endl;
    cout<<"or suggestions"<<endl;
    cout<<"currently supported options are "<<endl;
    cout<<"1) "<<" : specifies the template for index filenames. as there"<<endl;
    cout<<"\tare multiple files created, each file will have this name as the prefix"<<endl;
    cout<<"and be appended with the suitable suffix (0,1,2...)"<<endl;
    cout<<"2) "<<" : specifies the data file to be loaded in the index"<<endl;
    cout<<"\tdatafile must be in csv format."<<endl;
    cout<<"3) "<<" : name of the box query file"<<endl;
//  cout<<"\t\tbox query is specified as x1,x2,x3:y1,y2:z1"<<endl;
//  cout<<"\t\t\there, x1..x3 are characters range along first dimension, y1-y2 is range along second"<<endl;
//  cout<<"dimension and so on"<<endl
    cout<<"4) "<<" : name of the range query file"<<endl;
    cout<<"5) "<<" : radius of the range query"<<endl;
    cout<<"6) "<<" : number of records in the data file that should be skipped"<<endl;
    cout<<"7) "<<" : number of data records to be loaded in the database"<<endl;
    cout<<"8) "<<" : flag indicating that a new index file should be created"<<endl;
    cout<<"\tany existing index file will be deleted"<<endl;
}

#endif
