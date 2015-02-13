#include <getopt.h>
#include "ndtreeHelper.h"

struct option longopts[] = {
    { "dimension", required_argument, NULL,'k'},
    { "index", required_argument, NULL,'i'},
    { "data", required_argument, NULL,'d'},
    { "boxquery", required_argument, NULL,'b'},
    { "querydim", required_argument, NULL,'q'},
    { "rangequery", required_argument, NULL,'r'},
    { "aux", required_argument, NULL,'a'},
    { "record", required_argument, NULL,'c'},
    { "mode", required_argument, NULL,'m'},
    { "help", no_argument, NULL,'h'},
   {     0,    0,    0,    0},
};

//enum Action = {};

int main(int argc, char *argv[])
{
	ndtreeHelper ndtree;
	int query_dim = DIM;
	bool newTree = false;
	bool isInsert = false;
	bool isBoxQuery = false;
	int c;
    	while((c = getopt_long(argc, argv, 
		"k:i:d:b:q:r:s:a:c:m:h", longopts, NULL)) != -1){

		switch (c){
		case 'k':
		    printf("dimension is: %s\n", optarg);
		    break;
		case 'i':
		    printf("index file name is: %s\n", optarg);
		    globalIndexFilename = optarg;
		    break;
		case 'd':
		    printf("data file name is %s.\n", optarg);
		    globalDataFilename = optarg;
		    isInsert = true;
		    break;
		case 'b':
		    printf("box query file name is %s.\n", optarg);
		    globalBQFilename = optarg;
		    isBoxQuery = true;
		    break;
		case 'q':
		    printf("query dimension is %s.\n", optarg);
		    query_dim = atoi(optarg);
		    break;
		case 'r':
		    printf("range query file name is %s.\n", optarg);
		    globalRQFilename = optarg;
		    break;
		case 's':
		    printf("skip is %s.\n", optarg);
		    break;
		case 'a':
		    printf("aux file name is %s.\n", optarg);
		    globalAuxFilename = optarg;
		    break;
		case 'c':
		    printf("record file name is %s.\n", optarg);
		    globalRecordFilename = optarg;
		    break;
		case 'm':
		    if(!strcmp(optarg, "rebuild")){
			cout << "rebuild mode" << endl;
			newTree = true;

		    } else if(!strcmp(optarg, "append")) {
			cout << "append mode" << endl;
			newTree = false;
		    } else if(!strcmp(optarg, "default")) {
			cout << "deault mode" << endl;
			newTree = false;
		    } else {
			cout << "unkown mode" << endl;
			newTree = false;
		    }
		    break;
		case 'h':
		    printf("\nHelp information \n");
		    break;
		}
	}

	if(isInsert) {
		cout << "data file" << globalDataFilename << endl; 
	        ndtree.batchBuild_with_duplicate_record(newTree, LONG_MAX);
        }

	if( isBoxQuery ) {
		cout<<"Box query file "<<globalBQFilename<<endl;
		ndtree.batchRandomBoxQuery(query_dim);
	}
	
    
}
