#ifndef __ALIGN
#define __ALIGN



typedef struct Entry_{
	int cost;
	int matches;
	int origin;
}Entry_;

typedef struct result_{
	int start1;
	int stop1;
	int start2;
	int stop2;
	int matches;
	int errors;
}result_;

typedef struct adapter_{
	Entry_* columns;
	result_* result;
	int m;
	int n;
	int min_overlap;
	char s1[512];
	char* s2;
}adapter_;

typedef struct adaptline_{
	adapter_* adapter;
	gzFile* gzfp;
	char* id;
}adaptline_;


adapter_* initalign(const char* s1, const int min_overlap);
char match(adapter_* adapter);
void cutadapt_main(char* fq1, char* adapter_seq, int min_overlap);



#endif