#ifndef _FQTOOL
#define _FQTOOL
#define _out_buffer_ 1024*1024*4
#define Store_line_num 1024
//#define speedy_gzip_software "/annoroad/data1/bioinfo/PMO/malixiang/workdir/fqtools//speedy_gzip"
#include "cutadapt.h"
#include "hash.h"
#include "yarn.h"



struct STAT{
	unsigned int raw_reads;
	long long raw_bases;
	unsigned int clean_reads;
	long long clean_bases;
	float clean_reads_rate;
	unsigned int low_quality_reads;
	float low_quality_rate;
	unsigned int ns_reads;
	float Ns_rate;
	unsigned int adapter_reads;
	float adapter_rate;
	long long raw_Q30_bases;
	float raw_Q30_bases_rate;
	long long clean_Q30_bases;
	float clean_Q30_bases_rate;
	unsigned  int HighAT_reads;
	float HighAT_rate;
	unsigned  int HighGC_reads;
	float HighGC_rate;
	unsigned  int ployAT_reads;
	float ployAT_rate;
	unsigned int base1[1024][5];
	unsigned int base2[1024][5];
	unsigned int quality1[1024][50];
	unsigned int quality2[1024][50];
};

struct global{
	unsigned char single_or_pair;
	int CUT;
	int start;
	long QUOTA;
	int LOWQ;
	int HIGHQ;
	float LOW_QUALITY;
	float highN;
	float highAT;
	float highGC;
	int ployAT;
	int ASCII;
	int Nthread;
	char outraw;
	//----------------------------
	long read1Count; 
	long read2Count;
	//----------adapter-----------
	char* sequence1;
	char* sequence2;
	//float error_rate;
	int min_overlap;
	//adapter_* adapter;
	
	hashtable* table1;
	hashtable* table2;
	//--------out file name----------
	char* fq1; 
	char* fq2;
	char* outfq1;
	char* outfq2;
	char* rawcutfq1;
	char* rawcutfq2;
	char* adapter1;
	char* adapter2;
	gzFile* gzfp1;
	gzFile* gzfp2;
	FILE* fpout1;
	FILE* fpout2;
	FILE* fpout1_raw;
	FILE* fpout2_raw;
	char* speedy_gzip;
};


struct readline{
	int seqlen;
	int Q30;
	int A;
	int T;
	int GC;
	int low;
	int N;
	char highAT_sign;
	char reads[4][1024];
	char id[1024];
	adapter_* adapter;
};

struct readclass{
	struct readline* read1;
	struct readline* read2;
	struct STAT* stat; 
};

struct space{
	lock* use;
	struct readclass* readclass;
	struct space* next;
	struct pool* pool;
};


struct pool{
	lock* have;
	struct space* head;
	int limit;
	char term;
};


struct job {
    long seq;                   // sequence number
    int more;                   // true if this is not the last chunk
    struct space* in;
    lock* calc;                 // filter lock
    struct job* next;           // next job in the list (either list)
    struct job* next1;
    char* buffer1;
	int l1;
    char* buffer2;
    int l2;
    char* buffer1_raw;
    int l1_r;
    char* buffer2_raw;
    int l2_r;
};




#endif
