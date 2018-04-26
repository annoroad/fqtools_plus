#ifndef DUP_STAT
#define DUP_STAT
#define KMER 10
#define MAX 100
#define UINT_MAX_ 0xffffffff
#include "rb_tree.h"


typedef struct mat_unit_{
	int index;
	unsigned int count;
	char sign;
//	unit_* units;
}mat_unit_;

typedef struct DupStat_{
	mat_unit_** matrix_stat;
	RBroot* RB;
	Node* MinNode;
	gzFile* gzfp1;
	gzFile* gzfp2;
	char* seq1;
	char* seq2;
	unsigned int* stat;
	unsigned int statlen;
	int len1;
	int len2;
	int len;
	int ncol;
	int nrow; 
	int allreadnum;
	int TreeElementNum; //当前树的元素数目
}DupStat_;


#endif
