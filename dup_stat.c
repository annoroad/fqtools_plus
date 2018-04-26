#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <zlib.h>
#include "dup_stat.h"
#include "rb_tree.h"

extern RBroot* RB_create(void);
extern Node* RB_get(RBroot* root, unsigned int key);
extern void RB_insert(RBroot* root, unsigned int key, char* value);
extern int RB_delete(RBroot* root, unsigned int key);
extern Node* RB_min(RBroot* root);
extern queue* RB_keys(RBroot* root);
extern void RB_clean(RBroot* root);

static unsigned int kmer(char* seq);
static void stat_freq(DupStat_* dupclass);
static int func(const void* a, const void* b);
static void find_high_freq(DupStat_* dupclass);
static void seq_kmer(DupStat_* dupclass);
static void catseq(DupStat_* dupclass);
static char readfq(DupStat_* dupclass);
static DupStat_* init_dup(const char* fq1, const char* fq2);
static void mainfunc(const char* fq1 , const char* fq2);
static void finalstat(DupStat_* dupclass, char* outfile);

static const char* HELP="\
Usage:\n\
dup_stat <input fq1.gz> <input fq2.gz>\n\
";

static unsigned int kmer(char* seq){
	//A 0 T 1 G 2 C 3
	char* cur = seq;
	unsigned int i = 0;
	unsigned int j;
	unsigned int value = 0;
	while(i<KMER && *cur){
		switch (*cur){
			case 'A':
				j=i<<1;
				value += 0<<j;
				break;
			case 'a':
				j=i<<1;
				value += 0<<j;
				break;
			case 'T':
				j=i<<1;
				value += 1<<j;
				break;
			case 't':
				j=i<<1;
				value += 1<<j;
				break;
			case 'G':
				j=i<<1;
				value += 2<<j;
				break;
			case 'g':
				j=i<<1;
				value += 2<<j;
				break;
			case 'C':
				j=i<<1;
				value += 3<<j;
				break;
			case 'c':
				j=i<<1;
				value += 3<<j;
				break;
			case 'N':
				j=i<<1;
				value += 0<<j;
				break;
			case 'n':
				j=i<<1;
				value += 0<<j;
				break;
		}
		cur++;
		i++;
	}
	return(value);
}


static void stat_freq(DupStat_* dupclass){
	char* seq = dupclass->seq1;
	mat_unit_** matrix_stat = dupclass->matrix_stat;
	char curseq[KMER];
	int i;
	unsigned int KmerValue;
	unsigned int mincount=UINT_MAX_;
	char sign = 0;
	mat_unit_* curunit = NULL;
	for(i=0;i<dupclass->nrow;i++){
		memset(curseq,0,KMER);
		memcpy(curseq,seq,KMER);
		seq+=KMER;
		KmerValue = kmer(curseq);
		curunit = matrix_stat[i]+KmerValue;
		if (curunit->count < mincount){
			mincount = curunit->count;
			sign = curunit->sign;
		}
	}

	if (mincount > 1){
		if (dupclass->statlen < mincount){
			dupclass->stat = (unsigned int*)realloc(dupclass->stat,dupclass->statlen*2);
			dupclass->statlen*=2;
		}
		(dupclass->stat)[mincount]++;
	}

	if (sign){
		if (dupclass->TreeElementNum == MAX){
			if(dupclass->MinNode->key < mincount){
				if(RB_delete(dupclass->RB, dupclass->MinNode->key)){
					RB_insert(dupclass->RB, mincount, dupclass->seq1);
					dupclass->MinNode = RB_min(dupclass->RB);	
				}
			}
		}else{
			RB_insert(dupclass->RB, mincount, dupclass->seq1);
			dupclass->MinNode = RB_min(dupclass->RB);
			(dupclass->TreeElementNum)++;	
		}
	}
}

static int func(const void* a, const void* b){
	return((*(mat_unit_**)b)->count - (*(mat_unit_**)a)->count);
}
static void find_high_freq(DupStat_* dupclass){
	int i,j,m;
	mat_unit_** matrix_stat = dupclass->matrix_stat;
	mat_unit_*** matrix_final = (mat_unit_***)malloc(sizeof(mat_unit_**)*dupclass->nrow);
	for(i=0;i<dupclass->nrow;i++){
		matrix_final[i] = (mat_unit_**)malloc(sizeof(mat_unit_*)*dupclass->ncol);
		for (j=0;j<dupclass->ncol;j++){
			matrix_final[i][j]=&matrix_stat[i][j];
		}
		qsort(matrix_final[i],dupclass->ncol,sizeof(mat_unit_*),func);
		for (m=0;m<MAX;m++){
			matrix_final[i][m]->sign = 1;
		}
		free(matrix_final[i]);
	}
	free(matrix_final);
}


static void seq_kmer(DupStat_* dupclass){
	char* seq = dupclass->seq1;
	mat_unit_** matrix_stat = dupclass->matrix_stat;
	char curseq[KMER];
	int i;
	unsigned int KmerValue;
	for(i=0;i<dupclass->nrow;i++){
		memset(curseq,0,KMER);
		memcpy(curseq,seq,KMER);
		seq+=KMER;
		KmerValue = kmer(curseq);
		(matrix_stat[i][KmerValue].count)++;
	}
}

static void catseq(DupStat_* dupclass){
	char* cur=dupclass->seq1;
	char* cur1=dupclass->seq2;
	while(*cur++){}
	cur--;
	while(*cur1){
		*cur++ = *cur1++;
	}
	*cur='\0';
}

static char readfq(DupStat_* dupclass){
	char line[1024];
	if(gzgets(dupclass->gzfp1,line,1024)!=NULL){}
	else return(0);
	if(gzgets(dupclass->gzfp1,dupclass->seq1,1024)!=NULL) {
		(dupclass->seq1)[dupclass->len1]='\0';
	}else return(0);
	if(gzgets(dupclass->gzfp1,line,1024)!=NULL){}
	else return(0);
	if(gzgets(dupclass->gzfp1,line,1024)!=NULL){}
	else return(0);
	if (dupclass->gzfp2){
		if(gzgets(dupclass->gzfp2,line,1024)!=NULL){}
		else return(0);
		if(gzgets(dupclass->gzfp2,dupclass->seq2,1024)!=NULL) {
			(dupclass->seq2)[dupclass->len2]='\0';
		}else return(0);
		if(gzgets(dupclass->gzfp2,line,1024)!=NULL){}
		else return(0);
		if(gzgets(dupclass->gzfp2,line,1024)!=NULL){}
		else return(0);
	}
	return(1);
}


static DupStat_* init_dup(const char* fq1, const char* fq2){
	gzFile* gzfp1=gzopen(fq1,"rb");
	gzFile* gzfp2=NULL;
	char line[1024];
	if (*fq2){
		gzfp2=gzopen(fq2,"rb");
	}
	DupStat_* dupclass = (DupStat_*)malloc(sizeof(DupStat_));
	dupclass->gzfp1 = gzfp1;
	dupclass->gzfp2 = gzfp2;
	dupclass->seq1 = (char*)malloc(sizeof(char)*1024);
	dupclass->seq2 = (char*)malloc(sizeof(char)*1024);
	dupclass->len1=0;
	dupclass->len2=0;
	dupclass->allreadnum=0;
	gzgets(dupclass->gzfp1,line,1024);
	gzgets(dupclass->gzfp1,dupclass->seq1,1024);
	dupclass->len1 = strlen(dupclass->seq1)-1;
	gzseek(dupclass->gzfp1,0L,SEEK_SET);
	if (dupclass->gzfp2){
		gzgets(dupclass->gzfp2,line,1024);
		gzgets(dupclass->gzfp2,dupclass->seq2,1024);
		dupclass->len2 = strlen(dupclass->seq2)-1;
		gzseek(dupclass->gzfp2,0L,SEEK_SET);
	}
	dupclass->len = dupclass->len1 + dupclass->len2;
	int nrow = (int)ceil((double)dupclass->len/KMER);
	int ncol = 2<<((KMER<<1)-1); 
	dupclass->ncol=ncol;
	dupclass->nrow=nrow;
	dupclass->stat = (unsigned int*)malloc(sizeof(unsigned int)*5242880);
	memset(dupclass->stat,0,sizeof(unsigned int)*5242880);
	dupclass->statlen = 5242880-1;
	dupclass->matrix_stat = (mat_unit_**)malloc(sizeof(mat_unit_*)*nrow);
	int i,j;
	for(i=0;i<nrow;i++){
		(dupclass->matrix_stat)[i] = (mat_unit_*)malloc(sizeof(mat_unit_)*ncol);
		for (j=0;j<ncol;j++){
			(dupclass->matrix_stat)[i][j].index = j;
			(dupclass->matrix_stat)[i][j].count = 0;
			(dupclass->matrix_stat)[i][j].sign = 0;
		}
	}
	dupclass->RB = RB_create();
	dupclass->MinNode = NULL;
	dupclass->TreeElementNum = 0;
	return(dupclass);
}


static void mainfunc(const char* fq1 , const char* fq2){
	DupStat_* dupclass = init_dup(fq1,fq2);
	while(1){
		if(readfq(dupclass)==0) break;
		if(dupclass->gzfp2) catseq(dupclass);
		seq_kmer(dupclass);
		(dupclass->allreadnum)++;
	}

	find_high_freq(dupclass);
	gzclose(dupclass->gzfp1);
	dupclass->gzfp1 = gzopen(fq1,"rb");
	if (dupclass->gzfp2) {
		gzclose(dupclass->gzfp2);
		dupclass->gzfp2 = gzopen(fq2,"rb");
	}
	while(1){
		if(readfq(dupclass)==0) break;
		if(dupclass->gzfp2) catseq(dupclass);
		stat_freq(dupclass);
	}

	char* outfile = (char*)malloc(strlen(fq1)+10);
	sprintf(outfile,"%s.dup.stat",fq1);
	finalstat(dupclass, outfile);

}

static void finalstat(DupStat_* dupclass, char* outfile){
	FILE* fp=fopen(outfile,"w");
	int i;
	(dupclass->statlen)++;
	unsigned int DupAllCount=0;
	for(i=0;i<dupclass->statlen;i++){
		if ((dupclass->stat)[i]!=0){
			printf("dup:%d num:%u\n",i,(dupclass->stat)[i]);
			DupAllCount+=(unsigned int)((dupclass->stat)[i] - (dupclass->stat)[i]/(unsigned int)i);
		}
	}
	fprintf(fp, "Total Reads:\t%u\n",dupclass->allreadnum);
	fprintf(fp, "Dup Reads:\t%u\n",DupAllCount);
	fprintf(fp, "Dup Rate(%%):\t%.3f\n",100.0*DupAllCount/dupclass->allreadnum);
	fprintf(fp, "#Number\tReads\tcount\tRate(%%)\n");
	queue* curqueue = RB_keys(dupclass->RB);
	queue* tmpqueue = NULL;
	unsigned int curcount;
	char* curseq = NULL;
	i=1;
	while(curqueue){
		curcount = curqueue->node->key;
		curseq = curqueue->node->value;
		fprintf(fp, "%d\t%s\t%u\t%.3f\n",i,curseq,curcount,100.0*curcount/dupclass->allreadnum);
		tmpqueue = curqueue;
		curqueue = curqueue->after;
		free(tmpqueue);
		i++;
	}
}


int main(int argc, char** argv){
	time_t start,end;
	time(&start);
	fprintf(stdout,"%s\n", ctime(&start));
	if(argc<2) {fprintf(stderr, "%s\n",HELP);exit(0);}

	char fq1[1024];
	char fq2[1024];
	memset(fq1,0,1024);
	memset(fq2,0,1024);
	strcpy(fq1,argv[1]);
	if (argc==3) strcpy(fq2,argv[2]);
	if(*fq1=='\0') {fprintf(stderr, "%s\n",HELP);exit(0);}
	mainfunc(fq1,fq2);
	
	
	time(&end);
	fprintf(stdout,"%s\n", ctime(&end));

	return(0);
}