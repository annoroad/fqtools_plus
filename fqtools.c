#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <time.h>
#include <regex.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/wait.h>
#include "fqtools.h"

#define mode O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH
#define local static
//Author lixiangma
//email lixiangma@genome.cn

local const char* HELP="\
Usage:\n\
[options] \n\
	filter    filter reads\n\
	cutadapt  cutadapt\n\
";

local const char* HELP1="\
Usage:\n\
fqtools_plus filter [options]  <input fq1.gz> <input fq2.gz> <adapter_r1> <adapter_r2> <out clean fq1.gz> <out clean fq2.gz> \n \
[options] \n\
    --c       INT     cut sequence(截取片段长度), default [0 not cut sequence]\n\
    --start   INT     start site(起始位点，从0开始), default [0]\n\
    --q       INT     cut data cutoff(截取数据量大小), default [0 not cut data]   \n\
    --ql      INT     base quality lower limit, default [19]\n\
    --qh      INT     base quality higher limit, default [30]\n\
    --Q       FLOAT   low-quality base rate limit, default [0.5]\n\
    --n       FLOAT   N rate limit, default [0.05]\n\
    --AT      FLOAT   single pair read High AT content, default [0.8]\n\
    --GC      FLOAT   single pair read High GC content, default [0.8]\n\
    --ployAT  INT     ploy A/T length(多聚AT的长度), default [0]\n\
    --i       INT     ASCII value stands for qulity 0, default [33] \n\
    --m       INT     num thread, default 1.当使用多线程时，为了更好的性能：使用数据量截取参数时，必须>=4线程，否则默认为非多线程运行。不使用数据量截取参数时，必须>=3线程，否则默认为非多线程运行。\n\
    --adapt1   string  adapt r1 sequence, if use this parameter, not use <adapter_r1> <adapter_r2>\n\
    --adapt2   string  adapt r2 sequence, if use this parameter, not use <adapter_r1> <adapter_r2>\n\
    -min_overlap INT    default 5 Minimum overlap length. If the overlap between the read and the adapter is shorter than LENGTH, the read is not modified.\n\
    --zip     string  zip software, default gzip\n\
\n\
[example]\n\
fqtools_plus filter [options] test_R1.fq.gz test_R2.fq.gz test_R1.adapter.txt.gz test_R2.adapter.txt.gz test_R1.clean.fq.gz test_R2.clean.fq.gz \n\
fqtools_plus filter [options] --adapt1 GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapt2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA test_R1.fq.gz test_R2.fq.gz test_R1.clean.fq.gz test_R2.clean.fq.gz\n\
";

local const char* HELP2="\
Usage:\n\
filter cutadapt [options]  <input fq.gz>\n\
[options] \n\
    --adapt       string  adapt sequence\n\
    --min_overlap INT   default 5   Minimum overlap length. If the overlap between the read and the adapter is shorter than LENGTH, the read is not modified.\n\
[example]\n\
fqtools_plus cutadapt --adapt1 GATCGGAAGAGCACACGTCTGTACTCCAGTCAC test_R1.fq.gz\n\
";


#define fail() {fprintf(stderr,"wrong\n");;exit(0);}

extern hashtable* init_hash(void);
extern void insert_keyvalue(hashtable** ht, char* key, int value);
extern char get_index(hashtable* ht, char* key);
extern adapter_* initalign(const char* s1, const int min_overlap);
extern char match(adapter_* adapter);
extern void cutadapt_main(char* fq1, char* adapter_seq, int min_overlap);

typedef struct pool pool;
typedef struct space space;
typedef struct job job;
typedef struct readclass readclass;
typedef struct readline readline;
typedef struct STAT STAT;

struct STAT allstat;
struct global g;

// List of compress jobs (with tail for appending to list).
local lock *compress_have = NULL;   // number of compress jobs waiting
local struct job *compress_head, **compress_tail;

// List of write jobs.
local lock *write_first;            // lowest sequence number in list
local struct job *write_head;

// raw of write jobs
local lock *write_first_raw;
local struct job *write_head_raw;

struct pool in_pool;

local void merge_stat(struct pool* pool);
local int readn(struct space* space);
local void to_buffer_raw(struct readclass* readclass, struct job* job);
local inline void substr_(const char* line, char* str);
local void raw_Q30(struct readclass* readclass);
local void cutseq(struct readclass* readclass);
local void to_buffer(struct readclass* readclass, struct job* job);
local char checkadapter(struct readclass* readclass);
local void quality_base_percent(struct readclass* readclass);
local char check_ployAT(char* seq, int ployAT);
local void report(struct readclass* readclass);

local void new_pool(pool* pool, int limit){
	pool->have = new_lock(0);
	pool->head = NULL;
	pool->limit = limit;
	pool->term = 1;
}


local space* new_space(struct pool* pool){
	space* space = (struct space*)malloc(sizeof(struct space));
	space->use = new_lock(1);
	space->pool = pool;
	space->next = NULL;
	space->readclass = (struct readclass*)malloc(sizeof(struct readclass));
	space->readclass->read1 = NULL;
	space->readclass->read2 = NULL;
	space->readclass->stat = (struct STAT*)malloc(sizeof(struct STAT));
	memset(space->readclass->stat,0,sizeof(struct STAT));
	space->readclass->read1 = (struct readline*)malloc(sizeof(struct readline));
	memset(space->readclass->read1,0,sizeof(readline));
	space->readclass->read1->adapter = NULL;
	if (*g.sequence1) space->readclass->read1->adapter = initalign(g.sequence1, g.min_overlap);
	if (g.single_or_pair){
		space->readclass->read2 = (struct readline*)malloc(sizeof(struct readline));
		memset(space->readclass->read2,0,sizeof(struct readline));
		space->readclass->read2->adapter = NULL;
		if (*g.sequence2) space->readclass->read2->adapter = initalign(g.sequence2, g.min_overlap);
	}
	return(space);
}

local space* get_space(struct pool* pool){
	struct space* space;
	possess(pool->have);
	
	if (pool->limit == 0){
		wait_for(pool->have, NOT_TO_BE, 0);
	}
	
	if (pool->head != NULL){
		space = pool->head;
		possess(space->use);
		pool->head = space->next;
		twist(pool->have,BY,-1);
		twist(space->use,TO,1);
		return(space);
	}
	if (pool->limit > 0){
		pool->limit--;
	}
	release(pool->have);
	space = new_space(pool);
	return(space);
}


//local void use_space(struct space* space){
//	possess(space->use);
//	twist(space->use,BY,+1);
//}

local void drop_space(struct space* space){
	int use;
	struct pool* pool;
	if (space == NULL) return;
	possess(space->use);
	use = peek_lock(space->use);
	if (use == 1){
		pool = space->pool;
		possess(pool->have);
		space->next = pool->head;
		pool->head = space;
		twist(pool->have,BY,+1);
	}
	twist(space->use,BY,-1);
}


local void free_pool(struct pool* pool){
	int count = 0;
	struct space* space;
	possess(pool->have);
	while((space = pool->head) != NULL){
		pool->head = space->next;
		free_lock(space->use);
		free(space->readclass->stat);
		free(space->readclass->read1->adapter);
		if (g.single_or_pair) free(space->readclass->read2->adapter);
		free(space->readclass->read1);
		if (g.single_or_pair) free(space->readclass->read2);
		free(space);
		count++;
	}
	release(pool->have);
	free_lock(pool->have);
}


local void setup_jobs(void) {
	if (compress_have != NULL)
		return;
    // allocate locks and initialize lists
    compress_have = new_lock(0);
    compress_head = NULL;
    compress_tail = &compress_head;
    write_first = new_lock(-1);
    write_head = NULL;
    write_first_raw = new_lock(-1);
    write_head_raw = NULL;
    if (g.outraw) g.Nthread-=3;
    else g.Nthread-=2;
    new_pool(&in_pool, g.Nthread);
}

local void finish_jobs(void) {
    struct job job;
    // only do this once
    if (compress_have == NULL)
        return;

    // command all of the extant compress threads to return
    possess(compress_have);
    job.seq = -1;
    job.next = NULL;
    job.next1 = NULL;
    compress_head = &job;
    compress_tail = &(job.next);
    twist(compress_have, BY, +1);       // will wake them all up
    // join all of the compress threads, verify they all came back
    join_all();
    // free the resources
    merge_stat(&in_pool);
    free_pool(&in_pool);
    free_lock(write_first);
    free_lock(compress_have);
    if (g.outraw) free_lock(write_first_raw);
    compress_have = NULL;
    gzclose(g.gzfp1);
    if (g.single_or_pair) gzclose(g.gzfp2);
    fclose(g.fpout1);
    if (g.single_or_pair) fclose(g.fpout2);
    if (g.outraw){
    	fclose(g.fpout1_raw);
	    if (g.single_or_pair) fclose(g.fpout2_raw);
    }
}


local job* new_job(void){
	job* job = (struct job*)malloc(sizeof(struct job));
	job->next = NULL;
	job->next1 = NULL;
	job->buffer1 = NULL;
	job->buffer2 = NULL;
	job->buffer1_raw = NULL;
	job->buffer2_raw = NULL;
	if (g.outraw){
		job->calc = new_lock(0);
	}
	return(job);
}

local void free_job(struct job* job){
	if (job == NULL) return;
	if (g.outraw){
		possess(job->calc);
		wait_for(job->calc,TO,1);
		release(job->calc);
		free_lock(job->calc);
	}
	if (job->buffer1 != NULL) free(job->buffer1);
	if (g.single_or_pair) if (job->buffer2 != NULL) free(job->buffer2);
	if (job->buffer1_raw != NULL) free(job->buffer1_raw);
	if (g.single_or_pair) if (job->buffer2_raw != NULL) free(job->buffer2_raw);
	free(job);
}


local void push_job_raw(struct job* job){
	struct job *here, **prior;
	possess(write_first_raw);
	prior = &(write_head_raw);
	while((here = *prior) != NULL){
		if (here->seq > job->seq){
			break;
		}
		prior = &(here->next1);
	}
	job->next1 = here;
	*prior = job;
	twist(write_first_raw,TO,write_head_raw->seq);
}

local void push_job(struct job* job){
	struct job *here, **prior;
	possess(write_first);
	prior = &(write_head);
	while((here = *prior) != NULL){
		if (here->seq > job->seq){
			break;
		}
		prior = &(here->next);
	}
	job->next = here;
	*prior = job;
	twist(write_first,TO,write_head->seq);
}



void writen(void){
	struct job* job = NULL;
	long seq = 0;
	long clean_seq = 0;
	FILE* fpout1 = g.fpout1;
	FILE* fpout2 = g.fpout2;
	for(;;){
		possess(write_first);
		wait_for(write_first,TO_BE,seq);
		job = write_head;
		if (clean_seq == g.QUOTA) {
			possess(in_pool.have);
			in_pool.term = 0;
			release(in_pool.have);
		}
		drop_space(job->in);
		write_head = job->next;
		twist(write_first,TO,write_head == NULL? -1:write_head->seq);
		
		seq++;
		if (!job->more) break;
		if (job->buffer1 == NULL){
			free_job(job);
			continue;
		}

		fwrite(job->buffer1,sizeof(char),job->l1,fpout1);
		if (g.single_or_pair) fwrite(job->buffer2,sizeof(char),job->l2,fpout2);
		free_job(job);
		clean_seq++;
	}

	free_job(job);

	possess(compress_have);
	if (compress_head->more != 0)
		fail();
	release(compress_have);
	possess(write_first);
	if (write_head != NULL)
		fail();
	twist(write_first,TO,-1);
}


void writen_raw(void){
	struct job* job;
	long seq = 0;
	FILE* fpout1_raw = g.fpout1_raw;
	FILE* fpout2_raw = g.fpout2_raw;
	for(;;){
		possess(write_first_raw);
		wait_for(write_first_raw,TO_BE,seq);
		job = write_head_raw;
		write_head_raw = job->next1;
		twist(write_first_raw,TO,write_head_raw == NULL? -1:write_head_raw->seq);
		seq++;
		if (!job->more) break;
		fwrite(job->buffer1_raw,sizeof(char),job->l1_r,fpout1_raw);
		if (g.single_or_pair) fwrite(job->buffer2_raw,sizeof(char),job->l2_r,fpout2_raw);
		possess(job->calc);
		twist(job->calc,TO,1);
	}

	possess(job->calc);
	twist(job->calc,TO,1);

	possess(write_first_raw);
	if (write_head_raw != NULL)
		fail();
	twist(write_first_raw,TO,-1);
}


static void filter_pair(void){
	struct job* job;
	struct readclass* readclass;
	struct readline* read1=NULL;
	struct readline* read2=NULL;
	struct STAT* stat;
	for (;;){
		possess(compress_have);
		wait_for(compress_have,NOT_TO_BE,0);
		job = compress_head;
		//job != NULL
		if (job->seq == -1) break;
		if (job->next == NULL){
			compress_tail = &compress_head;
		}else{
			compress_head = job->next;
		}
		twist(compress_have,BY,-1);
		
		readclass = job->in->readclass;
		read1 = readclass->read1;
		if (g.single_or_pair) read2 = readclass->read2;
		stat = readclass->stat;
		//filter funtion---------------------------------------
		if (g.outraw){
			to_buffer_raw(readclass,job);
			push_job_raw(job);
		}

		if (!job->more) {
			push_job(job);
			continue;
		}

		if (g.single_or_pair){
			read1->seqlen = strlen(read1->reads[1])-1;
			read2->seqlen = strlen(read2->reads[1])-1;
			substr_((read1->reads)[0], read1->id);
			substr_((read2->reads)[0], read2->id);
			if (read1->seqlen < g.CUT || read2->seqlen < g.CUT) {
				push_job(job);
				continue;
			}
			raw_Q30(readclass);
			stat->raw_reads++;
			stat->raw_bases+= read1->seqlen + read2->seqlen;
			cutseq(readclass);
			read1->seqlen = g.CUT;
			read2->seqlen = g.CUT;
			if (!checkadapter(readclass)){
				stat->adapter_reads++;
				push_job(job);
				continue;
			}
			quality_base_percent(readclass);
			if (read1->low >= g.LOW_QUALITY || read2->low >= g.LOW_QUALITY){
				stat->low_quality_reads++;
				push_job(job);
				continue;
			}
			if (read1->N >= g.highN || read2->N >= g.highN){
			    stat->ns_reads++;
			  	push_job(job);
				continue;
			}
			if (check_ployAT((read1->reads)[1], g.ployAT) || check_ployAT((read2->reads)[1], g.ployAT)){
			  	stat->ployAT_reads++;
			}
			if (read1->GC >= g.highGC || read2->GC >= g.highGC){
			   	stat->HighGC_reads++;
			}
			if (read1->highAT_sign || read2->highAT_sign){
			   	stat->HighAT_reads++;
			}
			stat->clean_Q30_bases += read1->Q30 + read2->Q30;
			stat->clean_reads++;
			report(readclass);
		}else{
			read1->seqlen = strlen(read1->reads[1])-1;
			substr_((read1->reads)[0], read1->id);
			if (read1->seqlen < g.CUT) {
				push_job(job);
				continue;
			}
			raw_Q30(readclass);
			stat->raw_reads++;
			stat->raw_bases+= read1->seqlen;
			cutseq(readclass);
			read1->seqlen = g.CUT;
			if (!checkadapter(readclass)){
				stat->adapter_reads++;
				push_job(job);
				continue;
			}
			quality_base_percent(readclass);
			if (read1->low >= g.LOW_QUALITY){
				stat->low_quality_reads++;
				push_job(job);
				continue;
			}
			if (read1->N >= g.highN){
			    stat->ns_reads++;
			  	push_job(job);
				continue;
			}

			if (check_ployAT((read1->reads)[1], g.ployAT)){
			  	stat->ployAT_reads++;
			}
			if (read1->GC >= g.highGC){
			   	stat->HighGC_reads++;
			}
			if (read1->highAT_sign){
			   	stat->HighAT_reads++;
			}
			stat->clean_Q30_bases += read1->Q30;
			stat->clean_reads++;
			report(readclass);
		}
		//-----------------------------------------------------
		to_buffer(readclass,job);
		push_job(job);
	}
	release(compress_have);
}

void parallel_filter_main(void){
	setup_jobs();
	struct job* job = NULL;
	int more;
	long seq = 0;
	int cthreads = 0;
	thread* writeth = launch((void*)writen,NULL);
	thread* writeth_raw = NULL;
	if (g.outraw) {
		writeth_raw= launch((void*)writen_raw,NULL);
	}
	if (g.Nthread <= 0) {
		fprintf(stderr, "Nthread is wrong!\n");
		exit(0);
	}

	do{
		job = new_job();
		job->in = get_space(&in_pool);
		more = readn(job->in);
		job->more = more;
		job->seq = seq;

		if (cthreads++ < g.Nthread){
			launch((void*)filter_pair,NULL);
		}
		seq++;
		possess(compress_have);
		job->next = NULL;
		*compress_tail = job;
		compress_tail = &(job->next);
		twist(compress_have,BY,+1);
	}while(more);

	join(writeth);
	if (writeth_raw != NULL) join(writeth_raw);
	finish_jobs();
}

void single_filter_main(void){
    new_pool(&in_pool, -1);
	struct job* job = new_job();
	job->in = get_space(&in_pool);
	in_pool.head = job->in;
	struct readclass* readclass = job->in->readclass;
	struct readline* read1 = readclass->read1;
	struct readline* read2 = readclass->read2;
	struct STAT* stat = readclass->stat;

	for(;;){
		if (!readn(job->in)) break;
		if (g.outraw){
			to_buffer_raw(readclass,job);
			fwrite(job->buffer1_raw,sizeof(char),job->l1_r,g.fpout1_raw);
			if (g.single_or_pair) fwrite(job->buffer2_raw,sizeof(char),job->l2_r,g.fpout2_raw);
		}

		if (g.single_or_pair){
			read1->seqlen = strlen(read1->reads[1])-1;
			read2->seqlen = strlen(read2->reads[1])-1;
			substr_((read1->reads)[0], read1->id);
			substr_((read2->reads)[0], read2->id);
			if (read1->seqlen < g.CUT || read2->seqlen < g.CUT) {
				continue;
			}
			raw_Q30(readclass);
			stat->raw_reads++;
			stat->raw_bases+= read1->seqlen + read2->seqlen;	
			cutseq(readclass);
			read1->seqlen = g.CUT;
			read2->seqlen = g.CUT;
			if (!checkadapter(readclass)){
				stat->adapter_reads++;
				continue;
			}
			quality_base_percent(readclass);
			if (read1->low >= g.LOW_QUALITY || read2->low >= g.LOW_QUALITY){
				stat->low_quality_reads++;
				continue;
			}
			if (read1->N >= g.highN || read2->N >= g.highN){
			    stat->ns_reads++;
				continue;
			}
			if (check_ployAT((read1->reads)[1], g.ployAT) || check_ployAT((read2->reads)[1], g.ployAT)){
			  	stat->ployAT_reads++;
			}
			if (read1->GC >= g.highGC || read2->GC >= g.highGC){
			   	stat->HighGC_reads++;
			}
			if (read1->highAT_sign || read2->highAT_sign){
			   	stat->HighAT_reads++;
			}
			stat->clean_Q30_bases += read1->Q30 + read2->Q30;
			stat->clean_reads++;
			report(readclass);
		}else{
			read1->seqlen = strlen(read1->reads[1])-1;
			substr_((read1->reads)[0], read1->id);
			if (read1->seqlen < g.CUT) {
				continue;
			}
			raw_Q30(readclass);
			stat->raw_reads++;
			stat->raw_bases+= read1->seqlen;
			cutseq(readclass);
			read1->seqlen = g.CUT;
			if (!checkadapter(readclass)){
				stat->adapter_reads++;
				continue;
			}
			quality_base_percent(readclass);
			if (read1->low >= g.LOW_QUALITY){
				stat->low_quality_reads++;
				continue;
			}
			if (read1->N >= g.highN){
			    stat->ns_reads++;
				continue;
			}

			if (check_ployAT((read1->reads)[1], g.ployAT)){
			  	stat->ployAT_reads++;
			}
			if (read1->GC >= g.highGC){
			   	stat->HighGC_reads++;
			}
			if (read1->highAT_sign){
			   	stat->HighAT_reads++;
			}
			stat->clean_Q30_bases += read1->Q30;
			stat->clean_reads++;
			report(readclass);

		}
		//-----------------------------------------------------
		to_buffer(readclass,job);
		fwrite(job->buffer1,sizeof(char),job->l1,g.fpout1);
		if (g.single_or_pair) fwrite(job->buffer2,sizeof(char),job->l2,g.fpout2);
		if (stat->clean_reads == g.QUOTA) break;
	}


	if (g.outraw){
		possess(job->calc);
		twist(job->calc,TO,1);	
	}

	gzclose(g.gzfp1);
    if (g.single_or_pair) gzclose(g.gzfp2);
    fclose(g.fpout1);
    if (g.single_or_pair) fclose(g.fpout2);
    if (g.outraw){
    	fclose(g.fpout1_raw);
	    if (g.single_or_pair) fclose(g.fpout2_raw);
    }

	if (job != NULL) free_job(job);
	merge_stat(&in_pool);
	free_pool(&in_pool);
}


local void merge_stat(struct pool* pool){
	struct space* space = pool->head;
	struct STAT* stat;
	memset(&allstat,0,sizeof(allstat));
	possess(pool->have);
	int m,n;
	while(space != NULL){
		stat = space->readclass->stat;
		allstat.raw_reads += stat->raw_reads;
		allstat.clean_reads += stat->clean_reads;
		allstat.low_quality_reads += stat->low_quality_reads;
		allstat.ns_reads += stat->ns_reads;
		allstat.adapter_reads += stat->adapter_reads;
		allstat.raw_bases += stat->raw_bases;
		allstat.raw_Q30_bases += stat->raw_Q30_bases;
		allstat.clean_Q30_bases += stat->clean_Q30_bases;
		allstat.HighAT_reads += stat->HighAT_reads;
		allstat.HighGC_reads += stat->HighGC_reads;
		allstat.ployAT_reads += stat->ployAT_reads;
		for(m=0;m<g.CUT;m++){
			for(n=0;n<5;n++){
				(allstat.base1)[m][n] += (stat->base1)[m][n];
				if (g.single_or_pair) (allstat.base2)[m][n] += (stat->base2)[m][n];
			}
			for(n=0;n<50;n++){
				(allstat.quality1)[m][n] += (stat->quality1)[m][n];
				if (g.single_or_pair) (allstat.quality2)[m][n] += (stat->quality2)[m][n];
			}
		}
		space = space->next;
	}
	release(pool->have);
	
}
//filter-----------------------------------------------------

local int readn(struct space* space){
	if (!space->pool->term) return(0);
	readline* read1 = space->readclass->read1;
	if (!gzgets(g.gzfp1,read1->reads[0],1024)) return(0);
	else g.read1Count++;
	if (!gzgets(g.gzfp1,read1->reads[1],1024))	return(0);
	if (!gzgets(g.gzfp1,read1->reads[2],1024))	return(0);
	if (!gzgets(g.gzfp1,read1->reads[3],1024))	return(0);

	if (g.single_or_pair){
		readline* read2 = space->readclass->read2;
		if (!gzgets(g.gzfp2,read2->reads[0],1024)) return(0);
		else g.read2Count++;
		if (!gzgets(g.gzfp2,read2->reads[1],1024))	return(0);
		if (!gzgets(g.gzfp2,read2->reads[2],1024))	return(0);
		if (!gzgets(g.gzfp2,read2->reads[3],1024))	return(0);
	}
	return(1);
}


local void to_buffer_raw(struct readclass* readclass, struct job* job){
	if (job->buffer1_raw == NULL) {job->buffer1_raw = malloc(1024*4);memset(job->buffer1_raw,0,1024*4);}
	char* cur1 = job->buffer1_raw;
	char (*reads1)[1024] = readclass->read1->reads;
	int i,j;
	for (i=0;i<4;i++){
		j=0;
		while( (*cur1++ = reads1[i][j++]) != '\n' ){}
	}
	*cur1 = '\0';
	job->l1_r = strlen(job->buffer1_raw);
	if (g.single_or_pair){
		if (job->buffer2_raw == NULL) {job->buffer2_raw = malloc(1024*4);memset(job->buffer2_raw,0,1024*4);}
		char* cur2 = job->buffer2_raw;
		char (*reads2)[1024] = readclass->read2->reads;
		for (i=0;i<4;i++){
			j=0;
			while( (*cur2++ = reads2[i][j++]) != '\n' ){}
		}
		*cur2 = '\0';
		job->l2_r = strlen(job->buffer2_raw);
	}
	
}


local inline void substr_(const char* line, char* str){
	while(*line++ != '\n'){
		*str = *line;
		str++;
	}
	*str = '\0';
}


local void raw_Q30(struct readclass* readclass){
	int i=0;
	int curqual;
	readline* read1 = readclass->read1;
	
	STAT* stat = readclass->stat;
	char* qual=(read1->reads)[3];
	while (i < read1->seqlen){
		curqual=(int)*(qual + i) - g.ASCII;
		if (curqual >= g.HIGHQ){
			stat->raw_Q30_bases++;
		}
		i++;
	}
	if (g.single_or_pair){
		readline* read2 = readclass->read2;
		i = 0;
		qual=(read2->reads)[3];
		while (i < read2->seqlen){
			curqual=(int)*(qual + i) - g.ASCII;
			if (curqual >= g.HIGHQ){
				stat->raw_Q30_bases++;
			}
			i++;
		}
	}
		
	
}

local void cutseq(struct readclass* readclass){
	readline* read1 = readclass->read1;
	memcpy((read1->reads)[1],(read1->reads)[1] + g.start,g.CUT);
	*((read1->reads)[1] + g.CUT) = '\n';
	*((read1->reads)[1] + g.CUT + 1) = '\0';
	memcpy((read1->reads)[3],(read1->reads)[3] + g.start,g.CUT);
	*((read1->reads)[3] + g.CUT) = '\n';
	*((read1->reads)[3] + g.CUT + 1) = '\0';
	if (g.single_or_pair){
		readline* read2 = readclass->read2;
		memcpy((read2->reads)[1],(read2->reads)[1],g.CUT);
		*((read2->reads)[1] + g.CUT) = '\n';
		*((read2->reads)[1] + g.CUT + 1) = '\0';
		memcpy((read2->reads)[3],(read2->reads)[3],g.CUT);
		*((read2->reads)[3] + g.CUT) = '\n';
		*((read2->reads)[3] + g.CUT + 1) = '\0';
	}
}

local void to_buffer(struct readclass* readclass, struct job* job){
	if (job->buffer1 == NULL) job->buffer1 = malloc(1024*4);
	char* cur1 = job->buffer1;
	char (*reads1)[1024] = readclass->read1->reads;
	int i,j;
	for (i=0;i<4;i++){
		j=0;
		while( (*cur1++ = reads1[i][j++]) != '\n' ){}
	}
	*cur1 = '\0';
	job->l1 = strlen(job->buffer1);
	if (g.single_or_pair){
		if (job->buffer2 == NULL) job->buffer2 = malloc(1024*4);
		char* cur2 = job->buffer2;
		char (*reads2)[1024] = readclass->read2->reads;
		for (i=0;i<4;i++){
			j=0;
			while( (*cur2++ = reads2[i][j++]) != '\n' ){}
		}
		*cur2 = '\0';
		job->l2 = strlen(job->buffer2);
	}
		
}


local char checkadapter(struct readclass* readclass){
	char a,b;
	b = 1;
	readline* read1 = readclass->read1;
	if (read1->adapter){
		read1->adapter->n = read1->seqlen;
		read1->adapter->s2 = (read1->reads)[1];
		a = match(read1->adapter);
	}else{
		a= get_index(g.table1, read1->id);
	}

	if (g.single_or_pair){
		readline* read2 = readclass->read2;
		if (read2->adapter){
			read2->adapter->n = read2->seqlen;
			read2->adapter->s2 = (read2->reads)[1];
			b = match(read2->adapter);
		}else{
			b= get_index(g.table2, read2->id);
		}
	}
	return(a * b);
}

local void quality_base_percent(struct readclass* readclass){
	readline* read1 = readclass->read1;
    int i=0;
    int curqual;
    read1->highAT_sign=0;
    read1->low = 0;
	read1->N = 0;
    read1->Q30=0;
    read1->A=0;
    read1->T=0;
    read1->GC=0;
	char* seq=(read1->reads)[1];
	char* qual=(read1->reads)[3];
    while(i < g.CUT){
        switch (*(seq + i)){
            case 'A': read1->A++; break;
            case 'a': read1->A++; break;
            case 'T': read1->T++; break;
            case 't': read1->T++; break;
            case 'G': read1->GC++; break;
            case 'g': read1->GC++; break;
            case 'C': read1->GC++; break;
            case 'c': read1->GC++; break;
            case 'N': read1->N++;  break;
            case 'n': read1->N++;  break;
        }
        curqual = (int)*(qual + i) - g.ASCII;
        i++;
        if (curqual <= g.LOWQ){
           	read1->low++;
            continue;
        }

        if (curqual >= g.HIGHQ){
            read1->Q30++;
            continue;
        }
    }
    if (read1->A >= g.highAT || read1->T >= g.highAT){
    	read1->highAT_sign=1;
    }
    if (g.single_or_pair){
    	readline* read2 = readclass->read2;
    	i=0;
	    read2->highAT_sign=0;
	    read2->low = 0;
		read2->N = 0;
	    read2->Q30=0;
	    read2->A=0;
	    read2->T=0;
	    read2->GC=0;
	   	seq=(read2->reads)[1];
		qual=(read2->reads)[3];
	    while(i < g.CUT){
	        switch (*(seq + i)){
	            case 'A': read2->A++; break;
	            case 'a': read2->A++; break;
	            case 'T': read2->T++; break;
	            case 't': read2->T++; break;
	            case 'G': read2->GC++; break;
	            case 'g': read2->GC++; break;
	            case 'C': read2->GC++; break;
	            case 'c': read2->GC++; break;
	            case 'N': read2->N++;  break;
	            case 'n': read2->N++;  break;
	        }
	        curqual = (int)*(qual + i) - g.ASCII;
	        i++;
	        if (curqual <= g.LOWQ){
	            read2->low++;
	            continue;
	        }
	        if (curqual >= g.HIGHQ){
	            read2->Q30++;
	            continue;
	        }
	    }
		if (read2->A >= g.highAT || read2->T >= g.highAT){
		   	read2->highAT_sign=1;
		}
    }
    
}

local char check_ployAT(char* seq, int ployAT){
	if (ployAT>0){
		int seqlen = strlen(seq);
		int i;
		int count;
		if (seq[0] == 'A' || seq[0] == 'T') count = 1;
		for(i=1;i<seqlen;i++){
			if (seq[i]=='A' || seq[i]=='T'){
				if (seq[i] == seq[i-1]){
					count++;
					if (count == ployAT) return(1);
					continue;
				}
			}
			count = 0;
		}
	}
	return(0);
}


local void report(struct readclass* readclass){
	readline* read1 = readclass->read1;
	STAT* stat = readclass->stat;
    int i=0;
    int curqual;
    char* seq=(read1->reads)[1];
	char* qual=(read1->reads)[3];
    while (i < g.CUT){
        switch (*(seq + i)){
            case 'A': (*(*(stat->base1+i)+0))++; break;
            case 'a': (*(*(stat->base1+i)+0))++; break;
            case 'T': (*(*(stat->base1+i)+1))++; break;
            case 't': (*(*(stat->base1+i)+1))++; break;
            case 'C': (*(*(stat->base1+i)+2))++; break;
            case 'c': (*(*(stat->base1+i)+2))++; break;
            case 'G': (*(*(stat->base1+i)+3))++; break;
            case 'g': (*(*(stat->base1+i)+3))++; break;
            case 'N': (*(*(stat->base1+i)+4))++; break;
            case 'n': (*(*(stat->base1+i)+4))++; break;
        }
        curqual=(int)*(qual + i) - g.ASCII;
        (*(*(stat->quality1 + i) + curqual))++;
        i++;
    }

    if (g.single_or_pair){
    	readline* read2 = readclass->read2;
    	i = 0;
	    seq=(read2->reads)[1];
	    qual=(read2->reads)[3];
	    while (i < g.CUT){
			switch (*(seq + i)){
				case 'A': (*(*(stat->base2+i)+0))++; break;
				case 'a': (*(*(stat->base2+i)+0))++; break;
				case 'T': (*(*(stat->base2+i)+1))++; break;
				case 't': (*(*(stat->base2+i)+1))++; break;
				case 'C': (*(*(stat->base2+i)+2))++; break;
				case 'c': (*(*(stat->base2+i)+2))++; break;
				case 'G': (*(*(stat->base2+i)+3))++; break;
				case 'g': (*(*(stat->base2+i)+3))++; break;
				case 'N': (*(*(stat->base2+i)+4))++; break;
				case 'n': (*(*(stat->base2+i)+4))++; break;
			}
			curqual=(int)*(qual + i) - g.ASCII;
			(*(*(stat->quality2 + i) + curqual))++;
			i++;
		}
    }
	    
}
//--------------------------------------------------------------


local void prepare(void){
	char line[1024];
	g.gzfp1 = NULL;
	g.gzfp2 = NULL;
	g.gzfp1 = gzopen(g.fq1,"rb");
	if (g.single_or_pair) g.gzfp2 = gzopen(g.fq2,"rb");
	g.fpout1 = NULL;
	g.fpout2 = NULL;
	if ((g.fpout1 = fopen(g.outfq1,"w")) == NULL){
		fprintf(stderr, "%s can not open\n", g.outfq1);
		exit(0);
	}
	if (g.single_or_pair) {
		if ((g.fpout2 = fopen(g.outfq2,"w")) == NULL){
			fprintf(stderr, "%s can not open\n", g.outfq2);
			exit(0);
		}	
	}

	
	if (gzgets(g.gzfp1,line,1024)==NULL){
		exit(0);
	}
	if (gzgets(g.gzfp1,line,1024)==NULL){
		exit(0);
	}
	int readlen = strlen(line)-1-g.start;

	if (g.CUT == 0){
		g.CUT = readlen;
	}else{
		if (readlen < g.CUT){
			g.CUT = readlen;
		}
	}
	if (g.start + g.CUT > strlen(line)){
		fprintf(stderr, "start: %d  CUT: %d range out seqlen: %d\n", g.start, g.CUT, (int)strlen(line)-1);
		exit(0);
	}
	if (g.QUOTA == 0){
		g.QUOTA = LONG_MAX;
		g.outraw = 0;
	}else{
		if (g.single_or_pair==0) g.QUOTA /= g.CUT;
		if (g.single_or_pair==1) g.QUOTA /= g.CUT*2;
		g.outraw = 1;
	}

	g.fpout1_raw = NULL;
	g.fpout2_raw = NULL;
	if (g.outraw){
		if ((g.fpout1_raw = fopen(g.rawcutfq1,"w"))==NULL){
			fprintf(stderr, "%s can not open\n", g.rawcutfq1);
			exit(0);	
		}
		if (g.single_or_pair){
			if ((g.fpout2_raw = fopen(g.rawcutfq2,"w"))==NULL){

				fprintf(stderr, "%s can not open\n", g.rawcutfq2);
				exit(0);	
			}	
		}
	}


	g.LOW_QUALITY = g.LOW_QUALITY * g.CUT;
	g.highN = g.highN * g.CUT;
	g.highAT = g.highAT * g.CUT;
	g.highGC = g.highGC * g.CUT;
	gzseek(g.gzfp1,0L,SEEK_SET);
}

local void process(void){
	prepare();
	if (g.Nthread > 1){
		parallel_filter_main();
	}else{
		single_filter_main();
	}

	//判断read1 和 read2是否都读完？------------------------------
	if (g.single_or_pair){
		if (g.read1Count != g.read2Count){
			fprintf(stderr, "read1 num not equal read2 num\n");
			exit(-1);
		}
	}


	//gzip--------------------------------------------------------
	char s1[strlen(g.outfq1)+100];
	char s2[strlen(g.outfq1)+100];
	char s3[strlen(g.rawcutfq1)+100];
	char s4[strlen(g.rawcutfq2)+100];
	sprintf(s1,"%s %s",g.speedy_gzip,g.outfq1);
	if(g.single_or_pair) sprintf(s2,"%s %s",g.speedy_gzip,g.outfq2);
	if (g.QUOTA > 0 && g.QUOTA < LONG_MAX){
		sprintf(s3,"%s %s",g.speedy_gzip,g.rawcutfq1);
		if(g.single_or_pair) sprintf(s4,"%s %s",g.speedy_gzip,g.rawcutfq2);
	}
	char* shell1[]={"sh","-c",s1,NULL};
	char* shell2[]={"sh","-c",s2,NULL};
	char* shell3[]={"sh","-c",s3,NULL};
	char* shell4[]={"sh","-c",s4,NULL};

	pid_t jobs[4]={0,0,0,0};

	int n=0;
	pid_t curjob;
	curjob = fork();
	if(curjob==0){
		execv("/bin/sh",shell1);
	}
	jobs[n] = curjob;
	n++;

	if (g.single_or_pair){
		curjob = fork();
		if(curjob==0){
			execv("/bin/sh",shell2);	
		}
		jobs[n] = curjob;
		n++;
	}

	if (g.QUOTA > 0 && g.QUOTA < LONG_MAX){
		curjob = fork();
		if(curjob==0){
			execv("/bin/sh",shell3);
		}
		jobs[n] = curjob;
		n++;
		if (g.single_or_pair){
			curjob = fork();
			if(curjob==0){
				execv("/bin/sh",shell4);
			}
			jobs[n] = curjob;
			n++;
		}
	}
		
	n--;

	while(n > -1){
		curjob = waitpid(jobs[n],NULL,0);
		if (curjob == jobs[n]){
			n--;
		}else{
			fprintf(stderr, "%d is wrong! return %d\n",jobs[0],curjob);
			exit(-1);
		}
	} 

}


local void readadapter(void){
	gzFile* gzfp=gzopen(g.adapter1,"r");
	if (gzfp==NULL){
		perror("adapter1 ");
		exit(0);
	}
	char line[1024];
	char* cur;
	while (gzgets(gzfp,line,1024)!=NULL){
		if (strlen(line)==1 && strcmp(line,"\n")==0){
			continue;
		}
		cur = (char*)memchr(line,'\t',1024);
		if (cur){
			*cur = '\n';
			cur++;
			*cur = '\0';
			insert_keyvalue(&g.table1,line,1);
		}
	}
	gzclose(gzfp);
	if (g.single_or_pair){
		gzfp=gzopen(g.adapter2,"r");
		if (gzfp==NULL){
			perror("adapter2 ");
			exit(0);
		}
		while (gzgets(gzfp,line,1024)!=NULL){
			if (strlen(line)==1 && strcmp(line,"\n")==0){
				continue;
			}
			cur = (char*)memchr(line,'\t',1024);
			if (cur!=NULL){
				*cur = '\n';
				cur++;
				*cur = '\0';
				insert_keyvalue(&g.table2,line,1);
			}
		}
		gzclose(gzfp);
	}
}

local void outstat(const char* statfile){
	FILE* out=fopen(statfile,"w");
	//allstat.raw_bases = (long long)allstat.raw_reads * g.CUT * 4;
	int xx=1;
	if(g.single_or_pair) xx=2;
	allstat.clean_bases = (long long)allstat.clean_reads * g.CUT * xx;

	allstat.clean_reads_rate = 100.0 * allstat.clean_reads / allstat.raw_reads;
	
	allstat.low_quality_rate = 100.0 * allstat.low_quality_reads / allstat.raw_reads ;

	allstat.Ns_rate = 100.0 * allstat.ns_reads / allstat.raw_reads;

	allstat.adapter_rate = 100.0 * allstat.adapter_reads / allstat.raw_reads;

	allstat.raw_Q30_bases_rate = 100.0 * allstat.raw_Q30_bases / allstat.raw_bases;

	allstat.clean_Q30_bases_rate = 100.0 * allstat.clean_Q30_bases / allstat.clean_bases;

	allstat.ployAT_rate = 100.0 * allstat.ployAT_reads / allstat.raw_reads;

	allstat.HighAT_rate = 100.0 * allstat.HighAT_reads / allstat.raw_reads;

	allstat.HighGC_rate = 100.0 * allstat.HighGC_reads / allstat.raw_reads;

	fprintf(out,"Read Length(bp)\t%d\n",g.CUT);
	fprintf(out,"Original reads number\t%u\n",allstat.raw_reads*xx);
	fprintf(out,"Original bases number\t%lld\n",allstat.raw_bases);
	fprintf(out,"Clean reads number\t%u\n",allstat.clean_reads*xx);
	fprintf(out,"Clean bases number\t%lld\n",allstat.clean_bases);
	fprintf(out,"Clean reads rate(%%)\t%.3f\n",allstat.clean_reads_rate);
	fprintf(out,"Low-quality reads number\t%u\n",allstat.low_quality_reads*xx);
	fprintf(out,"Low-quality reads rate(%%)\t%.3f\n",allstat.low_quality_rate);
	fprintf(out,"Ns reads number\t%u\n",allstat.ns_reads*xx);
	fprintf(out,"Ns reads rate(%%)\t%.3f\n",allstat.Ns_rate);
	fprintf(out,"Adapter polluted reads number\t%u\n",allstat.adapter_reads*xx);
	fprintf(out,"Adapter polluted reads rate(%%)\t%.3f\n",allstat.adapter_rate);
	fprintf(out,"HighAT Reads(normal)\t%u\n",allstat.HighAT_reads*xx);
	fprintf(out,"HighAT Reads(normal) Rate(%%)\t%.3f\n",allstat.HighAT_rate);
	fprintf(out,"HighGC Reads\t%u\n",allstat.HighGC_reads*xx);
	fprintf(out,"HighGC Reads Rate(%%)\t%.3f\n",allstat.HighGC_rate);
	fprintf(out,"PloyAT Reads\t%u\n",allstat.ployAT_reads*xx);
	fprintf(out,"PloyAT Reads Rate(%%)\t%.3f\n",allstat.ployAT_rate);
	fprintf(out,"Original Q30 bases rate(%%)\t%.3f\n",allstat.raw_Q30_bases_rate);
	fprintf(out,"Clean Q30 bases rate(%%)\t%.3f\n",allstat.clean_Q30_bases_rate);
	
	fclose(out);
}

local void outreport(unsigned int (*base)[5], unsigned int (*quality)[50],  const char* outfile){
	FILE* fp=fopen(outfile,"w");
	fprintf(fp,"#\tA\tT\tC\tG\tN");
	int m,n;
	long long total_base = (long long)allstat.clean_reads * g.CUT;
	unsigned int q30=0;
	unsigned int q20=0;
	unsigned int N=0;
	unsigned int GC=0;

	float GC_bias=0;
	float AT_bias=0;
	float cur_GC_bias=0;
	float cur_AT_bias=0;

	for (m=0;m<42;m++){
		fprintf(fp,"\t%d",m);
	}
	fprintf(fp,"\n");
	for (n=0;n<g.CUT;n++){
		GC += base[n][2] + base[n][3];
		N += base[n][4];
		fprintf(fp,"%d\t%u\t%u\t%u\t%u\t%u",n+1,base[n][0],base[n][1],base[n][2],base[n][3],base[n][4]);
		cur_AT_bias=1.0 * abs(base[n][0] - base[n][1])/allstat.clean_reads;
		cur_GC_bias=1.0 * abs(base[n][2] - base[n][3])/allstat.clean_reads;
		if (cur_AT_bias > AT_bias){
			AT_bias = cur_AT_bias;
		}
		if (cur_GC_bias > GC_bias){
			GC_bias = cur_GC_bias;
		}
		for (m=0;m<42;m++){
			fprintf(fp,"\t%d",quality[n][m]);
			if (m >= 20){
				q20 += quality[n][m];
			}
			if (m >= 30){
				q30 += quality[n][m];
			}
		}	
		fprintf(fp,"\n");
	}
	fprintf(fp, "#total_reads\t%d\n",allstat.clean_reads);
	fprintf(fp, "#total_base\t%lld\n", total_base);
	fprintf(fp, "#GC_bias(%%)\t%.3f\n",(float)100.0*GC_bias);
	fprintf(fp, "#AT_bias(%%)\t%.3f\n",(float)100.0*AT_bias);
	fprintf(fp, "#N_percent(%%)\t%.3f\n",(float)100.0*(long long)N/total_base);
	fprintf(fp, "#Q20(%%)\t%.3f\n",(float)100.0*(long long)q20/total_base);
	fprintf(fp, "#Q30(%%)\t%.3f\n",(float)100.0*(long long)q30/total_base);
	fprintf(fp, "#GC_percent(%%)\t%.3f\n",(float)100.0*(long long)GC/total_base);
	fclose(fp);
}

local int parameter(const char* str, void* p, const char* type){
	int len=strlen(str);
	int i=0;

	if (strcmp(type,"int")==0){
		while(i<len){
			if(isdigit(*(str+i))==0){
				return(0);
			}
			i++;
		}
		int* x=(int*)p;
		*x=(int)atoi(str);
		return(1);
	}else if (strcmp(type,"float")==0){
		int dot=0;
		while(i<len){
			if(*(str+i) == '.' && dot == 0){
				i++;
				dot++;
				continue;
			}
			if(isdigit(*(str+i))==0){
				return(0);
			}
			i++;
		}
		float* x=(float*)p;
		*x=(float)atof(str);
		return(1);
	}else if (strcmp(type,"str")==0){
		char* x=*(char**)p;
		strcpy(x,str);
		return(1);
	}else if (strcmp(type,"long")==0){
		while(i<len){
			if(isdigit(*(str+i))==0){
				return(0);
			}
			i++;
		}
		long* x=(long*)p;
		*x=atol(str);
		return(1);
	}else{
		return(0);
	}
}

int main(int argc, char** argv){

	if (argc < 2){
		puts(HELP);
		exit(-1);
	}
	int function=0;
	char* HELP_output=NULL;
	if (strcmp(argv[1],"filter")==0){
		function = 1;
		HELP_output = (char*)HELP1;
	}else if(strcmp(argv[1],"cutadapt")==0){
		function = 2;
		HELP_output = (char*)HELP2;
	}else{
		fprintf(stderr, "%s\n",HELP);
		exit(0);
	}

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);
	time_t start,end;
	time(&start);
	//fprintf(stdout,"%s\n", ctime(&start));
	//-------------------------------------------------
	g.single_or_pair=0;
	g.CUT=0;
	g.start=0;
	g.QUOTA=0;
	g.LOWQ = 19;
	g.HIGHQ = 30;
	g.LOW_QUALITY=0.5;
	g.highN=0.05;
	g.highAT=0.8;
	g.highGC=0.8;
	g.ployAT=0;
	g.ASCII = 33;
	g.Nthread = 1;

	g.read1Count=0;
	g.read2Count=0;

	g.sequence1 = (char*)malloc(512);
	g.sequence2 = (char*)malloc(512);

	g.min_overlap = 5;
	g.table1=init_hash();
	g.table2=init_hash();
	g.fq1 = (char*)malloc(1024);
	g.fq2 = (char*)malloc(1024);
	g.rawcutfq1 = (char*)malloc(1024);
	g.rawcutfq2 = (char*)malloc(1024);
	g.adapter1 = (char*)malloc(1024);
	g.adapter2 = (char*)malloc(1024);
	g.outfq1 = (char*)malloc(1024);
	g.outfq2 = (char*)malloc(1024);
	g.speedy_gzip = (char*)malloc(1024);
	char* files1[3]={g.fq1, g.adapter1, g.outfq1};
	char* files2[6]={g.fq1, g.fq2, g.adapter1, g.adapter2, g.outfq1, g.outfq2};
	char* files3[2]={g.fq1, g.outfq1};
	char* files4[4]={g.fq1, g.fq2, g.outfq1, g.outfq2};
	char** selectfiles=NULL;
	int n=0;
	int i;
	char x=1;
	for (i=2; i < argc; i++){
		if (argv[i][0] == '-'){
			if (strcmp(argv[i],"--c")==0 || strcmp(argv[i],"-c")==0){
				if(parameter(argv[i+1],&g.CUT,"int")==0){
					puts("--c INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--start")==0 || strcmp(argv[i],"-start")==0){
				if(parameter(argv[i+1],&g.start,"int")==0){
					puts("--start INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if(strcmp(argv[i],"--q")==0 || strcmp(argv[i],"-q")==0){
				if(parameter(argv[i+1],&g.QUOTA,"long")==0){
					puts("--q INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--ql")==0 || strcmp(argv[i],"-ql")==0){
				if(parameter(argv[i+1],&g.LOWQ,"int")==0){
					puts("--ql INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--qh")==0 || strcmp(argv[i],"-qh")==0){
				if(parameter(argv[i+1],&g.HIGHQ,"int")==0){
					puts("--qh INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--Q")==0 || strcmp(argv[i],"-Q")==0){
				if(parameter(argv[i+1],&g.LOW_QUALITY,"float")==0){
					puts("--Q FLOAT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--n")==0 || strcmp(argv[i],"-n")==0){
				if(parameter(argv[i+1],&g.highN,"float")==0){
					puts("--n FLOAT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--AT")==0 || strcmp(argv[i],"-AT")==0){
				if(parameter(argv[i+1],&g.highAT,"float")==0){
					puts("--highAT FLOAT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--GC")==0 || strcmp(argv[i],"-GC")==0){
				if(parameter(argv[i+1],&g.highGC,"float")==0){
					puts("--highGC FLOAT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--ployAT")==0 || strcmp(argv[i],"-ployAT")==0){
				if(parameter(argv[i+1],&g.ployAT,"int")==0){
					puts("--ployAT INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--i")==0 || strcmp(argv[i],"-i")==0){
				if(parameter(argv[i+1],&g.ASCII,"int")==0){
					puts("--i INT\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--adapt1")==0 || strcmp(argv[i],"-adapt")==0){
				if(parameter(argv[i+1],&g.sequence1,"str")==0){
					puts("--adapt1 string\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--adapt2")==0 || strcmp(argv[i],"-adapt")==0){
				if(parameter(argv[i+1],&g.sequence2,"str")==0){
					puts("--adapt2 string\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"-min_overlap")==0 || strcmp(argv[i],"--min_overlap")==0){
				if(parameter(argv[i+1],&g.min_overlap,"int")==0){
					puts("-min_overlap int\n");
					puts(HELP_output);
					exit(-1);
				}
				i++;
				continue;
			}else if (strcmp(argv[i],"--m")==0 || strcmp(argv[i],"-m")==0){
				if(parameter(argv[i+1],&g.Nthread,"int")==0){
					puts("--m int\n");
					puts(HELP);
					exit(-1);
				}
				i++;
				continue;
			}else{
				fprintf(stderr,"invalid option %s\n",argv[i]);
				exit(0);
			}
		}else {
			if (function == 2){
				strcpy(g.fq1,argv[i]);
				break;	
			}
			if (x){
				if (argc - i == 2){
					g.single_or_pair = 0;
					if (*g.sequence1) selectfiles = files3;
					else selectfiles = files1;
				}else if(argc - i == 3){
					g.single_or_pair = 0;
					if (*g.sequence1) selectfiles = files3;
					else selectfiles = files1;
				}else if(argc - i == 4){
					g.single_or_pair = 1;
					if (*g.sequence1) selectfiles = files4;
					else selectfiles = files2;
				}else if(argc - i == 6){
					g.single_or_pair = 1;
					if (*g.sequence1) selectfiles = files4;
					else selectfiles = files2;
				}else{
					fprintf(stderr, "wrong !\n%s\n",HELP_output);
				}
			}
			x=0;
			strcpy(selectfiles[n],argv[i]);
			n++;
		}
	}

	if (*g.fq1 == '\0') {printf("%s\n",HELP_output);exit(0);}
	if (g.single_or_pair && *g.fq2 == '\0') {printf("%s\n",HELP_output);exit(0);}
	if (*g.outfq1 == '\0' && function == 1) {printf("%s\n",HELP_output);exit(0);}
	if (g.single_or_pair && *g.outfq2 == '\0') {printf("%s\n",HELP_output);exit(0);}

	int out1len = strlen(g.outfq1);
	int out2len = strlen(g.outfq2);
	if (strcmp(g.outfq1+out1len-3,".gz")==0) *(g.outfq1+out1len-3) = '\0';
	if (g.single_or_pair && strcmp(g.outfq2+out2len-3,".gz")==0) *(g.outfq2+out2len-3) = '\0';

	if (*g.sequence1 == '\0'){
		if (function==1) {
			if (*g.adapter1 == '\0') {printf("%s\n",HELP_output);exit(0);}
			if (g.single_or_pair && *g.adapter2 == '\0') {printf("%s\n",HELP_output);exit(0);}
		}else if (function==2){
			printf("%s\n",HELP_output);exit(0);
		}
	}else{
		if (g.single_or_pair){
			if (*g.sequence2 == '\0'){fprintf(stderr, "%s\n",HELP_output);exit(0);}
		}
	}

	if (g.QUOTA < 0){
		puts("-q num > 0");
		exit(-1);
	}else if(g.QUOTA > 0){
		sprintf(g.rawcutfq1,"%s.cut.raw",g.outfq1);
		if (g.single_or_pair) sprintf(g.rawcutfq2,"%s.cut.raw",g.outfq2);
	}
	
	if (g.CUT < 0){
		puts("--c num >= 0");
		exit(-1);
	}

	if (g.Nthread < 3){
		g.Nthread = 1;
	}else{
		if (g.QUOTA != 0){
			if (g.Nthread < 4){
				fprintf(stderr, "--m 不正确\n%s\n",HELP1);
				exit(0);
			}
		}
	}

	//-----------------------------------------
	if (function ==1){
		char statfile[strlen(g.outfq1)+20];
		sprintf(statfile,"%s.gz.stat",g.outfq1);
		char reportfile1[strlen(g.outfq1)+20];
		char reportfile2[strlen(g.outfq2)+20];
		sprintf(reportfile1,"%s.gz.report",g.outfq1);
		sprintf(reportfile2,"%s.gz.report",g.outfq2);
		if (*g.sequence1 == '\0'){
			readadapter();
		}
		if (g.Nthread >= 3){
			sprintf(g.speedy_gzip,"%s -f -p 3 ",speedy_gzip_software);
		}else{
			sprintf(g.speedy_gzip,"%s -f -p 1 ",speedy_gzip_software);	
		}
		//-----------------main func---------------
		process();
		//---------------outstat-------------------
		outstat(statfile);
		outreport(allstat.base1, allstat.quality1, reportfile1);
		if (g.single_or_pair) outreport(allstat.base2, allstat.quality2, reportfile2);
		//-----------------------------------------
	}else if (function == 2){
		cutadapt_main(g.fq1,g.sequence1,g.min_overlap);
	}
		
	time(&end);
	//fprintf(stdout,"%s\n", ctime(&end));
	return(0);
}

