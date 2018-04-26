#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

//lixiangma

enum FLAG {
	START_WITHIN_SEQ1 = 1,
	START_WITHIN_SEQ2 = 2,
	STOP_WITHIN_SEQ1 = 4,
	STOP_WITHIN_SEQ2 = 8,
	SEMIGLOBAL = 15,
	ALLOW_WILDCARD_SEQ1 = 1,
	ALLOW_WILDCARD_SEQ2 = 2,
};

#define GAPCHAR '\0'

#define DELETION_COST 1
#define INSERTION_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0
#define degenerate 1
#define flags 14
#define error_rate 0.1
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#include "cutadapt.h"


static void alignread(adapter_* adapter);
static void printmatch(adaptline_* adaptline);
static void func(adaptline_* adaptline);
static char readfq(adaptline_* adaptline);


adapter_* initalign(const char* s1, const int min_overlap){
	adapter_* adapter = (adapter_*)malloc(sizeof(adapter_));
	adapter->min_overlap = min_overlap;
	strcpy(adapter->s1,s1);
	adapter->m = strlen(s1);
	adapter->s2 = NULL;
	adapter->columns = (Entry_*)malloc(sizeof(Entry_)*(adapter->m+1));
	adapter->result = (result_*)malloc(sizeof(result_));
	return(adapter);
}


static void alignread(adapter_* adapter)
{

	int m = adapter->m;
	int n = adapter->n;
	Entry_* column = adapter->columns;

	int i, j, best_i, best_j, best_cost, best_matches, best_origin;

	for (i = 0; i <= m; i++) {
		column[i].matches = 0;
		column[i].cost = i * DELETION_COST;
		column[i].origin = 0;
	}

	best_i = m;
	best_j = 0;
	best_cost = column[m].cost;
	best_matches = 0;
	best_origin = column[m].origin;

	// maximum no. of errors
	int k = error_rate * m;
	int last = k + 1;
	int match;
	int cost_diag;
	int cost_deletion;
	int cost_insertion;
	int origin, cost, matches;
	for (j = 1; j <= n; j++) {
		Entry_ tmp_entry = column[0];
		column[0].cost = 0;
		column[0].origin = j;
		column[0].matches = 0;
		for (i = 1; i <= last; ++i) {

			match = (((adapter->s1)[i-1] == (adapter->s2)[j-1]) || ((adapter->s1)[i-1] == 'N'));
			cost_diag = tmp_entry.cost + (match?MATCH_COST : MISMATCH_COST);
			cost_deletion = column[i].cost + DELETION_COST;
			cost_insertion = column[i-1].cost + INSERTION_COST;

			if (cost_diag <= cost_deletion && cost_diag <= cost_insertion) {
				// MATCH or MISMATCH
				cost = cost_diag;
				origin = tmp_entry.origin;
				matches = tmp_entry.matches + match;
			} else if (cost_insertion <= cost_deletion) {
				// INSERTION
				cost = cost_insertion;
				origin = column[i-1].origin;
				matches = column[i-1].matches;
			} else {
				// DELETION
				cost = cost_deletion;
				origin = column[i].origin;
				matches = column[i].matches;
			}
			//memcpy(&tmp_entry,column[i],l);
			tmp_entry = column[i];
			// remember current cell for next iteration
			column[i].cost = cost;
			column[i].origin = origin;
			column[i].matches = matches;

		}

		while (column[last].cost > k) {
			last--;
		}
		if (last < m) {
			last++;
		} else {
		
			int length = m + min(column[m].origin, 0);
			cost = column[m].cost;
			matches = column[m].matches;
			if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
				// update
				best_matches = matches;
				best_cost = cost;
				best_origin = column[m].origin;
				best_i = m;
				best_j = j;
			}

		}
		// column finished
	}


	for (i = 0; i <= m; ++i) {
		int length = i + min(column[i].origin, 0);
		int cost = column[i].cost;
		int matches = column[i].matches;
		if (cost <= length * error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
			best_matches = matches;
			best_cost = cost;
			best_origin = column[i].origin;
			best_i = i;
			best_j = n;
		}
	}


	if (best_origin >= 0) {
		adapter->result->start1 = 0;
		adapter->result->start2 = best_origin;
	} else {
		adapter->result->start1 = -best_origin;
		adapter->result->start2 = 0;
	}
	adapter->result->stop1 = best_i;
	adapter->result->stop2 = best_j;
	adapter->result->matches = best_matches;
	adapter->result->errors = best_cost;
	// return (start1, stop1, start2, stop2, matches, errors)
}


char match(adapter_* adapter){
	alignread(adapter);
	if (adapter->result->stop1 - adapter->result->start1 < adapter->min_overlap){
		return(1);
	}
	//printf("%d\t%d\t%d\t%d\t%d\t%d\n",adapter->result->start1,adapter->result->stop1,adapter->result->start2,adapter->result->stop2,adapter->result->matches,adapter->result->errors);	
	return(0);
}

static void printmatch(adaptline_* adaptline){
	char* seq = adaptline->adapter->s2;
	int len = adaptline->adapter->n;
	result_* align = adaptline->adapter->result;
	char seq1[1024];
	char seq2[1024];
	char seq3[1024];
	memset(seq1,0,1024);
	memset(seq2,0,1024);
	memset(seq3,0,1024);
	memcpy(seq1,seq,align->start2);
	memcpy(seq2,seq+align->start2,align->stop2-align->start2);
	memcpy(seq3,seq+align->stop2, len - align->stop2);
	fprintf(stdout,"%s\t%d\t%d\t%d\t%s\t%s\t%s\n",adaptline->id,align->errors,align->start2,align->stop2,seq1,seq2,seq3);
		
}

static char readfq(adaptline_* adaptline){
	char line[1024];
	char* cur = NULL;
	char* cur1 = NULL;
	if (gzgets(adaptline->gzfp,line,1024)!=NULL){
		cur = line+1;
		cur1 = adaptline->id;
		while(*cur!='\n'){
			*cur1++ = *cur++;
		}
		*cur1 = '\0';
	}else return(1);
	if (gzgets(adaptline->gzfp,line,1024)!=NULL){
		cur = line;
		cur1 = adaptline->adapter->s2;
		int i=0;
		while(*cur!='\n'){
			*cur1++ = *cur++;
			i++;
		}
		*cur1 = '\0';
		adaptline->adapter->n = i;
	}else return(1);
	if (gzgets(adaptline->gzfp,line,1024)!=NULL){}
	else return(1);
	if (gzgets(adaptline->gzfp,line,1024)!=NULL){}
	else return(1);
	return(0);
}


static void func(adaptline_* adaptline){
	while(1){
		if (readfq(adaptline))break;
		if (match(adaptline->adapter)==0){
			printmatch(adaptline);
		}
	}
}


void cutadapt_main(char* fq1, char* adapter_seq, int min_overlap){
	gzFile* gzfp1 = gzopen(fq1,"rb");
	if (gzfp1 == NULL){
		fprintf(stderr, "%s not open\n",fq1 );
		exit(1);
	}
	adaptline_* adaptline1 = (adaptline_*)malloc(sizeof(adaptline_));
	adaptline1->adapter = initalign(adapter_seq,min_overlap);
	adaptline1->adapter->s2 = (char*)malloc(sizeof(char)*1024);
	adaptline1->id = (char*)malloc(sizeof(char)*1024);
	adaptline1->gzfp = gzfp1;
	func(adaptline1);
}

/*
int main(){
	char* fq1="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/Hic/ngs_bioinfo/hic-13/malixiang/reseach/fqtools/E-1_R1.fq.gz";
	char* adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
	int min_overlap= 5 ;
	cutadapt_main(fq1,adapter,min_overlap);
	
	return(1);
}
*/