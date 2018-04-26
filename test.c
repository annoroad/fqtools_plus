#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char check_ployAT(char* seq, int ployAT){
	if (ployAT>0){
		int seqlen = strlen(seq);
		int i,j;
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


int main(){
	char* string ="GTCAGCTAAAAACTGTGTTTTTTTTGTTTTAAAAA";
	if (check_ployAT(string,9)){
		puts("yes");
	}

	return(1);

}