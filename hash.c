#include "hash.h"

unsigned int hash_string(const char* key, int htsize){
	unsigned long hash = 5381;
	while(*key){
		hash += (hash<<5) + (int)(*key++);
		if (hash >= htsize){
			hash = hash % htsize;
		}
	}
	return(hash);
}

hashtable* init_hash(void){
	hashtable* ht;
	ht=(hashtable*)malloc(sizeof(hashtable));
	ht->htsize=HT_SIZE;
	ht->bucketlengths=(int*)malloc(sizeof(int)*ht->htsize);
	ht->blist=(keyvalue**)malloc(sizeof(keyvalue*)*ht->htsize);
	int i;
	for(i=0;i<ht->htsize;i++){
		(ht->blist)[i]=NULL;
	}
	return(ht);
}

void insert_keyvalue(hashtable** ht, char* key, int value){
	unsigned int hash_value=hash_string(key,(*ht)->htsize);
	keyvalue* newhash=(keyvalue*)malloc(sizeof(keyvalue));
	newhash->key=(char*)malloc(sizeof(char)*(strlen(key)+1));
	strcpy(newhash->key,key);
	newhash->key[strlen(key)]='\0';
	newhash->value=value;
	newhash->next=(*ht)->blist[hash_value];
	(*ht)->blist[hash_value]=newhash;
	(*ht)->bucketlengths[hash_value]++;
}


char get_index(hashtable* ht, char* key){
	unsigned int hash_value=hash_string(key,ht->htsize);
	keyvalue* curkeyvalue=ht->blist[hash_value];
	while(curkeyvalue != NULL){
		if(strcmp(curkeyvalue->key, key)==0){
			return(0);
			break;
		}
		curkeyvalue=curkeyvalue->next;
	}
	return(1);
}
