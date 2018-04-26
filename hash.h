#ifndef MY_HASH
#define MY_HASH
#define HT_SIZE 1048576
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

typedef struct keyvalue{
	char* key;
	int value;
	struct keyvalue* next;
} keyvalue;


typedef struct hashtable{
	int htsize;
	int* bucketlengths;
	keyvalue** blist;
}hashtable;

unsigned int hash_string(const char* str,int htsize);
hashtable* init_hash(void);
void insert_keyvalue(hashtable** ht, char* key, int value);
char get_index(hashtable* ht, char* key);

#endif