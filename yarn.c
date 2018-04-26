#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "yarn.h"

#define fail() {fprintf(stderr,"yarn wrong!\n");exit(0);}

#define local static

struct lock_s{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    long value;
};


struct thread_s{
    pthread_t id;
    struct thread_s* next;
    int done;
};

local lock threads_lock = {
    PTHREAD_MUTEX_INITIALIZER,
    PTHREAD_COND_INITIALIZER,
    0
};

local thread* threads = NULL;


lock* new_lock(long inital){
    lock* bolt = (lock*)malloc(sizeof(lock));
    if (pthread_mutex_init(&(bolt->mutex),NULL) ||
        pthread_cond_init(&(bolt->cond),NULL)){
        fail();
    }
    bolt->value = inital;
    return(bolt);
}

void free_lock(lock* bolt){
    if (bolt){
        if (pthread_cond_destroy(&(bolt->cond)) ||
            pthread_mutex_destroy(&(bolt->mutex))){
            fail();
        }
    }
    free(bolt);
}

void possess(lock* bolt){
    if (pthread_mutex_lock(&(bolt->mutex))){
        fail();
    }
}

void release(lock* bolt){
    if (pthread_mutex_unlock(&(bolt->mutex))){
        fail();
    }   
}

struct capsule{
    void (*probe)(void*);
    void* payload;
};


local void rtn(void){
    pthread_t me;
    thread *match, **prior;

    me = pthread_self();

    possess(&(threads_lock));
    prior = &(threads);
    while((match = *prior) != NULL){
        if (pthread_equal(match->id,me)) 
            break;
        prior = &(match->next);
    }
    if (match==NULL) fail();
    match->done = 1;
    if (match != threads){
        *prior = match->next;
        match->next = threads;
        threads = match;
    }
    twist(&(threads_lock),BY,+1);
}

local void ignition(void* args){
    struct capsule* capsule = args; //!!!
    pthread_cleanup_push((void*)rtn, NULL);
    capsule->probe(capsule->payload);
    free(capsule);
    pthread_cleanup_pop(1);
}


thread* launch(void (*probe)(void*), void* args){ //!!!
    struct capsule* capsule;
    thread* th;
    pthread_attr_t attr;

    possess(&(threads_lock));
    capsule = (struct capsule*)malloc(sizeof(struct capsule));
    capsule->probe = probe;
    capsule->payload = args;
    th = (thread*)malloc(sizeof(thread));
    if (pthread_attr_init(&attr) ||
        pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE) ||
        pthread_create(&(th->id),&attr, (void*)ignition,capsule) ||
        pthread_attr_destroy(&attr)){
        fail();
    }
    th->done = 0;
    th->next = threads;
    threads = th;
    release(&(threads_lock));
    return(th);
}


void twist(lock* bolt, enum twist_op op, long val){
    if (op == TO){
        bolt->value = val;
    }else if(op == BY){
        bolt->value += val;
    }
    if (pthread_cond_broadcast(&(bolt->cond)) ||
        pthread_mutex_unlock(&(bolt->mutex))){
        fail();
    }
}

#define until(a) while(!(a))

void wait_for(lock* bolt, enum wait_for op, long val){
    switch (op){
        case TO_BE:
            until(bolt->value == val)
                if (pthread_cond_wait(&(bolt->cond),&(bolt->mutex)))
                    fail();
            break;
        case NOT_TO_BE:
            until(bolt->value != val)
                if (pthread_cond_wait(&(bolt->cond),&(bolt->mutex)))
                    fail();
            break;
        case TO_MORE_THAN:
            until(bolt->value >= val)
                if (pthread_cond_wait(&(bolt->cond),&(bolt->mutex)))
                    fail();
            break;
        case TO_LESS_THAN:
            until(bolt->value <= val)
                if (pthread_cond_wait(&(bolt->cond),&(bolt->mutex)))
                    fail();
            break;      
    }
}

void join(thread* ally){
    thread *match,**prior;
    if (pthread_join(ally->id,NULL)) fail();
    possess(&(threads_lock));
    prior = &(threads);
    while((match = *prior) != NULL){
        if (match == ally)
            break;
        prior = &(match->next);
    }
    if (match == NULL) fail();
    if (match->done)
        threads_lock.value--;
    *prior = match->next;
    release(&(threads_lock));
    free(ally);
}

void join_all(void){
    thread *match,**prior;
    possess(&(threads_lock));
    while(threads != NULL){
        wait_for(&(threads_lock),NOT_TO_BE,0);
        prior = &(threads);
        while((match = *prior) != NULL){
            if (match->done)
                break;
            prior = &(match->next);
        }
        if (pthread_join(match->id, NULL) != 0)
            fail();
        threads_lock.value--;
        *prior = match->next;
        free(match);
    }
    release(&(threads_lock));
}


long peek_lock(lock *bolt)
{
    return(bolt->value);
}

void destruct(thread* off_course){
    if (pthread_cancel(off_course->id))
        fail();
    join(off_course);
}












