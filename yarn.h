#ifndef __YARN
#define __YARN



typedef struct lock_s lock;
typedef struct thread_s thread;

enum twist_op { TO, BY };
enum wait_for {TO_BE,NOT_TO_BE,TO_MORE_THAN,TO_LESS_THAN};


lock* new_lock(long inital);
void free_lock(lock* bolt);
void possess(lock* bolt);
void release(lock* bolt);
thread* launch(void (*probe)(void*), void* args);
void twist(lock* bolt, enum twist_op op, long val);
void wait_for(lock* bolt, enum wait_for op, long val);
void join(thread* ally);
void join_all(void);
long peek_lock(lock *bolt);
void destruct(thread* off_course);


#endif