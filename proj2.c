#include <semaphore.h>
#include <sys/types.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include <stdbool.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdarg.h>
#include <sys/wait.h>

// declare macro
#define MMAP(pointer) { (pointer) = mmap(NULL, sizeof(*(pointer)), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0); }
#define UNMAP(pointer) { munmap((pointer), sizeof((pointer))); }
#define mysleep(max) { usleep(1000 * (rand() % (max + 1))); }

// atoms stuff
sem_t *atMut = NULL;
char *atMut_name = "xpekni01-atMut-sem";
int *oxygen;
int *hydrogen;
sem_t *oxyQueue = NULL;
char *oxyQueue_name = "xpekni01-oxyQueue-sem";
sem_t *hydroQueue = NULL;
char *hydroQueue_name = "xpekni01-hydroQueue-sem";

// barrier stuff
sem_t *turnstile = NULL;
char *turnstile_name = "xpekni01-turnstile-sem";
sem_t *turnstile2 = NULL;
char *turnstile2_name = "xpekni01-turnstile2-sem";
sem_t *barrMut = NULL;
char *barrMut_name = "xpekni01-barrMut-sem";

// barrier2 stuff
sem_t *notTurnstile = NULL;
char *notTurnstile_name = "xpekni01-notTurnstile-sem";
sem_t *notTurnstile2 = NULL;
char *notTurnstile2_name = "xpekni01-notTurnstile2-sem";
sem_t *notBarrMut = NULL;
char *notBarrMut_name = "xpekni01-notBarrMut-sem";

// global mutex
sem_t *outputMut = NULL;
char *outputMut_name = "xpekni01-outputMut-sem";
sem_t *molMut = NULL;
char *molMut_name = "xpekni01-molMut-sem";
sem_t *sleepMut = NULL;
char *sleepMut_name = "xpekni01-sleepMut-sem";

// args
int *NO;
int *NH;
int *TI;
int *TB;

// vars
FILE *file;
int *line;
int *counter;
int *notCounter;
int *notCnt;
int *molecule;
int *notEnough;
int *notEnoughBeg;
int *atomsCreatingMol;

// functions
bool mmaping();
bool check_args(int argc, char ** argv);
bool semaphores_open();
void semaphores_destroy();
void shutdown();
void error(char *errmsg);
bool create_atoms();
void create_oxygen(int id);
void create_hydrogen(int id);
void output(const char * format, ...);
void barrier();
void notBarrier();

int main (int argc, char **argv)
{
    // checks cargs
    if (!check_args(argc, argv)) { error("ERROR: Wrong parameters."); return 1; }

    // checks for fails at the booting of the program
    if (semaphores_open()) { error("ERROR: Semaphore failed to open."); return 1; }
    
    if (mmaping()) { error("ERROR: MMAP failed."); return 1; }

    if((file = fopen("proj2.out", "w")) == NULL) { error("ERROR: File could have not been opened."); return 1; }
    

    // loads args into vars and initialize vars
    *NO = atoi(argv[1]);
    *NH = atoi(argv[2]);
    *TI = atoi(argv[3]);
    *TB = atoi(argv[4]);
    *line = 1;
    *counter = 0;
    *notCounter = 0;
    *oxygen = 0;
    *hydrogen = 0;
    *molecule = 1;
    *notEnough = false;
    *notEnoughBeg = false;
    if (*NO < 1 || *NH < 2)
    {
        *notEnoughBeg = true;
        *notCnt = *NO + *NH;
    }
    

    // starts creation of atoms and molecules
    if (!create_atoms()) return 1;
    while (wait(NULL) > 0);

    shutdown();
    return 0;
}

bool create_atoms()
{
    // hydrogen generator
    for (int i = 1; i <= *NH; i++)
    {
        pid_t id = fork();
        if (id == 0)
        {
            create_hydrogen(i);
        }
        else if (id == -1)
        {
            error("ERROR: Fork failed.");
            return false;
        }
    }
    // oxigen generator
    for (int i = 1; i <= *NO; i++)
    {
        pid_t id = fork();
        if (id == 0)
        {
            create_oxygen(i);
            return false;
        }
        else if (id == -1)
        {
            error("ERROR: Fork failed.");
        }
    }
    return true;
}

void create_oxygen(int id)
{
    // sets rand seed
    srand(time(NULL) * getpid());

    output("O %d: started", id);

    mysleep(*TI);

    output("O %d: going to queue", id);

    if (*notEnoughBeg) { notBarrier(); output("O %d: not enough H", id); exit(0); }

    // begins molecule creation
    sem_wait(atMut);
    *oxygen += 1;
    if (*hydrogen >= 2)
    {
        sem_post(hydroQueue);
        sem_post(hydroQueue);
        *hydrogen -= 2;
        sem_post(oxyQueue);
        *oxygen -= 1;
    }
    else
    {
        sem_post(atMut);
    }
    sem_wait(oxyQueue);
    
    if (*notEnough) { output("O %d: not enough H", id); exit(0); }
    
    // bond
    sem_wait(molMut);
    if (*atomsCreatingMol >= 3) { *molecule += 1; *atomsCreatingMol = 0; }   
    *atomsCreatingMol += 1; 
    sem_post(molMut);
    int mole = *molecule;
    output("O %d: creating molecule %d", id, mole);

    mysleep(*TB);
    sem_post(sleepMut);
    sem_post(sleepMut);

    barrier();
    output("O %d: molecule %d created", id, mole);

    // checks wether there is enough atoms to create other molecule
    if (*molecule >= *NO || (*molecule * 2 >= *NH) || (*NH - *molecule * 2 <= 1))
    {
        *notEnough = true;
        for (int i = 0; i < *NO; i++)
        {
            sem_post(oxyQueue);
        }
        for (int i = 0; i < *NH; i++)
        {
            sem_post(hydroQueue);
        }
    }
    
    sem_post(atMut);

    exit(0);
}

void create_hydrogen(int id)
{
    // sets rand seed
    srand(time(NULL) * getpid());

    output("H %d: started", id);

    mysleep(*TI);

    output("H %d: going to queue", id);

    if (*notEnoughBeg) { notBarrier(); output("H %d: not enough O or H", id); exit(0); }

    // begins molecule creation
    sem_wait(atMut);
    *hydrogen += 1;
    if (*hydrogen >= 2 && *oxygen >= 1)
    {
        sem_post(hydroQueue);
        sem_post(hydroQueue);
        *hydrogen -= 2;
        sem_post(oxyQueue);
        *oxygen -= 1;
    }
    else
    {
        sem_post(atMut);
    }
    sem_wait(hydroQueue);
    
    if (*notEnough) { output("H %d: not enough O or H", id); exit(0); }

    // bond
    sem_wait(molMut);
    if (*atomsCreatingMol >= 3) { *molecule += 1; *atomsCreatingMol = 0; }    
    *atomsCreatingMol += 1; 
    sem_post(molMut);
    int mole = *molecule;
    output("H %d: creating molecule %d", id, mole);

    barrier();

    sem_wait(sleepMut);
    output("H %d: molecule %d created", id, mole);
    
    exit(0);
}

void barrier()
{
    // barrier
    // rendezvous
    sem_wait(barrMut);
    *counter += 1;
    if (*counter == 3)
    {
        sem_wait(turnstile2);
        sem_post(turnstile);
    }
    sem_post(barrMut);
    
    sem_wait(turnstile);
    sem_post(turnstile);

    // critical point
    sem_wait(barrMut);
    *counter -= 1;
    if (*counter == 0)
    {
        sem_wait(turnstile);
        sem_post(turnstile2);
    }

    sem_post(barrMut);
    
    sem_wait(turnstile2);
    sem_post(turnstile2);
}

void notBarrier()
{
    // barrier
    // rendezvous
    sem_wait(notBarrMut);
    *notCounter += 1;
    if (*notCounter == *notCnt)
    {
        sem_wait(notTurnstile2);
        sem_post(notTurnstile);
    }
    sem_post(notBarrMut);
    
    sem_wait(notTurnstile);
    sem_post(notTurnstile);

    // critical point
    sem_wait(notBarrMut);
    *notCounter -= 1;
    if (*notCounter == 0)
    {
        sem_wait(notTurnstile);
        sem_post(notTurnstile2);
    }

    sem_post(notBarrMut);
    
    sem_wait(notTurnstile2);
    sem_post(notTurnstile2);
}

void output (const char * format, ...)
{
    sem_wait(outputMut);
    va_list args;
    va_start (args, format);
    fprintf (file, "%d: ", *line);
    vfprintf (file, format, args);
    fprintf (file, "\n");
    fflush(file);
    va_end (args);
    *line += 1;
    sem_post(outputMut);
}

bool mmaping()
{
    MMAP(oxygen);
    MMAP(hydrogen);
    MMAP(NO);
    MMAP(NH);
    MMAP(TI);
    MMAP(TB);
    MMAP(file);
    MMAP(line);
    MMAP(counter);
    MMAP(molecule);
    MMAP(notEnough);
    MMAP(notEnoughBeg);
    MMAP(atomsCreatingMol);
    MMAP(notCounter);
    MMAP(notCnt);

    
    return (oxygen == MAP_FAILED || hydrogen == MAP_FAILED || NO == MAP_FAILED || NH == MAP_FAILED || TB == MAP_FAILED || TI == MAP_FAILED || file == MAP_FAILED || line == MAP_FAILED || counter == MAP_FAILED || molecule == MAP_FAILED || notEnough
 == MAP_FAILED || notEnoughBeg == MAP_FAILED || atomsCreatingMol == MAP_FAILED || notCounter == MAP_FAILED || notCnt == MAP_FAILED);
}

void unmaping()
{
    UNMAP(oxygen);
    UNMAP(hydrogen);
    UNMAP(NO);
    UNMAP(NH);
    UNMAP(TI);
    UNMAP(TB);
    UNMAP(file);
    UNMAP(line);
    UNMAP(counter);
    UNMAP(molecule);
    UNMAP(notEnough);
    UNMAP(notEnoughBeg);
    UNMAP(atomsCreatingMol);
    UNMAP(notCounter);
    UNMAP(notCnt);

}


bool check_args(int argc, char ** argv)
{
    // num of params not 4
    if (argc != 5) {
        return false;
    }
    
    // verify NO
    char *scrapNO;
    long NO = strtol(argv[1], &scrapNO, 10);
    if (scrapNO[0] != '\0' || NO <= 0) {
        return false;
    }
    // verify NH
    char *scrapNH;
    long NH = strtol(argv[2], &scrapNH, 10);
    if (scrapNH[0] != '\0' || NH <= 0) {
        return false;
    }
    // verify TI
    long TI;
    char *scrapTI;
    TI = strtol(argv[3], &scrapTI, 10);
    if (*scrapTI != '\0' || scrapTI == argv[3] || !(TI >= 0 && TI <= 1000)) {
        return false;
    }
    // verify TB
    long TB;
    char *scrapTB;
    TB = strtol(argv[4], &scrapTB, 10);
    if (*scrapTB != '\0' || scrapTB == argv[4] || !(TB >= 0 && TB <= 1000)) {
        return false;
    }
    return true;
}

void error(char *errmsg)
{
    fprintf(stderr, "%s\n", errmsg);
    shutdown();
}

void shutdown()
{
    semaphores_destroy();
    unmaping();
    
    if (file != NULL)
    {
        fclose(file);
    }
}

bool semaphores_open()
{
    atMut = sem_open(atMut_name, O_CREAT, 0666, 1);
    oxyQueue = sem_open(oxyQueue_name, O_CREAT, 0666, 0);
    hydroQueue = sem_open(hydroQueue_name, O_CREAT, 0666, 0);
    turnstile = sem_open(turnstile_name, O_CREAT, 0666, 0);
    turnstile2 = sem_open(turnstile2_name, O_CREAT, 0666, 1);
    barrMut = sem_open(barrMut_name, O_CREAT, 0666, 1);
    outputMut = sem_open(outputMut_name, O_CREAT, 0666, 1);
    molMut = sem_open(molMut_name, O_CREAT, 0666, 1);
    sleepMut = sem_open(sleepMut_name, O_CREAT, 0666, 1);
    notTurnstile = sem_open(notTurnstile_name, O_CREAT, 0666, 0);
    notTurnstile2 = sem_open(notTurnstile2_name, O_CREAT, 0666, 1);
    notBarrMut = sem_open(notBarrMut_name, O_CREAT, 0666, 1);

    return (atMut == SEM_FAILED || oxyQueue == SEM_FAILED || hydroQueue == SEM_FAILED || turnstile == SEM_FAILED || turnstile2 == SEM_FAILED || barrMut == SEM_FAILED || outputMut == SEM_FAILED || molMut == SEM_FAILED || sleepMut == SEM_FAILED || notBarrMut == SEM_FAILED || notTurnstile == SEM_FAILED || notTurnstile2 == SEM_FAILED);
}

void semaphores_destroy()
{
    sem_close(atMut);
    sem_close(oxyQueue);
    sem_close(hydroQueue);
    sem_close(turnstile);
    sem_close(turnstile2);
    sem_close(barrMut);
    sem_close(outputMut);
    sem_close(molMut);
    sem_close(sleepMut);
    sem_close(notBarrMut);
    sem_close(notTurnstile);
    sem_close(notTurnstile2);
    
    sem_unlink(atMut_name);
    sem_unlink(oxyQueue_name);
    sem_unlink(hydroQueue_name);
    sem_unlink(turnstile_name);
    sem_unlink(turnstile2_name);
    sem_unlink(barrMut_name);
    sem_unlink(outputMut_name);
    sem_unlink(molMut_name);
    sem_unlink(sleepMut_name);
    sem_unlink(notBarrMut_name);
    sem_unlink(notTurnstile_name);
    sem_unlink(notTurnstile2_name);
    
}
