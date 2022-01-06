#ifndef WORKER_H
#define WORKER_H

#include "Ring.h"
#include "FeasibleFuture.h"

class Memory
{
    public:
        Memory();
        void runCycle();
    protected:
        Ring m_states;
    private:
        Ring& m_pastStates;
        FeasibleFuture m_feasibleFuture;
};

void worker(void* args);

#endif
