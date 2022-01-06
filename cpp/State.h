#ifndef STATE_H
#define STATE_H

// https://github.com/jtlap/nt2

class State
{
    public:
        State();
        void runCycle();
    protected:
        Ring m_states;
    private:
        Ring& m_pastStates;
        FeasibleFuture m_feasibleFuture;
};

void worker(void* args);

#endif
