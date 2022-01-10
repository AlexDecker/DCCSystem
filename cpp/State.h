#ifndef STATE_H
#define STATE_H

#include "ProblemParameters.h"

// https://github.com/jtlap/nt2

class State
{
    public:
        /**
         * @brief Create the initial state based on the problem parameters.
         */
        State(const ProblemParameters& params);
        /**
         * @brief Create the next state based on a state from the last time-slot,
         *        the transmitting voltages and the problem parameters.
         */
        State(const /*decidir o tipo das cargas*/, const ProblemParameters& params,
              const int timeSlot);
};

void worker(void* args);

#endif
