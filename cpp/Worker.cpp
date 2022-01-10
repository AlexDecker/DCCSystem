#include "Worker.h"

Memory::Memory()
{}

void Memory::runCycle()
{
    // Pick a state from the last FeasibleFuture
    std::optional<State> state = m_pastStates.next();
    
    // There is no new initial state now
    if (!state)
    {
        // Skip to next thread
        return;
    }
    
    int ttl = m_sampleSize;
    while (--ttl && !m_states.full())
    {
        // Random state which is reachable from state at this timeslot
        const State nextState(state->charges, m_problemParameters, m_timeSlot);
        
        // Store only the next states which generate a substantial contribution
        // to the next FeasibleFuture
        if (m_feasibleFuture.contribution(nextState) > m_minimalContribution) {
            m_feasibleFuture.add(nextState);
            m_states.add(nextState);
        }
    }
}

void worker(void* args)
{
    Memory* memory = static_cast<Memory*> (args);
}
