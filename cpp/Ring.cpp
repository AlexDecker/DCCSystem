#include "Ring.h"

Ring::Ring(const size_t ringSize, const size_t vectorSize)
{}

void Ring::push()
{
    if (!full())
    {
        ++m_last;
        if (m_last == m_data.end())
        {
            m_last = m_data.begin();
        }
    }
}

void Ring::pop()
{
    if (!empty())
    {
        ++m_first;
        if (m_first == m_data.end())
        {
            m_first = m_data.begin();
        }
    }
}

const State& Ring::first() const
{
    return *m_first;
}
        
State& Ring::last()
{
    return *m_last;
}

bool Ring::full()
{}
        
bool Ring::empty()
{}

static bool test()
{
    const int ringSize = 4;
    const int vectorSize = 3;
    Ring ring(ringSize, vectorSize);
    
    for (int i = 0; i <= 10;
}