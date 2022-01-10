#ifndef RING_H
#define RING_H

#include State.h

class Ring
{
    public:
        Ring(const size_t ringSize, const size_t vectorSize);

        /**
         * @brief Creates a new state to be the "last()".
         */
        void push();
        
        /**
         * @brief Eliminates the oldest state ("first()").
         */
        void pop();

        /**
         * @brief The state to be read.
         */
        const State& first() const;
        
        /**
         * @brief The state to be written.
         */
        State& last();
        
        bool full();
        
        bool empty();
        
        /**
         * @brief Unit test for this class.
         */
        static bool test();

    private:
        std::vector<State>::iterator m_first;
        std::vector<State>::iterator m_last;
        std::vector<State> m_data;
};

#endif