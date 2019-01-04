#include "root_distribution.h"
#include <numeric>
#include <algorithm>
#include <random>

extern std::mt19937 randomizer_engine; // seeding random number engine

//! Used when simulating gene families (-s)
void root_distribution::vector(const std::vector<int>& dist)
{
    vectorized_dist = dist;
    max_calculated = false;
}

void root_distribution::vectorize(const std::map<int, int>& rootdist) {

    for (auto it = rootdist.begin(); it != rootdist.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            vectorized_dist.push_back(it->first);
        }
    }
    max_calculated = false;
}

void root_distribution::vectorize_uniform(int max) {
    vectorized_dist.resize(max);
    std::fill(vectorized_dist.begin(), vectorized_dist.end(), 1);
    max_calculated = false;

}

void root_distribution::vectorize_increasing(int max) {
    vectorized_dist.resize(max);

    int n = 0;
    std::generate(vectorized_dist.begin(), vectorized_dist.end(), [&n]() mutable { return n++; });
    max_calculated = false;
}

int root_distribution::max() const
{
    if (!max_calculated)
    {
        max_value = *std::max_element(vectorized_dist.begin(), vectorized_dist.end());
        max_calculated = true;
    }

    return max_value;
}

int root_distribution::sum() const
{
    if (vectorized_dist.empty())
        throw std::runtime_error("Root distribution not created yet");

    return std::accumulate(vectorized_dist.begin(), vectorized_dist.end(), 0);
}

int root_distribution::at(size_t val) const
{
    if (val >= vectorized_dist.size())
        throw std::out_of_range("Root distribution value out of range");

    return vectorized_dist.at(val);
}


/* START: Drawing random root size from uniform */
template<typename itr, typename random_generator>
itr select_randomly(itr start, itr end, random_generator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g)); // advances iterator (start) by dis(g), where g is a seeded mt generator
    return start;
}

template<typename itr>
itr select_randomly(itr start, itr end) {
    static std::random_device rd; // randomly generates uniformly distributed ints (seed)
    static std::mt19937 gen(rd()); // seeding Mersenne Twister generator
    return select_randomly(start, end, gen); // plug in mt generator to advance our container and draw random element from it
}
/* END: Drawing random root size from uniform */


int root_distribution::select_randomly() const
{
    return *::select_randomly(vectorized_dist.begin(), vectorized_dist.end()); // getting a random root size from the provided (core's) root distribution
}

void root_distribution::pare(size_t new_size)
{
    if (vectorized_dist.size() < new_size) return;

    shuffle(vectorized_dist.begin(), vectorized_dist.end(), randomizer_engine);
    vectorized_dist.erase(vectorized_dist.begin()+new_size, vectorized_dist.end());
    sort(vectorized_dist.begin(), vectorized_dist.end());
    max_calculated = false;
}
