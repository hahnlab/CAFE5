#ifndef ROOT_DISTRIBUTION_H
#define ROOT_DISTRIBUTION_H

#include <map>
#include <vector>

class root_distribution
{
    std::vector<int> vectorized_dist;

    mutable bool max_calculated = false;
    mutable int max_value;
public:
    std::size_t size() const {
        return vectorized_dist.size();
    }
    int max() const;
    bool empty() const { return vectorized_dist.empty();  }

    void vector(const std::vector<int>& dist);
    void vectorize(const std::map<int, int>& rootdist);
    void vectorize_uniform(int max);
    void vectorize_increasing(int max);

    int sum() const;
    int at(size_t index) const;

    int select_randomly() const;

    void pare(size_t new_size);
};

#endif
