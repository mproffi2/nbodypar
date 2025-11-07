#ifndef OMP_LOOP_HPP
#define OMP_LOOP_HPP

#include <omp.h>
#include <functional>
#include <cstddef>

class OmpLoop {
private:
    int nbThread = 1;
    int granularity = 1; 

public:
    OmpLoop() = default;

    void setNbThread(int n) {
        nbThread = n > 0 ? n : 1;
        omp_set_num_threads(nbThread);
    }

    void setGranularity(int g) {
        granularity = g > 0 ? g : 1;
    }

    void parfor(size_t begin, size_t end, const std::function<void(size_t)> &func) const {
    #pragma omp parallel for schedule(static, 1) num_threads(nbThread)
        for (size_t i = begin; i < end; i++) {
            func(i);
        }
    }
};

#endif 