#include "helper_functions.hpp"

#include <random>
#include "../../include/pcg_random.hpp"


double getRandom(double min, double max) 
{
    // Seed from true rng once
    static pcg_extras::seed_seq_from<std::random_device> seed_source;
    
    // Make a random engine
    static thread_local pcg32_fast rng(seed_source);
    // std::uniform_real_distribution<> distrib(min, max);
    // return distrib(rng);


    std::uniform_real_distribution<> distrib(0, 1);
    return (max - min) * distrib(rng) + min;
}

void loadingBar(int Nbins, int i, int j)
{
    int totalBins = Nbins*Nbins;

    int binNum = i*Nbins + j;

    std::cout << "Finished " << binNum << "/" << totalBins << "\r";
}