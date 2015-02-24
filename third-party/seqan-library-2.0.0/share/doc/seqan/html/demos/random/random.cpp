#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/random.h>

using namespace seqan;

const int SEED = 42;

int main()
{
    Rng<MersenneTwister> rng(SEED);
    std::cout << "pickRandomNumber(rng) == " << pickRandomNumber(rng) << std::endl;

    typedef Value<Rng<MersenneTwister> >::Type TMTValue;
    TMTValue value = pickRandomNumber(rng);

    Pdf<Uniform<double> > uniformDouble(0, 1);
    std::cout << "pickRandomNumber(rng, uniformDouble) == " << pickRandomNumber(rng, uniformDouble) << std::endl;

    Pdf<Uniform<int> > uniformInt(0, 42);
    std::cout << "pickRandomNumber(rng, uniformInt) == " << pickRandomNumber(rng, uniformInt) << std::endl;

    Pdf<Normal> normal(0, 1);
    std::cout << "pickRandomNumber(rng, normal) == " << pickRandomNumber(rng, normal) << std::endl;

    Pdf<LogNormal> logNormal(0, 1);
    std::cout << "pickRandomNumber(rng, logNormal) == " << pickRandomNumber(rng, logNormal) << std::endl;

    Pdf<LogNormal> logNormal2(0, 1, MuSigma());
    std::cout << "pickRandomNumber(rng, logNormal2) == " << pickRandomNumber(rng, logNormal2) << std::endl;

    Pdf<LogNormal> logNormal3(0.1, 1, MeanStdDev());
    std::cout << "pickRandomNumber(rng, logNormal3) == " << pickRandomNumber(rng, logNormal3) << std::endl;

    CharString container = "Hello World!";
    shuffle(container, rng);
    std::cout << "shuffle(\"Hello World!\") == " << container << std::endl;

    return 0;
}
