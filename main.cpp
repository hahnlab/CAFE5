#include <random>

std::random_device rd;
std::mt19937 randomizer_engine(rd()); // seeding random number engine


int cafexp(int argc, char *const argv[]);

int main(int argc, char *const argv[]) {
    cafexp(argc, argv);
}
