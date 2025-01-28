#include <random>
#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

std::random_device rd;
std::mt19937 randomizer_engine(rd()); // seeding random number engine


int cafe5(int argc, char *const argv[]);

int main(int argc, char *const argv[]) {

    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Global, el::ConfigurationType::Format, "%msg");
    defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, "WARNING: %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);

    cafe5(argc, argv);
}
