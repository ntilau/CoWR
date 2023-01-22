#include "model.h"
#include "project.h"
#include <iostream>

int main(int argc, char* argv[]) {
    timer t;
    try{
        project p (argc, argv);
    }
    catch(std::string err)
    {
        std::cout << err << std::endl;
        std::cout << "Wall time: " << t.strtoc() << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Wall time: " << t.strtoc() << std::endl;
    return EXIT_SUCCESS;
}
