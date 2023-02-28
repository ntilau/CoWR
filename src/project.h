#ifndef PROJECT_H
#define PROJECT_H

#define FES_VERSION 0.0.0.1

#include <string>
#include <time.h>
#include <filesystem>

#include "model.h"

class timer
{
public:
    timer()
    {
        tic();
    }
    ~timer() {}
    void tic()
    {
        lc = clock();
    }
    double toc()
    {
        clock_t cc = clock();
        return (cc - lc) / 1000.0;
    }
    std::string strtoc()
    {
        clock_t cc = clock();
        std::stringstream timing;
        timing << (cc - lc) / 1000.0 << " s";
        return timing.str();
    }

private:
    clock_t lc;
};

class logger : public std::ostringstream
{
public:
    template <typename T>
    logger &operator<<(T a)
    {
        oss << a;
        return *this;
    }
private:
    std::ostringstream oss;
};

class project : private std::filesystem::path,
                private timer
{
public:
    project(const std::string& name);
    ~project();
    logger log;
    // all the model information is stored here
    mdl_core model;

    std::string get_stats(timer &); // Statistics
    std::string get_info();         // User and computer names, Max RAM, CPU cores and threads number
    std::string get_proc_mem();
    std::string get_loc_time();

private:
    std::string get_var(const std::string name);
    int get_int(const std::string name);
};

#endif // PROJECT_H
