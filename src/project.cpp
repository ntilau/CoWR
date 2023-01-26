#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
#include <thread>
#include <unistd.h>

#elif _WIN32
#include <iphlpapi.h>
#include <psapi.h>
#include <windows.h>
#include <winsock2.h>

#else
#error "OS not supported!"
#endif

#include "project.h"
#include "solver.h"

project::project(int argc, char **argv) // parse and execute
{
  if (argc > 1) {
    full_path_name = argv[1];
    name = full_path_name;
    task = LOAD_CDNS;
    execute_task();
  }

  // task = LOAD_POLY;
  // execute_task();
  // int ret = system(std::string("tetgen -pqA " + full_path_name +
  // ".poly").data()); task = RUN_TETGEN; execute_task();
  // task = SAVE_FES;
  // execute_task();

  // model.msh.save_vtk_mesh(name);
  model.msh.get_mesh_statistics();

  // std::cout << get_info() << std::endl;
  // std::cout << get_proc_mem() << std::endl;
}

project::~project() {}

void project::execute_task() {
  try {
    switch (task) {
    case LOAD_FES:
      std::cout << "Loading " << name << ".fes\n";
      model.frm.read_prj_file(full_path_name);
      model.sld.read_prj_file(full_path_name);
      model.msh.read_prj_file(full_path_name);
      break;
    case SAVE_FES:
      std::cout << "Saving " << name << ".fes\n";
      model.frm.write_prj_file(full_path_name);
      model.sld.write_prj_file(full_path_name);
      model.msh.write_prj_file(full_path_name);
      break;
    case LOAD_POLY:
      std::cout << "Loading " << name << ".poly\n";
      model.sld.read_poly_file(full_path_name);
      break;
    case SAVE_POLY:
      break;
    case LOAD_STL:
      std::cout << "Loading " << name << ".stl\n";
      model.sld.read_stl_file(full_path_name);
      model.create_tri_mesh();
      break;
    case LOAD_HFSS:
      std::cout << "Loading " << name << ".hfss\n";
      model.wrap_hfss(full_path_name, aux_path);
      break;
    case LOAD_AEDT:
      std::cout << "Loading " << name << ".aedt\n";
      model.wrap_aedt(full_path_name, aux_path);
      break;
    case LOAD_CDNS:
      std::cout << "Loading " << name << ".mesh\n";
      model.wrap_cdns(full_path_name, aux_path);
      break;
    case RUN_TETGEN:
      std::cout << "Loading " << name << ".poly products\n";
      model.msh.read_tetgen_files(full_path_name);
      break;
    case RUN_TRIANGLE:
      std::cout << "Loading " << name << ".poly products\n";
      model.msh.read_triangle_files(full_path_name);
      break;
    case REFINE_HOMOGENEOUSLY:
      std::cout << "Refining homogeneously\n";
      model.msh.refine_homogeneous();
      break;
    case NONE:
      std::cout << "Nothing to do\n";
      break;
    }
    if (task == ANALYZE) {
      std::cout << "Running solver\n";
      // solver sol(model);
    }
  } catch (std::string &str) {
    std::cout << "Error: " << str << "\n";
  }
}

void project::parser(int, char **) {}

std::string project::get_info() {
  std::stringstream tag;
  std::string name, cores, threads;
#ifdef __linux__
  char hostname[HOST_NAME_MAX];
  char username[LOGIN_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  getlogin_r(username, LOGIN_NAME_MAX);
  name = std::string(hostname);
  cores = std::to_string(std::thread::hardware_concurrency());
  threads = std::to_string(std::thread::hardware_concurrency());
#elif _WIN32
  name = get_var("COMPUTERNAME");
  cores = get_var("NUMBER_OF_PROCESSORS");
  threads = get_var("OMP_NUM_THREADS");
#endif
  if (name.size() != 0) {
    tag << "COMPUTERNAME         = " << name;
  } else {
    tag << "COMPUTERNAME         = ?";
  }
  tag << "\n";
  if (cores.size() != 0) {
    tag << "NUMBER_OF_PROCESSORS = " << cores;
  } else {
    tag << "NUMBER_OF_PROCESSORS = ?";
  }
  tag << "\n";
  if (threads.size() != 0) {
    tag << "OMP_NUM_THREADS      = " << threads;
  } else {
    tag << "OMP_NUM_THREADS      = "
        << cores; // automatically setting to the number of cores
  }
  return tag.str();
}

int project::get_num_proc() {
  unsigned int processor_count = 0;
#ifdef __linux__
  processor_count = std::thread::hardware_concurrency();
#elif _WIN32
  processor_count = atoi(get_var("NUMBER_OF_PROCESSORS").data());
#endif
  return processor_count;
}

#ifdef _WIN32
std::string project::get_var(const std::string name) {
  char *ptr = getenv(name.c_str());
  std::string ret;
  if (ptr == NULL) {
    ret = std::string("");
  } else {
    ret = std::string(ptr);
  }
  return ret;
}

int project::get_int(const std::string name) {
  const std::string data = get_var(name);
  int ret = -1;
  if (data.size() != 0) {
    ret = atoi(data.c_str());
  }
  return ret;
}

std::string project::set_priority(unsigned int lvl) {
  // setpriority(PRIO_PROCESS, 0, -20);
  HANDLE process = GetCurrentProcess();
  switch (lvl) {
  case 0:
    SetPriorityClass(process, NORMAL_PRIORITY_CLASS);
    break;
  case 1:
    SetPriorityClass(process, HIGH_PRIORITY_CLASS);
    break;
  case 2:
    SetPriorityClass(process, REALTIME_PRIORITY_CLASS);
    break;
  default:
    SetPriorityClass(process, NORMAL_PRIORITY_CLASS);
  }
  return get_priority();
}

std::string project::get_priority() {
  DWORD dwPriClass = GetPriorityClass(GetCurrentProcess());
  std::stringstream out;
  if (dwPriClass == REALTIME_PRIORITY_CLASS) {
    out << "REALTIME";
  } else if (dwPriClass == HIGH_PRIORITY_CLASS) {
    out << "HIGH";
  } else if (dwPriClass == NORMAL_PRIORITY_CLASS) {
    out << "NORMAL";
  } else {
    out << "level = " << dwPriClass;
  }
  out << " priority process";
  return out.str();
}
#endif

std::string project::get_proc_mem() {
  std::stringstream out;
  double physiPeak, physiPres;
#ifdef __linux__
  long pages = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret = getrusage(who, &usage);
  physiPeak = pages * page_size / 1048576;
  physiPres = usage.ru_maxrss / 1024;
#elif _WIN32
  HANDLE hProcess = GetCurrentProcess();
  if (NULL == hProcess) {
    out << "Memory stats: Failed to acquire process handle";
    return out.str();
  }
  PROCESS_MEMORY_COUNTERS pmc;
  if (!GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
    out << "Memory stats: Failed to acquire process memory information";
  } else {
    physiPeak = (double)pmc.PeakWorkingSetSize / 1048576;
    physiPres = (double)pmc.WorkingSetSize / 1048576;
  }
  CloseHandle(hProcess);
#endif
  out << "Used RAM: " << physiPres << " MB |" << physiPeak << "|";
  return out.str();
}

std::string project::get_sys_mem() {
  std::stringstream out;
#ifdef __linux__
  out << "linux" << std::endl;
#elif _WIN32
  MEMORYSTATUSEX statex;
  statex.dwLength = sizeof(statex);
  if (GlobalMemoryStatusEx(&statex)) {
    out << "Free RAM: ";
    out << statex.ullAvailPhys / 1048576;
    out << " MB";
  } else {
    out << "No memory information available";
  }
  out << "\n";
#else
#error "OS not supported!"
#endif
  return out.str();
}

std::string project::get_loc_time() {
  time_t ct = time(NULL);
  return std::string(asctime(localtime(&ct)));
}

std::string project::get_gmt() {
  time_t ct = time(NULL);
  return std::string(asctime(gmtime(&ct)));
}

std::string project::get_stats(timer &t) {
  return "- " /*+ get_proc_mem() + " - " */ + t.strtoc() + "\n";
}
