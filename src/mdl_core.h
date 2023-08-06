#ifndef MDL_CORE_H
#define MDL_CORE_H

#include <cstring>
#include <string>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>


#include "mdl_frm.h"
#include "mdl_sld.h"
#include "mdl_msh.h"
 
using namespace std;

class mdl_core
{
public:
    mdl_frm frm;
    mdl_sld sld;
    mdl_msh msh;
    void create_tri_mesh();
    void clear()
    {
        sld.clear();
        frm.clear();
        msh.clear();
    }

protected:
    void removeCharsFromString(string &str, const char *charsToRemove)
    {
        for (unsigned int i = 0; i < strlen(charsToRemove); ++i)
        {
            str.erase(remove(str.begin(), str.end(), charsToRemove[i]), str.end());
        }
    }
    char find_SI_factor(string &str)
    {
        for (unsigned int i = 0; i < strlen(SI_chars); ++i)
        {
            size_t found = str.find(SI_chars[i]);
            if (found <= str.size())
                return str[found];
        }
        return 0;
    }
    double set_factor(char fact)
    {
        if (fact == 'm')
            return 1e-3;
        else if (fact == 'u')
            return 1e-6;
        else if (fact == 'n')
            return 1e-9;
        else if (fact == 'p')
            return 1e-12;
        else
            return 1;
    }
    const char *SI_chars = "munp";
};

#endif // MDL_CORE_H
