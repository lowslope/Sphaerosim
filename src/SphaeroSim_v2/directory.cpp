#include "directory.h"

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>

#ifdef __linux__

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#else

#include <tchar.h>

#endif

#include "Exception.h"

Directory::Directory() : elements() {
    current_directory = _T("");
}

#ifdef __linux__

Directory::Directory(const char *path) : elements(), it() {
    current_directory = _T("");
    std::string param(path);
    if (!param.empty()) {
        if (param[param.size() - 1] != '/') {
            param += '/';
        }
    }
    DIR *dirp;
    dirent *entry;
    dirp = opendir(param.c_str());
    if (dirp != nullptr) {
        entry = readdir(dirp);
        while (entry != nullptr) {
            this->elements.push_back(*entry);
            entry = readdir(dirp);
        }
        closedir(dirp);
        it = elements.begin();
        current_directory = param;
    }
}

#else

Directory::Directory(const TCHAR* path)
{
  current_directory = _T("");

  std::basic_string<TCHAR> param;
  this->make_windows_compatible(path, param);
  HANDLE dirp = INVALID_HANDLE_VALUE;
  WIN32_FIND_DATA entry;
  dirp = FindFirstFile(param.c_str(), &entry);
  if(dirp != INVALID_HANDLE_VALUE)
  {
    elements.push_back(entry);
    while(FindNextFile(dirp, &entry) != 0)
      elements.push_back(entry);
    FindClose(dirp);
    param.erase(param.size() - 1, 1);
    current_directory.assign(param);
  }
}

#endif

Directory::Directory(const Directory &rhs) : elements(rhs.elements) {
    current_directory = rhs.current_directory;
}

Directory::~Directory() = default;

Directory &Directory::operator=(const Directory &rhs) {
    if (this == &rhs) {
        return *this;
    }
    elements = rhs.elements;
    return *this;
}

#ifdef __linux__

bool Directory::open_dir(const char *path) {
    std::string param(path);
    if (!param.empty()) {
        if (param[param.size() - 1] != '/') {
            param += '/';
        }
    }
    DIR *dirp;
    dirent *entry;
    dirp = opendir(param.c_str());
    if (dirp != nullptr) {
        entry = readdir(dirp);
        while (entry != nullptr) {
            this->elements.push_back(*entry);
            entry = readdir(dirp);
        }
        closedir(dirp);
        it = elements.begin();
        current_directory = param;
        return true;
    }
    return false;
}

#else

bool Directory::open_dir(const TCHAR* path)
{
  //erst schliessen
  elements.clear();
  std::basic_string<TCHAR> param;
  this->make_windows_compatible(path, param);
  HANDLE dirp = INVALID_HANDLE_VALUE;
  WIN32_FIND_DATA entry;
  dirp = FindFirstFile(param.c_str(), &entry);
  if(dirp != INVALID_HANDLE_VALUE)
  {
    elements.push_back(entry);
    while(FindNextFile(dirp, &entry) != 0)
      elements.push_back(entry);
    FindClose(dirp);
    param.erase(param.size() - 1, 1);
    current_directory.assign(param);
    return true;
  }
  return false;
}
#endif // __linux__


#ifdef __linux__

bool Directory::compare(dirent entry_one, dirent entry_two) {
    return (strcmp(entry_one.d_name, entry_two.d_name) < 0);
}

#else
bool Directory::compare(const WIN32_FIND_DATA &entry_one, const WIN32_FIND_DATA &entry_two)
{
  return (strcmp(entry_one.cFileName, entry_two.cFileName) < 0);
}
#endif


#ifdef __linux__

bool Directory::compare_length(dirent entry_one, dirent entry_two) {
    return strlen(entry_one.d_name) < strlen(entry_two.d_name);
}

#else
bool Directory::compare_length(const WIN32_FIND_DATA &entry_one, const WIN32_FIND_DATA &entry_two)
{
  return _tcslen(entry_one.cFileName) < _tcslen(entry_two.cFileName);
}
#endif

#ifndef __linux__
void Directory::make_windows_compatible(const TCHAR *path, std::basic_string<TCHAR>& out)
{
  out.assign(path);
  if(out.size() < 2)
    return;
  if(out[out.size() - 1] != '*')
  {
    if(out[out.size() - 1] != '\\')
    {
      out += '\\';
    }
    out += '*';
    return;
  }
}
#endif


bool Directory::file_exists(const std::string &filename) {
    std::ifstream file_s(filename.c_str());
    bool ret_val = file_s.is_open();
    if (ret_val) {
        file_s.peek();
        ret_val = file_s.good();
        if (ret_val) {
            file_s.close();
        }
    }
    return ret_val;
}


#ifdef __linux__

bool Directory::create_directory(const std::string &path) {
    //const mode_t dir_mode = 0775;
    //system(("mkdir -p "+path).c_str());
    return (system(("mkdir -p " + path).c_str()) == 0);
}

#else
bool Directory::create_directory(const std::string &path)
{
  const LPSECURITY_ATTRIBUTES attributes = 0;

  // create all the subfoler
  std::vector<std::string> subfolder;
  std::size_t start_pos = path.find('/') + 1;
  std::size_t end_pos;
  do {
    end_pos = path.find('/', start_pos);
    if (end_pos != std::string::npos)
      subfolder.push_back(path.substr(0, end_pos));
    else
      subfolder.push_back(path.substr(0));
    if (end_pos != std::string::npos)
      start_pos = end_pos + 1;
  } while (end_pos != std::string::npos);

  for (std::size_t cur = 0; cur < subfolder.size(); ++cur) {
    if (Directory::folderExists(subfolder[cur]) == false) {
      if (CreateDirectory(subfolder[cur].c_str(), attributes) == 0)
        throw SphaeroSim::Exception("Failed to create the directory.",
                                    path);
    }
  }
  return true;
}
#endif


// return the filename for an element
std::string Directory::getFileNameFromIndex(int iIndex) {
#ifdef __linux__
    return elements[iIndex].d_name;
#else
    return elements[iIndex].cFileName;
#endif
}

// returns true if the element iIndex is a folder
bool Directory::isFolder(int iIndex) {
#ifdef __linux__
    struct stat st;
    lstat(elements[iIndex].d_name, &st);
    return (S_ISDIR(st.st_mode));
#else
    if ((elements[iIndex].dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0)
      return false;
    return true;
#endif
}

// returns true if pathname exists as a folder
bool Directory::folderExists(const std::string &strPathname) {
#ifdef __linux__
    DIR *dir = opendir(strPathname.c_str());
    if (dir) {
        // Directory exists.
        closedir(dir);
        return true;
    } else if (ENOENT == errno) {
        // Directory does not exist.
        return false;
    } else {
        //opendir() failed for some other reason.
        throw "Internal IO Error";
    }
#else
    DWORD dwAttributes = GetFileAttributes(strPathname.c_str());
    if(dwAttributes == INVALID_FILE_ATTRIBUTES)
      return false; // does not exist
    if ((dwAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0)
      return false;
    return true;
#endif
}
