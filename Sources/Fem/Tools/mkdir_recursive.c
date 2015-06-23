
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <libgen.h>

int mkdir_recursive(char *path)
{
    char *subpath, *fullpath;
    mode_t mode = 0775;
    int ok = 0;
    struct stat fileStat;

    fullpath = strdup(path);
    subpath = dirname(path);
    if (strlen(subpath) > 1)
        ok = mkdir_recursive(subpath);

    if(stat(fullpath,&fileStat) < 0)
        if (mkdir(fullpath, mode) < 0)
          ok = -1;

    free(fullpath);
    return ok;
}


