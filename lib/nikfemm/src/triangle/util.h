#ifndef UTIL_H
#define UTIL_H

extern "C" {
#include "../../lib/triangle/triangle.h"
}

void syntax();
void info();

void report(struct triangulateio * io, int markers, int reporttriangles, int reportneighbors, int reportsegments, int reportedges, int reportnorms);

#endif /* UTIL_H */