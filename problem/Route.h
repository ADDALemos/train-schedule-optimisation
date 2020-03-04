//
// Created by Alexandre Lemos on 15/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_ROUTE_H
#define TRAIN_SCHEDULE_OPTIMISATION_ROUTE_H



#include <string>
#include <list>
#include "route_path.h"

class Route {


public:
    std::string id;
    int totalSeq;
    std::list<route_path> route_path;
};


#endif //TRAIN_SCHEDULE_OPTIMISATION_ROUTE_H
