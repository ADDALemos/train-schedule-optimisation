//
// Created by Alexandre Lemos on 17/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_ROUTE_PATH_H
#define TRAIN_SCHEDULE_OPTIMISATION_ROUTE_PATH_H

#include <string>
#include <list>
#include "route_section.h"

class route_path {
public:
    std::string id;
    std::list<route_section> route_section;

};


#endif //TRAIN_SCHEDULE_OPTIMISATION_ROUTE_PATH_H
