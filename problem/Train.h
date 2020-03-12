//
// Created by Alexandre Lemos on 12/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
#define TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
#include "Requirement.h"
#include "route_section.h"
#include <vector>


class Train {

public:
    std::string id;
    std::string route;
    std::vector<Requirement*> t;
};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
