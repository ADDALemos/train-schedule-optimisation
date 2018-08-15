//
// Created by Alexandre Lemos on 12/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
#define TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
#include "Requirement.h"
#include <list>


class Train {

public:
    int id;
    int route;
    std::list<Requirement> t;
};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TRAIN_H
