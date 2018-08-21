//
// Created by Alexandre Lemos on 11/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
#define TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
#include <iostream>
#include <list>
#include "Train.h"
#include "Route.h"


class Instance {
public:
    int hash;
    std::string label;
    std::string maxBandabweichung;
    std::list<Train> train;
    std::list<Route> route;
    std::list<Resource> resource;




};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
