//
// Created by Alexandre Lemos on 11/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
#define TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include "Train.h"
#include "Route.h"

using namespace std;
class Instance {
public:
    int hash;
    std::string label;
    std::string maxBandabweichung;
    std::vector<Train> train;
    std::map<std::string,Route> route;
    std::list<Resource> resource;


    vector<map<string,vector<pair<route_section, route_section>>>> pair;
};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
