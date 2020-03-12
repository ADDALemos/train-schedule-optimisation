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
#include "train_run_sections.h"

class Instance {
public:
    int hash;
    std::string label;
    std::string maxBandabweichung;
    std::vector<Train> train;
    std::map<std::string,Route> route;
    std::list<Resource> resource;
    std::map<std::string,std::map<int,route_section*>> sectionMap;
    std::map<std::string,std::vector<route_section*>> entryMap;
    std::map<std::string,std::vector<route_section*>> exitMap;
    std::map<std::string,std::vector<route_section*>> markerMap;
    std::map<std::string, std::map<int,std::vector<route_section*>>> end;
    std::map<std::string, double > route_pen;//
    std::map<std::string,std::map<int,train_run_sections*>> results;



};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
