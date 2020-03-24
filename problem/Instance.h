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
    std::map<std::string,std::map<int,route_section*>> sectionMap;//train id section id section
    std::map<std::string,std::map<std::string,std::map<int,route_section*>>> pathMap;//train id path id section
    std::map<std::string,std::vector<Resource*>> markerRMap;//Resource marker Resource
    std::map<std::string,std::vector<route_section*>> entryMap;//entry marker+trainID section
    std::map<std::string,std::vector<route_section*>> exitMap;//exit marker+trainID section
    std::map<std::string,std::vector<route_section*>> markerMap;//Resource marker section
    std::map<std::string, std::map<int,std::vector<route_section*>>> end;//train end nodes sections
    std::map<std::string, double > route_pen;//

    std::map<std::string,std::map<int,train_run_sections*>> results;

    //solution
    int solution_hash;



};


#endif //TRAIN_SCHEDULE_OPTIMISATION_TIMETABLE_H
