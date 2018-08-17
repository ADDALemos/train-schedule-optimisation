//
// Created by Alexandre Lemos on 17/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_ROUTE_SECTION_H
#define TRAIN_SCHEDULE_OPTIMISATION_ROUTE_SECTION_H
#include <string>
#include <list>
#include "Resource.h"



class route_section {
public:
    int sequence_number;
    std::list<std::string> route_alternative_marker_at_entry;
    std::list<std::string> section_marke;
    std::list<Resource> resource_occupations;
    double penalty=0;
    std::string starting_point;
    std::string minimum_running_time;
    std::string ending_point;


};


#endif //TRAIN_SCHEDULE_OPTIMISATION_ROUTE_SECTION_H
