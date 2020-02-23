//
// Created by Alexandre Lemos on 15/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_RESOURCE_H
#define TRAIN_SCHEDULE_OPTIMISATION_RESOURCE_H

#include <iostream>
#include <string>


class Resource {

public:
    Resource(std::string id, std::string release_time, bool following_allowed) : id(id), release_time(release_time.substr(2,2)),
                                                                                 following_allowed(following_allowed) {
    }
    Resource(std::string id, std::string occupation_direction) : id(id), occupation_direction(occupation_direction){
    }
    Resource(std::string id) : id(id){
    }
    Resource(){
    }
    friend std::ostream &operator<<(std::ostream &os, const Resource &Resource) {
        os << "id: " << Resource.id << " release_time: " << Resource.release_time << " following_allowed: "
           << Resource.following_allowed;
        return os;
    }


private:
    std::string id;
    std::string release_time;
    bool following_allowed = false;
    std::string occupation_direction;


};


#endif //TRAIN_SCHEDULE_OPTIMISATION_RESOURCE_H
