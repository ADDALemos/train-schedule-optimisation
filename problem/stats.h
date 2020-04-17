//
// Created by Alexandre Lemos on 13/03/2020.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_STATS_H
#define TRAIN_SCHEDULE_OPTIMISATION_STATS_H

#include "Instance.h"

void stat(Instance instance, int diff){
    int res=0;int sec=0;
    for (int j = 0; j < instance.train.size(); ++j) {
        int s = 0;
        res += instance.train[j].t.size();
        for (route_path rp: instance.route[instance.train[j].route].route_path) {
            for (route_section *rs: rp.route_section) {
                sec++;
            }
        }
    }
    printf("Number of Trains: %d\n",instance.train.size());
    printf("Number of Sections: %d\n",sec);
    printf("Number of Resources: %d\n",res);
    printf("MaxTime: %d\n",diff);
    printf("MaxBandwitch: %s\n",instance.maxBandabweichung.c_str());
    printf("Number of Resources: %d\n",instance.resource.size());



}


#endif //TRAIN_SCHEDULE_OPTIMISATION_STATS_H
