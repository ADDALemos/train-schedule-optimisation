//
// Created by Alexandre Lemos on 13/03/2020.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_STATS_H
#define TRAIN_SCHEDULE_OPTIMISATION_STATS_H

#include "Instance.h"

void size(Instance *i){
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
    printf("Number of Sections: %d\n",sec);
    printf("Number of Resources: %d\n",res);

}


#endif //TRAIN_SCHEDULE_OPTIMISATION_STATS_H
