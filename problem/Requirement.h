//
// Created by Alexandre Lemos on 12/08/2018.
//

#ifndef TRAIN_SCHEDULE_OPTIMISATION_REQUIREMENT_H
#define TRAIN_SCHEDULE_OPTIMISATION_REQUIREMENT_H

# include <string>

#include<time.h>
#include <stdio.h>
#include <string>
#include <ostream>
#include <list>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

class connection {
public:
    int id;
    std::string onto_section_marker;
    std::string min_connection_time;

    connection(int id, const std::string &onto_section_marker, const std::string &min_connection_time) : id(id),
                                                                                                         onto_section_marker(
                                                                                                                 onto_section_marker),
                                                                                                         min_connection_time(
                                                                                                                 min_connection_time) {}

    friend std::ostream &operator<<(std::ostream &os, const connection &connection1) {
        os << "id: " << connection1.id << " onto_section_marker: " << connection1.onto_section_marker
           << " min_connection_time: " << connection1.min_connection_time;
        return os;
    }
    std::string toString(){
        std::string s="";
        s.append("id: "+id);s.append(" onto_section_marker: ");s.append( onto_section_marker);
        s.append(" min_connection_time: " + min_connection_time);
        return s;
    }
};

class Requirement {
public:

    std::string id;
    std::string section_marker;
    std::string  type;
    std::string min_stopping_time;
    Requirement(std::string id, const std::string &section_marker,  const std::string & type, const std::string &min_stopping_time,
                const std::string &entry_earliest, std::string entry_delay_weight,
                const std::string &exit_earliest,const std::string &entry_latest,const std::string &exit_latest) : id(id),
                                                                                    section_marker(section_marker),
                                                                                    type(type), min_stopping_time(
                    min_stopping_time/*.substr(2,2)*/), entry_earliest(entry_earliest), entry_delay_weight(entry_delay_weight),
                                                                                    exit_earliest(exit_earliest),entry_latest(entry_latest),exit_latest(exit_latest)
                                                                                     {
                                                                                         if(entry_earliest!=""){
                                                                                             int h=0,m=0,s=0;
                                                                                             sscanf(entry_earliest.c_str(), "%d:%d:%d", &h, &m,&s);
                                                                                             sec_entry_earliest = (h * 60*60) + (m*60)+s;
                                                                                         }
                                                                                         if(entry_latest!=""){
                                                                                             int h=0,m=0,s=0;
                                                                                             sscanf(entry_latest.c_str(), "%d:%d:%d", &h, &m,&s);
                                                                                             sec_entry_latest = (h * 60*60) + (m*60)+s;
                                                                                         }
                                                                                         if(exit_latest!=""){
                                                                                             int h=0,m=0,s=0;
                                                                                             sscanf(exit_latest.c_str(), "%d:%d:%d", &h, &m,&s);
                                                                                             sec_exit_latest = (h * 60*60) + (m*60)+s;
                                                                                         }
                                                                                         if(exit_earliest!=""){
                                                                                             int h=0,m=0,s=0;
                                                                                             sscanf(exit_earliest.c_str(), "%d:%d:%d", &h, &m,&s);
                                                                                             sec_exit_earliest = (h * 60*60) + (m*60)+s;
                                                                                         }





                                                                                     }
    std::string entry_earliest;
    std::string entry_delay_weight;
    std::string exit_earliest;
    std::string  exit_latest;
    std::string  entry_latest;
    std::list<connection> connections;
    int sec_entry_earliest;
    int sec_exit_earliest,sec_entry_latest,sec_exit_latest;

    const std::list<connection, std::allocator<connection> > &getConnections() {
        return connections;
    }

    friend std::ostream &operator<<(std::ostream &os, const Requirement &requirement) {
        os << "id: " << requirement.id << " section_marker: " << requirement.section_marker << " type: "
           << requirement.type << " min_stopping_time: " << requirement.min_stopping_time << " entry_earliest: "
           << requirement.entry_earliest << " entry_delay_weight: " << requirement.entry_delay_weight
           << " exit_earliest: " << requirement.exit_earliest << " exit_latest: " << requirement.exit_latest
           << " entry_latest: " << requirement.entry_latest
           << " sec_entry_earliest: " << requirement.sec_entry_earliest << " sec_exit_earliest: "
           << requirement.sec_exit_earliest << " sec_entry_latest: " << requirement.sec_entry_latest
           << " sec_exit_latest: " << requirement.sec_exit_latest;
        return os;
    }

    void toString()
    {
        for (std::list<connection>::iterator it=connections.begin(); it != connections.end(); ++it)
            std::cout<<(*it).toString();
        std::cout<<std::endl;
    }


};



#endif //TRAIN_SCHEDULE_OPTIMISATION_REQUIREMENT_H
