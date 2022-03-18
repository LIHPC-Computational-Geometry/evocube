#pragma once

/**
 * @brief Logging utilities
 * @date 2021-10-12
 * 
 */

#include <string>
#include <vector>
#include <map>
#include <chrono>

using std::chrono::milliseconds;
using std::chrono::duration_cast;

void coloredPrint(std::string text, std::string color);

void resetLogFile(std::string filepath);

void fillLogInfo(std::string tag, std::string filepath, std::string value);
void fillLogInfo(std::string tag1, std::string tag2, std::string filepath, std::string value);
void fillLogInfo(std::string tag1, std::string tag2, std::string filepath, double value);

void removeLogInfo(std::string tag, std::string filepath);

std::string readLogInfo(std::string tag, std::string filepath);

std::map<std::string, std::string> readAllLogInfo(std::string filepath);

std::vector<std::string> getLogTags(std::string filepath);

void printLog(std::string filepath);

double measureTime(std::chrono::time_point<std::chrono::steady_clock> t1,
                   std::chrono::time_point<std::chrono::steady_clock> t2);

class Chronometre {
public:
    Chronometre(std::string log_file)
            : log_file_(log_file) {
        initial_time_ = std::chrono::steady_clock::now();
        last_lap_ = std::chrono::steady_clock::now();
    };

    void newTimestamp(std::string lap_name){
        std::chrono::steady_clock::time_point new_lap = std::chrono::steady_clock::now();
        int64_t lap_time = duration_cast<milliseconds>(new_lap - last_lap_).count();
        fillLogInfo(lap_name, log_file_, std::to_string(lap_time));
        last_lap_ = new_lap;
    }

    void totalTime(){
        std::chrono::steady_clock::time_point new_lap = std::chrono::steady_clock::now();
        int64_t total_time = duration_cast<milliseconds>(new_lap - initial_time_).count();
        fillLogInfo("TimeTotal", log_file_, std::to_string(total_time));
    }
    
private:
    std::string log_file_;
    std::chrono::steady_clock::time_point last_lap_;
    std::chrono::steady_clock::time_point initial_time_;
};