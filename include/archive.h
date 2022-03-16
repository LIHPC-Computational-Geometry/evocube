#pragma once

#pragma once

#include <queue>
#include <random>
#include <iostream>


template<typename T> 
class Archive {
private:
    std::deque<std::pair<double, T>> archive_;
    int archive_max_size_;
    double max_accepted_score_ = 10e10;

public:
    Archive(int archive_size);
    Archive(){};
    double requiredScore() const;
    bool insert(T individual, double score);
    void print() const;
    void clear();
    int probabilisticIndividual() const;
    void insertAll(const Archive<T> other_archive);

    // getters & setters
    T bestIndividual() const;
    double bestScore() const;
    int getSize() const;
    T getIndiv(int id) const;
    std::pair<T, double> popIndiv(int id);
    double getScore(int id) const;
    void setMaxAcceptScore(double max_accepted_score);
    int maxSize() const {return archive_max_size_;};
};

template<typename T> 
Archive<T>::Archive(int archive_size)
    : archive_max_size_(archive_size) {

}

template<typename T> 
void Archive<T>::setMaxAcceptScore(double max_accepted_score){
    max_accepted_score_ = max_accepted_score;
}

template<typename T> 
double Archive<T>::requiredScore() const {
    if (archive_.size() < archive_max_size_) return max_accepted_score_; //TODO
    return archive_.front().first;
}

template<typename T> 
bool Archive<T>::insert(T individual, double score){
    if (score >= requiredScore() && !(archive_.size() == 0)) return false;
    
    auto it = archive_.begin();
    while (it != archive_.end() && score <= (*it).first){
        if (individual == (*it).second) return false;
        it ++;
    }
    archive_.insert(it, std::pair<double, T>(score, individual));

    if (archive_.size() > archive_max_size_) archive_.pop_front();
    return true;
}

template<typename T> 
void Archive<T>::print() const {
    for (int i=0; i<archive_.size(); i++){
        std::cout << archive_[i].first << " "; // << archive_[i].second << "  ";
    }
    std::cout << std::endl;
}

template<typename T> 
void Archive<T>::insertAll(const Archive<T> other_archive){
    int s = other_archive.getSize();
    for (int i=0; i<s; i++){
        insert(other_archive.getIndiv(i), other_archive.getScore(i));
    }
}

template<typename T> 
void Archive<T>::clear() {
    archive_.clear();
}

template<typename T> 
T Archive<T>::bestIndividual() const {
    if (archive_.size() < 1) coloredPrint("Error: trying to get best individual from empty archive", "red");
    return archive_.back().second;
}

template<typename T> 
double Archive<T>::bestScore() const {
    if (archive_.size() < 1) coloredPrint("Error: trying to get best score from empty archive", "red");
    return archive_.back().first;
}

template<typename T> 
int Archive<T>::probabilisticIndividual() const {
    // Inspiration: Roulette wheel selection section from
    // "Towards Automatic Blocking of Shapes using Evolutionary Algorithm"

    if (archive_.size() == 0) {
        coloredPrint("Error: called probabilisticIndividual but archive is empty", "red");
    }

    double total_props = 0;
    int PN = archive_.size();
    for (int i=1; i<PN+1; i++) total_props += i;

    double random = (double) std::rand() / (RAND_MAX);
    double prob_sum = 0; // (do the math)
    for (int i=PN-1; i>0; i--) {
        if (random < prob_sum + static_cast<double>(i + 1) / total_props){
            return i;
            //return archive_[i].second;
        } 
        else {
            prob_sum += static_cast<double>(i + 1) / total_props;
        }
    }

    return 0;
    //return archive_.front().second;
}


template<typename T> 
int Archive<T>::getSize() const{
    return archive_.size();
}

template<typename T> 
T Archive<T>::getIndiv(int id) const {
    return archive_[id].second;
}

template<typename T> 
std::pair<T, double> Archive<T>::popIndiv(int id) {
    if (id >= archive_.size()){
        coloredPrint("Error: archive popping id too large", "red");
        return std::pair<T, double>();
    }
    std::pair<T, double> res = std::make_pair(archive_[id].second, archive_[id].first);
    archive_.erase(archive_.begin() + id); 
    return res;
}

template<typename T> 
double Archive<T>::getScore(int id) const {
    return archive_[id].first;
}