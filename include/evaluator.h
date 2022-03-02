#pragma once
#include <memory>

#include "labeling_individual.h"
#include "evocube.h"
#include "quick_label_ev.h"


class Evaluator {
public:
    Evaluator(std::shared_ptr<Evocube> evo, 
              std::shared_ptr<const QuickLabelEv> qle) 
        : evo_(evo), qle_(qle){};

    double evaluate(const LabelingIndividual& indiv) const {
        double fast_poly_score = qle_->evaluate(indiv.getLabeling()); // avg disto + flipped triangles
        double invalid_score = indiv.invalidChartsScore();
        return  fast_poly_score + 10.0 * invalid_score; // + 10000.0 * invalidChartsScore(adj);
    }
private:
    std::shared_ptr<Evocube> evo_;
    std::shared_ptr<const QuickLabelEv> qle_;
};