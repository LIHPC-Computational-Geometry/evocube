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

    double evaluate(const LabelingIndividual& indiv, double& fast_poly_score, double& invalid_score, int& n_fail_invert) const {
        fast_poly_score = qle_->evaluate(indiv.getLabeling(), n_fail_invert); // avg disto + flipped triangles
        invalid_score = indiv.invalidityScore();
        return  fast_poly_score + 100.0 * invalid_score; // + 10000.0 * invalidChartsScore(adj);
    }

    double evaluate(const LabelingIndividual& indiv) const {
        double fast_poly_score, invalid_score;
        int n_fail_invert;
        return evaluate(indiv, fast_poly_score, invalid_score, n_fail_invert);
    }
private:
    std::shared_ptr<Evocube> evo_;
    std::shared_ptr<const QuickLabelEv> qle_;
};