#pragma once

#include "date_utils.h"

namespace voltlbx
{
    
    class BusinessTimeMeasure : public Pimpl<BusinessTimeMeasure>
    {
    public:
        BusinessTimeMeasure(std::shared_ptr<Calendar> calendar,
                            double close_weight,
                            double yearly_nb_open);

        double distance(const chrono::DateTime& t0, const chrono::DateTime& t1) const;
    };

}