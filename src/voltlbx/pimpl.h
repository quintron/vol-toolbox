#pragma once
#include <memory>

namespace voltlbx
{
    template<typename T>
    class Pimpl
    {
    public:
        struct Implementation;

        Implementation& implementation()
        {
            return *impl.get();
        }

        const Implementation& implementation() const
        {
            return *impl.get();
        }

    protected:
        std::shared_ptr<Implementation> impl;

        Pimpl() = default;

        template<typename... Args>
        Pimpl(Args&&... args)
            : impl(new Implementation(std::forward<Args>(args)...))
        {}
    };

}