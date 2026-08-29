#pragma once
#include <cstdint>

using ObjectIdT = std::uint64_t;

inline ObjectIdT NextObjectId() {
    static ObjectIdT counter = 0;
    return counter++;
}