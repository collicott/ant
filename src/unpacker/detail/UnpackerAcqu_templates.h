#pragma once

#include "UnpackerAcqu_detail.h"
#include "UnpackerAcqu_legacy.h"

#include <memory>
#include <vector>
#include <type_traits>

// only to be included in implementation files...

namespace ant {
namespace unpacker {

// inspectHeader can actually be implemented for Mk1 and Mk2...
template<typename T>
bool inspectHeaderMk1Mk2_(const T* h, std::true_type) {
    return h->fMk2 != acqu::EHeadBuff;
}

template<typename T>
bool inspectHeaderMk1Mk2_(const T*, std::false_type) {
    return false;
}

template<typename T>
bool inspectHeaderMk1Mk2(const std::vector<uint32_t>& buffer) {
    // ensure the correct T at compile time
    static_assert(std::is_same<T, acqu::AcquExptInfo_t>::value
                  || std::is_same<T, acqu::AcquMk2Info_t>::value,
                  "T can only be either Mk1 or Mk2 header struct");

    if(buffer[0] != acqu::EHeadBuff)
        return false;

    // try to interpret the buffer as some Mk2 header
    const T* h = reinterpret_cast<const T*>(buffer.data()+1);

    // Mk2 has some additional head marker in the struct,
    // but the other type does not have it
    // we use tag-based dispatching here:
    // if T is AcquMk2Info_t, then call the checkMk2 which inspects the struct,
    // if T is something else, then call the empty
    using tag = std::integral_constant<bool, std::is_same<T, acqu::AcquMk2Info_t>::value >;
    if(inspectHeaderMk1Mk2_(h, tag{}))
        return false;

    /// \todo Improve the Mk1/Mk2 header checking here

    if(std_ext::string_sanitize(h->fTime).length() != 24)
        return false;

    if(std_ext::string_sanitize(h->fOutFile).empty())
        return false;

    return true;
}

}} // namespace ant::unpacker
