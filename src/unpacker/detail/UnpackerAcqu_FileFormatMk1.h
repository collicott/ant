#pragma once

#include "UnpackerAcqu_detail.h"

namespace ant {

struct TDAQError;

namespace unpacker {
namespace acqu {

class FileFormatMk1 : public FileFormatBase {
public:
    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

protected:

    std::list<unsigned> ScalerBlockSizes;

    virtual size_t SizeOfHeader() const override;
    virtual bool InspectHeader(const std::vector<std::uint32_t>& buffer) const override;
    virtual void FillInfo(reader_t& reader, buffer_t& buffer, Info& info) override;
    virtual void FillFirstDataBuffer(reader_t& reader, buffer_t& buffer) const override;
    virtual bool UnpackDataBuffer(queue_t& queue, it_t& it, const it_t& it_endbuffer) noexcept override;

    void FindScalerBlocks(const std::vector<std::string>& scaler_modnames);

    void UnpackEvent(queue_t& queue, it_t it, const it_t& it_end, bool& good) noexcept;
    void HandleDAQError(std::vector<TDAQError>& errors,
                        it_t& it, const it_t& it_end, bool& good) const noexcept;
    void HandleScalerBuffer(scalers_t& scalers,
                            it_t& it, const it_t& it_end, bool& good) const noexcept;
};

}}} // namespace ant::unpacker::acqu
