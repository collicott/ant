#pragma once

#include <memory>
#include <list>

#include "Reconstruct_traits.h"

namespace ant {

struct TTaggerHit;

namespace reconstruct {
class CandidateBuilder;
class Clustering_traits;
class UpdateableManager;
}

class Reconstruct : public Reconstruct_traits {
public:
    Reconstruct();

    // this method converts a TDetectorRead
    // into a calibrated TEvent
    virtual void DoReconstruct(TEventData& reconstructed) override;

    ~Reconstruct();

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };



protected:

    bool initialized = false;
    bool includeIgnoredElements = false;

    virtual void Initialize(const TID& tid);

    using sorted_readhits_t = ReconstructHook::Base::readhits_t;
    sorted_readhits_t sorted_readhits;

    void ApplyHooksToReadHits(std::vector<TDetectorReadHit>& detectorReadHits);

    template<typename T>
    using sorted_bydetectortype_t = std::map<Detector_t::Type_t, std::vector< T > >;

    void BuildHits(sorted_bydetectortype_t<TClusterHit>& sorted_clusterhits,
            std::vector<TTaggerHit>& taggerhits
            );

    void HandleTagger(const std::shared_ptr<TaggerDetector_t>& taggerdetector,
            const std::vector<std::reference_wrapper<TDetectorReadHit>>& readhits,
            std::vector<TTaggerHit>& taggerhits);

    using sorted_clusterhits_t = ReconstructHook::Base::clusterhits_t;
    using sorted_clusters_t = ReconstructHook::Base::clusters_t;
    void BuildClusters(const sorted_clusterhits_t& sorted_clusterhits,
                       sorted_clusters_t& sorted_clusters);

    // little helper class which stores the upcasted versions of shared_ptr
    // to Detector_t instances
    struct detector_ptr_t {
        std::shared_ptr<Detector_t> Detector;
        std::shared_ptr<TaggerDetector_t> TaggerDetector; // might be nullptr
        std::shared_ptr<ClusterDetector_t> ClusterDetector; // might be nullptr

        detector_ptr_t(const std::shared_ptr<Detector_t>& detector) :
            Detector(detector),
            TaggerDetector(std::dynamic_pointer_cast<TaggerDetector_t>(detector)),
            ClusterDetector(std::dynamic_pointer_cast<ClusterDetector_t>(detector))
        {
            if(TaggerDetector != nullptr && ClusterDetector != nullptr) {
                throw Exception("Found detector which is both clustering and tagging, not supported");
            }
        }
        // implicit conversion to simple base class Detector pointer
        operator std::shared_ptr<Detector_t>() const { return Detector; }
    };
    struct sorted_detectors_t : std::map<Detector_t::Type_t,  detector_ptr_t > {
        // use the implicit conversions of detector_ptr_t
        operator std::map<Detector_t::Type_t,  std::shared_ptr<Detector_t> >() const {
            return std::map<Detector_t::Type_t,  std::shared_ptr<Detector_t> >(begin(), end());
        }
    };
    sorted_detectors_t sorted_detectors;

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;
    shared_ptr_list<ReconstructHook::DetectorReadHits> hooks_readhits;
    shared_ptr_list<ReconstructHook::ClusterHits>      hooks_clusterhits;
    shared_ptr_list<ReconstructHook::Clusters>         hooks_clusters;
    shared_ptr_list<ReconstructHook::EventData>        hooks_eventdata;


    std::unique_ptr<const reconstruct::CandidateBuilder>  candidatebuilder;
    std::unique_ptr<const reconstruct::Clustering_traits> clustering;
    std::unique_ptr<reconstruct::UpdateableManager> updateablemanager;
};

}
