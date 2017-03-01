#pragma once

#include "base/interval.h"
#include "base/BinSettings.h"

#include <string>
#include <vector>

class TDirectory;
class TH1D;
class TH2D;
class TH3D;
class TTree;
class TGraph;
class TGraphErrors;

namespace ant {
namespace analysis {

class HistogramFactory {
private:

    TDirectory* my_directory = nullptr;
    void goto_dir() const;



    std::string title_prefix;

    std::string MakeTitle(const std::string& title) const;

    /**
     * @brief create a new TDirectory. If a TDirectory already exists, append a number (1, 2, 3, ...) at the end
     * @param name Name of the directory to create
     * @param rootdir where to create
     * @return new created directory
     */
    static TDirectory* mkDirNumbered(const std::string& name, TDirectory* rootdir);

    mutable unsigned n_unnamed = 0;
    std::string GetNextHistName(const std::string& name) const;


public:
    struct DirStackPush {
    private:
        TDirectory* dir;
    public:
        DirStackPush(const HistogramFactory& hf);
        ~DirStackPush();
    };


    HistogramFactory(const std::string& directory_name, TDirectory* root=nullptr, const std::string& title_prefix_ = "");
    HistogramFactory(const std::string& directory_name, const HistogramFactory &parent, const std::string& title_prefix_ = "");

    void SetRootDir(TDirectory* root_dir=nullptr);
    void SetTitlePrefix(const std::string& title_prefix_);
    std::string GetTitlePrefix() const;
    void SetDirDescription(const std::string& desc);

    TH1D* makeTH1D(
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& bins,
            const std::string& name="",
            bool  sumw2 = false) const;

    TH2D* makeTH2D(
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& xbins,
            const BinSettings& ybins,
            const std::string& name="",
            bool  sumw2 = false) const;

    TH3D* makeTH3D(const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const std::string& zlabel,
            const BinSettings& xbins,
            const BinSettings& ybins,
            const BinSettings& zbins,
            const std::string& name="",
            bool  sumw2 = false) const;

    TGraph* makeGraph(
            const std::string& title,
            const std::string& name="") const;

    TGraphErrors* makeGraphErrors(
            const std::string& title,
            const std::string& name="") const;

    TTree* makeTTree(const std::string& name) const;

    template<class T, typename... Args>
    T* make(Args&&... args) const {
        // save current dir and cd back to it on exit
        DirStackPush dirstack(*this);
        return new T(std::forward<Args>(args)...);
    }

    template <typename T>
    T* clone(const T* obj, const std::string& newName) const {
        DirStackPush dirstack(*this);
        return dynamic_cast<T*>(obj->Clone(newName.c_str()));
    }

};

}
}
