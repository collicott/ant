#pragma once

#include "base/vec3.h"

// ROOT compat
#include "TLorentzVector.h"


namespace ant {

/**
 * @brief A Lorentz vector (x,y,z,E), diag(-1,-1,-1, 1)
 */
struct LorentzVec {
    vec3   p = {};
    double E = {};

    LorentzVec() noexcept = default;
    LorentzVec(const LorentzVec&) noexcept = default;
    LorentzVec(LorentzVec&&) noexcept = default;
    LorentzVec(const vec3& p_, const double E_) noexcept:
        p(p_), E(E_) {}
    LorentzVec(const double px, const double py, const double pz, const double _E) noexcept :
        p(px,py,pz), E(_E) {}

    LorentzVec& operator=(const LorentzVec&) noexcept = default;
    LorentzVec& operator=(LorentzVec&&) noexcept = default;


    // ====== TLorentzVector interface ======

    operator TLorentzVector() const {
        return TLorentzVector(p, E);
    }

    LorentzVec(const TLorentzVector& other) noexcept:
        p(other.Vect()), E(other.E()) {}

    LorentzVec& operator=(const TLorentzVector other) noexcept {
        p = other.Vect();
        E = other.E();
        return *this;
    }

    // ======================================

    static LorentzVec EPThetaPhi(const double E, const double p, const double theta, const double phi) noexcept {
        return LorentzVec(vec3::RThetaPhi(p,theta,phi), E);
    }

    LorentzVec& operator+=(const LorentzVec& other) noexcept {
        p += other.p;
        E += other.E;
        return *this;
    }

    LorentzVec& operator-=(const LorentzVec& other) noexcept {
        p -= other.p;
        E -= other.E;
        return *this;
    }

    LorentzVec operator+(const LorentzVec& other) const noexcept {
        return LorentzVec(*this) += other;
    }

    LorentzVec operator-(const LorentzVec& other) const noexcept {
        return LorentzVec(*this) -= other;
    }

    double M2() const noexcept {
        return E*E - p.R2();
    }

    double M() const {
        const auto mm = this->M2();
        return mm < 0.0 ? -sqrt(-mm) : sqrt(mm);
    }

    double Theta() const {
        return p.Theta();
    }

    double Phi() const {
        return p.Phi();
    }

    double P() const {
        return p.R();
    }

    bool operator==(const LorentzVec& other) const noexcept {
        return p == other.p && E == other.E;
    }

    bool operator!=(const LorentzVec& other) const noexcept {
        return !(*this == other);
    }

    LorentzVec& operator*=(const double a) noexcept {
        p *= a;
        E *= a;
        return *this;
    }

    LorentzVec& operator/=(const double a) noexcept {
        p /= a;
        E /= a;
        return *this;
    }

    LorentzVec operator* (const double a) const noexcept {
        return LorentzVec(*this)*=a;
    }

    LorentzVec operator/ (const double a) const noexcept {
        return LorentzVec(*this)/=a;
    }

    double Beta() const noexcept {
        return p.R() / E;
    }

    /**
     * @brief Gamma = 1 / sqrt(1-beta^2)
     * @return
     */
    double Gamma() const noexcept {
        return 1.0 / sqrt(1 - (p.R2()/(E*E)));
    }

    double Dot(const LorentzVec& other) const noexcept {
        return E*other.E - p.Dot(other.p);
    }

    vec3 BoostVector() const noexcept {
        return p / E;
    }

    double Angle(const vec3& other) const {
        return this->p.Angle(other);
    }

    double Angle(const LorentzVec& other) const {
        return this->p.Angle(other.p);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(p, E);
    }

};

}