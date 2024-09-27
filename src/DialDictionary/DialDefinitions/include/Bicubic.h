#ifndef GUNDAM_BICUBIC_INTERPOLATION_H
#define GUNDAM_BICUBIC_INTERPOLATION_H

#include "DialBase.h"
#include "DialUtils.h"
#include "DialInputBuffer.h"

#include <vector>
#include <utility>

/// Manage a bicubic interpolation for a grid of points.
class Bicubic : public DialBase {

public:
  Bicubic() = default;
  ~Bicubic() override = default;

  [[nodiscard]] std::unique_ptr<DialBase> clone() const override {
       return std::make_unique<Bicubic>(*this);
  }

  [[nodiscard]] std::string getDialTypeName() const override {return {"Bicubic"};}

  /// Allow extrapolation of the data.  The default is to
  /// forbid extrapolation.
  virtual void setAllowExtrapolation(bool allow_) override;
  [[nodiscard]] virtual bool getAllowExtrapolation() const override;

  // Return the dial response for the input parameters.  The DialInputBuffer
  // should contain two input parameters.
  [[nodiscard]] double evalResponse(const DialInputBuffer& input_) const override;

  /// Pass information to the dial so that it can build it's internal
  /// information.
  virtual void buildDial(const TH2& h2, const std::string& option_="") override;

  [[nodiscard]] const std::vector<double>& getDialData() const override {return _splineData_;}

protected:
  bool _allowExtrapolation_{false};

  // A block of data to calculate the spline values.  This must be filled for
  // the Cache::Manager to work, and provides the input for spline calculation
  // functions that can be shared between the CPU and the GPU.
  std::vector<double> _splineData_{};

  // The vector of input parameter bounds.
  std::vector<DialUtils::Range> _splineBounds_;
};

typedef CachedDial<Bicubic> BicubicCache;


#endif //GUNDAM_BICUBIC_SPLINE_H
