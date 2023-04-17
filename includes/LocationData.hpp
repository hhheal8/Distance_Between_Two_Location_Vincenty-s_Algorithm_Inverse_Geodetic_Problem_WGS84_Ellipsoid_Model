#ifndef LOCATIONDATA_HPP
#define LOCATIONDATA_HPP

#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>
#include <utility>
#include <unordered_map>
#include <unordered_set>

namespace WGS84 {
  constexpr double PI{3.14159265358979323846};
  constexpr double EARTH_RADIUS_WGS84{6378137.0};                   //NOTE - WGS-84 Ellipsoid Model in meters 
  constexpr double FLATTENING_WGS84_ELLIPSOID{1.0 / 298.257223563}; //NOTE - Formula for flattening the WGS-84 ellipsoid model
};

class LocationData {

  private:
    std::string m_location_name{};
    double m_longitude{}, m_latitude{};

    auto to_radians(double degrees) -> double;
    auto vincenty_algorithm_inverse_geodetic_problem_WGS84(const LocationData &point) -> double;

  public:
    LocationData() = default;
    LocationData(std::string location_name, double longitude, double latitude);

    LocationData(const LocationData &other);
    LocationData(LocationData &&other) noexcept;

    LocationData &operator=(const LocationData &other);
    LocationData &operator=(LocationData &&other) noexcept;

    //SECTION - setter
    auto set_location_name(const std::string &location_name) -> void;
    auto set_longitude(const double &longitude) -> void;
    auto set_latitude(const double &latitude) -> void;

    //SECTION - getter
    auto get_location_name() const -> const std::string&;
    auto get_longitude() const -> const double&;
    auto get_latitude() const -> const double&;

    auto distance_to(const LocationData &point) -> double;

    ~LocationData();

};

#endif