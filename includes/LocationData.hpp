#ifndef LOCATIONDATA_HPP
#define LOCATIONDATA_HPP

#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace WGS84 {
  constexpr long double PI{3.141592653589793238462643};
  constexpr long double EARTH_RADIUS_WGS84{6378137.0};                   //NOTE - WGS-84 Ellipsoid Model in meters 
  constexpr long double FLATTENING_WGS84_ELLIPSOID{1.0 / 298.257223563}; //NOTE - Formula for flattening the WGS-84 ellipsoid model
};

class LocationData {

  private:
    std::string m_location_name{};
    long double m_latitude{}, m_longitude{};

    auto to_radians(long double degrees) -> long double;
    auto vincenty_algorithm_inverse_geodetic_problem_WGS84(const LocationData &point) -> long double;

  public:
    LocationData() = default;
    LocationData(std::string location_name, long double latitude, long double longitude);

    LocationData(const LocationData &other);
    LocationData(LocationData &&other) noexcept;

    LocationData &operator=(const LocationData &other);
    LocationData &operator=(LocationData &&other) noexcept;

    //SECTION - Setter
    auto set_location_name(std::string location_name) -> void;
    auto set_latitude(long double latitude) -> void;
    auto set_longitude(long double longitude) -> void;

    //SECTION - Getter
    auto get_location_name() const -> const std::string&;
    auto get_latitude() const -> const long double;
    auto get_longitude() const -> const long double;

    auto distance_to(const LocationData &point) -> long double;

    ~LocationData();

};

#endif