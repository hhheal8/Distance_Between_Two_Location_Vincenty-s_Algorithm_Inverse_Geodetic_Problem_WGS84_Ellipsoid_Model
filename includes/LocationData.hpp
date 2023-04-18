#ifndef LOCATIONDATA_HPP
#define LOCATIONDATA_HPP

#include <cmath>
#include <limits>
#include <stdexcept>
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

    const std::string error_location_name{
      "\nLocationData(!!!INVALID FORMAT/INPUT!!!, " + std::to_string(get_latitude()) + ", " + std::to_string(get_longitude()) + ")\nLocation name should not be empty.\n"
    };
    const std::string error_latitude{
      "\nLocationData(" + get_location_name() + ", !!!INVALID FORMAT/INPUT!!!, " + std::to_string(get_longitude()) + ")\nLatitude value should be within the range of -90 to 90 degrees.\n"
    };
    const std::string error_longitude{
      "\nLocationData(" + get_location_name() + ", " + std::to_string(get_latitude()) + ", !!!INVALID FORMAT/INPUT!!!)\nLongitude value should be within the range of -180 to 180 degrees.\n"
    };

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