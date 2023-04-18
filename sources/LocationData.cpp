#include "LocationData.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <utility>

LocationData::LocationData(std::string location_name, long double latitude, long double longitude):
m_location_name{location_name}, m_latitude{latitude}, m_longitude{longitude} {
  //TODO - Constructor
}

LocationData::LocationData(const LocationData &other):
m_location_name{other.m_location_name}, m_latitude{other.m_latitude}, m_longitude{other.m_longitude} {
  //TODO - Copy Constructor
}
LocationData::LocationData(LocationData &&other) noexcept:
m_location_name{std::move(other.m_location_name)}, m_latitude{std::move(other.m_latitude)}, m_longitude{std::move(other.m_longitude)} {
  //TODO - Move Constructor
}

LocationData &LocationData::operator=(const LocationData &other) {
  //NOTE - Copy Assignment
  if(this != &other) {
    m_location_name = other.m_location_name;
    m_latitude      = other.m_latitude;
    m_longitude     = other.m_longitude;
  }
  return *this;
}
LocationData &LocationData::operator=(LocationData &&other) noexcept {
  //TODO - Move Assignment
  if(this != &other) {
    m_location_name = std::move(other.m_location_name);
    m_latitude      = std::move(other.m_latitude);
    m_longitude     = std::move(other.m_longitude);
  }
  return *this;  
}

//SECTION - Setters
auto LocationData::set_location_name(std::string location_name) -> void {
  m_location_name = std::move(location_name);
}
auto LocationData::set_latitude(long double latitude) -> void {
  m_latitude = latitude;
}
auto LocationData::set_longitude(long double longitude) -> void {
  m_longitude = longitude;
}

//SECTION - Getters
auto LocationData::get_location_name() const -> const std::string& {
  return m_location_name;
}
auto LocationData::get_latitude() const -> const long double {
  return m_latitude;
}
auto LocationData::get_longitude() const -> const long double {
  return m_longitude;
}

auto LocationData::to_radians(long double degrees) -> long double {
  return degrees * WGS84::PI / 180.0;
}

auto LocationData::vincenty_algorithm_inverse_geodetic_problem_WGS84(const LocationData &point) -> long double {

  long double f{WGS84::FLATTENING_WGS84_ELLIPSOID};                                //NOTE - Flattening of the ellipsoid
  long double a{WGS84::EARTH_RADIUS_WGS84};                                        //NOTE - Length of semi-major axis of the ellipsoid (radius at equator)
  long double b{a * (1 - f)};                                                      //NOTE - Length of semi-minor axis of the ellipsoid (radius at the poles)
  long double L{to_radians((point.get_longitude() - get_longitude()))};            //NOTE - Difference in longitude of the points on the auxiliary sphere
  long double U1{std::atan((1 - f) * std::tan(to_radians(get_latitude())))};       //NOTE - Reduced latitude (latitude on the auxiliary sphere)
  long double U2{std::atan((1 - f) * std::tan(to_radians(point.get_latitude())))}; //NOTE - Reduced latitude (latitude on the auxiliary sphere)
  long double sinU1{std::sin(U1)};
  long double sinU2{std::sin(U2)};
  long double cosU1{std::cos(U1)};
  long double cosU2{std::cos(U2)};
  long double lambda{L};
  long double lambda_p{2 * WGS84::PI};

  long double sin_lambda{}, 
              cos_lambda{}, 
              sin_sigma{}, 
              cos_sigma{}, 
              sigma{}, 
              sin_alpha{}, 
              cos_sq_alpha{},
              cos2_sigma_M{};

  uint32_t iteration_limit{100};

  while(std::abs(lambda - lambda_p) > std::numeric_limits<long double>::epsilon() && --iteration_limit > 0) {

    sin_lambda = std::sin(lambda);
    cos_lambda = std::cos(lambda);

    sin_sigma = std::sqrt(
      std::pow(cosU2 * sin_lambda, 2) + std::pow(cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda, 2)
    );

    //REVIEW - Co-incident Points
    if(sin_sigma == 0) {
      return 0.0L;
    }

    cos_sigma    = sinU1 * sinU2 + cosU1 * cosU2 * cos_lambda;
    sigma        = std::atan2(sin_sigma, cos_sigma);
    sin_alpha    = cosU1 * cosU2 * sin_lambda / sin_sigma;
    cos_sq_alpha = 1 - std::pow(sin_alpha, 2);

    //REVIEW - Equatorial Line
    (cos_sq_alpha == 0) ? cos2_sigma_M = 0 : cos2_sigma_M = cos_sigma - 2 * sinU1 * sinU2 / cos_sq_alpha;

    long double C = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha));
    lambda_p = lambda;
    lambda   = L + (1 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos2_sigma_M + C * cos_sigma * (-1 + 2 * std::pow(cos2_sigma_M, 2))));

  }

  if(iteration_limit == 0) {
    return std::numeric_limits<long double>::quiet_NaN();
  }

  long double u_sq{cos_sq_alpha * (a * a - b * b) / (b * b)};
  long double A{1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))};
  long double B{u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))};
  long double delta_sigma{
    B * sin_sigma * (
      cos2_sigma_M + B / 4 * (
        cos_sigma * (-1 + 2 * std::pow(cos2_sigma_M, 2)) - B / 6 * cos2_sigma_M * (-3 + 4 * std::pow(sin_sigma, 2)) * (-3 + 4 * std::pow(cos2_sigma_M, 2))
      )
    )
  };

  return (b * A * (sigma - delta_sigma));

}

auto LocationData::distance_to(const LocationData &point) -> long double {
  return vincenty_algorithm_inverse_geodetic_problem_WGS84(point) / 1e6L;
}

LocationData::~LocationData() {
  //TODO - Destructor empty
}
