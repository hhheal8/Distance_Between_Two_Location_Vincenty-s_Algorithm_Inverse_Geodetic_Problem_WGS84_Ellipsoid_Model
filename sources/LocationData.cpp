#include "LocationData.hpp"

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <utility>

LocationData::LocationData(std::string location_name, double longitude, double latitude):
m_location_name{location_name}, m_longitude{longitude}, m_latitude{latitude} {

}

LocationData::LocationData(const LocationData &other):
m_location_name{other.m_location_name}, m_longitude{other.m_longitude}, m_latitude{other.m_latitude} {

}
LocationData::LocationData(LocationData &&other) noexcept:
m_location_name{std::move(other.m_location_name)}, m_longitude{std::move(other.m_longitude)}, m_latitude{std::move(other.m_latitude)} {

}

LocationData &LocationData::operator=(const LocationData &other) {
  if(this != &other) {
    m_location_name = other.m_location_name;
    m_longitude     = other.m_longitude;
    m_latitude      = other.m_latitude;
  }
  return *this;
}
LocationData &LocationData::operator=(LocationData &&other) noexcept {
  if(this != &other) {
    m_location_name = std::move(other.m_location_name);
    m_longitude     = std::move(other.m_longitude);
    m_latitude      = std::move(other.m_latitude);
  }
  return *this;  
}

auto LocationData::set_location_name(const std::string &location_name) -> void {
  m_location_name = std::move(location_name);
}
auto LocationData::set_longitude(const double &longitude) -> void {
  m_longitude = std::move(longitude);
}
auto LocationData::set_latitude(const double &latitude) -> void {
  m_latitude = std::move(latitude);
}

auto LocationData::get_location_name() const -> const std::string& {
  return m_location_name;
}
auto LocationData::get_longitude() const -> const double& {
  return m_longitude;
}
auto LocationData::get_latitude() const -> const double& {
  return m_latitude;
}

auto LocationData::to_radians(double degrees) -> double {
  return degrees * WGS84::PI / 180.0;
}

auto LocationData::vincenty_algorithm_inverse_geodetic_problem_WGS84(const LocationData &point) -> double {

  double f{WGS84::FLATTENING_WGS84_ELLIPSOID};                                //NOTE - Flattening of the ellipsoid
  double a{WGS84::EARTH_RADIUS_WGS84};                                        //NOTE - Length of semi-major axis of the ellipsoid (radius at equator)
  double b{a * (1 - f)};                                                      //NOTE - Length of semi-minor axis of the ellipsoid (radius at the poles)
  double L{to_radians((point.get_longitude() - get_longitude()))};            //NOTE - Difference in longitude of the points on the auxiliary sphere
  double U1{std::atan((1 - f) * std::tan(to_radians(get_latitude())))};       //NOTE - Reduced latitude (latitude on the auxiliary sphere)
  double U2{std::atan((1 - f) * std::tan(to_radians(point.get_latitude())))}; //NOTE - Reduced latitude (latitude on the auxiliary sphere)
  double sinU1{std::sin(U1)};
  double sinU2{std::sin(U2)};
  double cosU1{std::cos(U1)};
  double cosU2{std::cos(U2)};
  double lambda{L};
  double lambda_p{2 * WGS84::PI};

  double sin_lambda{}, 
         cos_lambda{}, 
         sin_sigma{}, 
         cos_sigma{}, 
         sigma{}, 
         sin_alpha{}, 
         cos_sq_alpha{},
         cos2_sigma_M{};

  uint32_t iteration_limit{100};

  while(std::abs(lambda - lambda_p) > std::numeric_limits<double>::epsilon() && --iteration_limit > 0) {

    sin_lambda = std::sin(lambda);
    cos_lambda = std::cos(lambda);

    sin_sigma = std::sqrt(
      std::pow(cosU2 * sin_lambda, 2) + std::pow(cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda, 2)
    );

    //REVIEW - Co-incident Points
    if(sin_sigma == 0) {
      return 0;
    }

    cos_sigma    = sinU1 * sinU2 + cosU1 * cosU2 * cos_lambda;
    sigma        = std::atan2(sin_sigma, cos_sigma);
    sin_alpha    = cosU1 * cosU2 * sin_lambda / sin_sigma;
    cos_sq_alpha = 1 - std::pow(sin_alpha, 2);

    //REVIEW - Equatorial Line
    (cos_sq_alpha == 0) ? cos2_sigma_M = 0 : cos2_sigma_M = cos_sigma - 2 * sinU1 * sinU2 / cos_sq_alpha;

    double C = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha));
    lambda_p = lambda;
    lambda   = L + (1 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos2_sigma_M + C * cos_sigma * (-1 + 2 * std::pow(cos2_sigma_M, 2))));

  }

  if(iteration_limit == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double u_sq{cos_sq_alpha * (a * a - b * b) / (b * b)};
  double A{1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))};
  double B{u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))};
  double delta_sigma{B * sin_sigma * (cos2_sigma_M + B / 4 * (cos_sigma * (-1 + 2 * std::pow(cos2_sigma_M, 2)) - B / 6 * cos2_sigma_M * (-3 + 4 * std::pow(sin_sigma, 2)) * (-3 + 4 * std::pow(cos2_sigma_M, 2))))};

  return b * A * (sigma - delta_sigma);

}

auto LocationData::distance_to(const LocationData &point) -> double {
  return vincenty_algorithm_inverse_geodetic_problem_WGS84(point);
}

LocationData::~LocationData() {

}
