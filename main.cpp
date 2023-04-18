#include <iostream>
#include "LocationData.hpp"
#include <matplot/matplot.h>
#include <string>
#include <vector>

//REVIEW - Vincenty's Formulae Algorithm Inverse Geodetic Problem
//LINK - https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Inverse_problem
//LINK - https://community.esri.com/t5/coordinate-reference-systems-blog/distance-on-an-ellipsoid-vincenty-s-formulae/ba-p/902053#:~:text=Vincenty%27s%20formulae%20were%20published%20in%201975,perhaps%2C%20more%20appropriate%20geometry%20to%20use&text=Vincenty%27s%20formulae%20were%20published,appropriate%20geometry%20to%20use&text=were%20published%20in%201975,perhaps%2C%20more%20appropriate%20geometry

//REVIEW - Latitude and Longitude of Tokyo, Japan
//LINK - https://www.latlong.net/place/tokyo-japan-8040.html

//REVIEW - Latitude and Longitude of Manila, Philippines
//LINK - https://www.latlong.net/place/manila-philippines-9339.html

auto main(int argc, char **argv) -> decltype(argc) {

  LocationData *japan_capital{new LocationData{"Tokyo, Japan", 35.652832L, 139.839478L}};
  LocationData *philippines_capital{new LocationData{"Manila, Philippines", 14.599512L, 120.984222L}};

  const long double distance{japan_capital->distance_to(*philippines_capital)};

  std::vector<long double> japan_philippines_longitude{japan_capital->get_longitude(), philippines_capital->get_longitude()};
  std::vector<long double> japan_philippines_latitude{japan_capital->get_latitude(), philippines_capital->get_latitude()};

  matplot::geoplot(japan_philippines_longitude, japan_philippines_latitude);
  matplot::hold(matplot::on);
  matplot::plot(japan_philippines_longitude, japan_philippines_latitude);
  matplot::title("Tokyo, Japan to Manila, Philippines Approximate Distance: " + std::to_string(distance) + " meters");
  
  matplot::show();

  matplot::save("../data/distance_tokyo_japan_manila_philippines.png");
  matplot::save("../data/distance_tokyo_japan_manila_philippines.svg");

  std::cout << "Approximate Distance: " << distance << " meters\n";

  delete japan_capital;
  japan_capital = nullptr;

  delete philippines_capital;
  philippines_capital = nullptr;

  return 0;
    
}
