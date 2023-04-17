#include <iostream>
#include "LocationData.hpp"

//REVIEW - Vincenty's Formulae Algorithm Inverse Geodetic Problem
//LINK - https://en.wikipedia.org/wiki/Vincenty%27s_formulae#Inverse_problem
//LINK - https://community.esri.com/t5/coordinate-reference-systems-blog/distance-on-an-ellipsoid-vincenty-s-formulae/ba-p/902053#:~:text=Vincenty%27s%20formulae%20were%20published%20in%201975,perhaps%2C%20more%20appropriate%20geometry%20to%20use&text=Vincenty%27s%20formulae%20were%20published,appropriate%20geometry%20to%20use&text=were%20published%20in%201975,perhaps%2C%20more%20appropriate%20geometry

//REVIEW - Latitude and Longitude of Tokyo, Japan
//LINK - https://www.latlong.net/place/tokyo-japan-8040.html

//REVIEW - Latitude and Longitude of Manila, Philippines
//LINK - https://www.latlong.net/place/manila-philippines-9339.html

auto main(int argc, char **argv) -> decltype(argc) {

  LocationData *japan_capital{new LocationData{"Tokyo, Japan", 139.839478, 35.652832}};
  LocationData *philippines_capital{new LocationData{"Manila, Philippines", 120.984222, 14.599512}};

  double distance{japan_capital->distance_to(*philippines_capital)};

  std::cout << "Distance: " << distance << " meters\n";

  delete japan_capital;
  japan_capital = nullptr;

  delete philippines_capital;
  philippines_capital = nullptr;

  return 0;
    
}
