#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>

#include "dlib/matrix.h"
  
using dlib::abs;
using dlib::identity_matrix;
using dlib::is_finite;
using dlib::length;
using dlib::length_squared;
using dlib::matrix;
using dlib::normalize;
using dlib::ones_matrix;
using dlib::randm;
using dlib::squared;
using dlib::sum;
using dlib::pointwise_multiply;
using dlib::dot;
using dlib::tmp;
using dlib::trans;
using dlib::zeros_matrix;

using std::abs;
using std::array;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::get;
using std::getline;
using std::isfinite;
using std::isinf;
using std::isnan;
using std::make_tuple;
using std::min_element;
using std::round;
using std::sort;
using std::stod;
using std::stoull;
using std::streamsize;
using std::string;
using std::stringstream;
using std::swap;
using std::to_string;
using std::transform;
using std::tuple;
using std::vector;

template< class ForwardIt, class UnaryPredicate >
constexpr ForwardIt move_if(
    ForwardIt first,
    ForwardIt last,
    UnaryPredicate p ) {
  for (auto it=first;it!=last;) {
    if (p(*it)) {
      --last;
      swap(*it,*last);
    } else {
      ++it;
    }
  }
  return last;
}

struct filter_with_pivot{
  const matrix<double,2,1>& x; // the piovt
  const matrix<double,2,1>& d; // the derivative
  const bool operator()(
    const tuple<
      double,
      matrix<double,2,1>,
      matrix<double,2,1>,
      double
    > & v
      ) const {
    const matrix<double,2,1>& y = get<1>(v);
    return (trans(d)*(y-x))>0;
  }
};

tuple<double,matrix<double,2,1>> slow_centerpoint(
    const vector<tuple<double,matrix<double,2,1>>>& points
    ) {
  if (points.size()==0) {
    matrix<double,2,1> z = zeros_matrix<double>(2,1);
    return tuple<double,matrix<double,2,1>>{0.0,z};
  }
  double min_dist = std::numeric_limits<double>::infinity();
  auto median = points.begin();
  for (auto it0=points.begin();it0!=points.end();++it0) {
    double dist = 0; 
    for (auto it1=points.begin();it1!=points.end();++it1) {
      dist += get<0>(*it1)*length(get<1>(*it1)-get<1>(*it0));
    }
    if (dist<min_dist) {
      median = it0;
      min_dist = dist;
    }
  }
  return *median;
}

tuple<double,matrix<double,2,1>> fast_centerpoint(
    const vector<tuple<double,matrix<double,2,1>>>& points
    ) {
  if (points.size()==0) {
    return {0.0,zeros_matrix<double>(2,1)};
  }
  vector<
    tuple<
      double,
      matrix<double,2,1>,
      matrix<double,2,1>,
      double
    >
  > data;
  for (auto it=points.begin();it!=points.end();++it) {
    const double             w = get<0>(*it);
    const matrix<double,2,1> x = get<1>(*it);
    const matrix<double,2,1> d{0.0,0.0};
    data.emplace_back(w,x,d,0.0);
  }
  auto pivot = data.begin();
  auto pivot_end = data.begin()+1;
  auto end = data.end();
  while (true) {
    get<2>(*pivot) = zeros_matrix<double>(2,1);
    for (auto it=data.begin()+1;it!=data.end();++it) {
      get<2>(*pivot) += get<0>(*it)
                       *(get<1>(*pivot)-get<1>(*it))
                       /length(get<1>(*pivot)-get<1>(*it));
      get<3>(*pivot) +=get<0>(*it)*length(get<1>(*pivot)-get<1>(*it));
    }
    // filter the previous pivot elements and throw out the ones that are
    // not suitable pivots any more because they are clearly not median elements
    pivot_end = move_if(
        data.begin()+1,
        pivot_end,
        filter_with_pivot{get<1>(*data.begin()),get<2>(*data.begin())}
        );
    if (pivot_end==end) break;
    // filter all the other elements with the new pivot
    end = move_if(
        pivot_end,
        end,
        filter_with_pivot{get<1>(*data.begin()),get<2>(*data.begin())}
        );
    if (pivot_end==end) break;
    // move the pivot element to the front
    swap(*pivot_end,*data.begin());
    ++pivot_end;
    pivot = data.begin();
  }
  auto median = *min_element(
      data.begin(),
      pivot_end,[](auto & l,auto & r){
        return get<3>(l)<get<3>(r);
      });
  return {get<0>(median),get<1>(median)};
}

matrix<double,2,1> geometric_median(
    const vector<tuple<double,matrix<double,2,1>>>& points
    ) {
  if (points.size()==0) {
    return {0.0,zeros_matrix<double>(2,1)};
  }
  matrix<double,2,1> x = get<1>(fast_centerpoint(points));
  for (size_t i=0;i!=16;++i) {
    matrix<double,2,1> d1 = zeros_matrix<double>(2,1);
    matrix<double,2,2> d2 = zeros_matrix<double>(2,2);
    double sum = 0;
    double w = 0;
    for (auto it=points.begin();it!=points.end();++it) {
      const matrix<double,2,1> y = get<1>(*it);
      const double l2 = length_squared(x-y);
      const double l  = sqrt(l2);
      sum += get<0>(*it)*l;
      if (l2==0) {
        w = get<0>(*it);
        continue;
      }
      d1  += get<0>(*it)*(x-y)/l;
      d2  += get<0>(*it)*(identity_matrix<double>(2)-(x-y)*trans(x-y)/l2)/l;
    }
    //cout << endl;
    //cout << trans(d1);
    //cout << d2;
    //d1(0)-=w;
    //d1(1)-=w;
    //cout << trans(x);
    //cout << trans(d1);
    //cout << trans(inv(d2)*d1);
    //cout << length(d1) << endl;
    if (length(d1)<w) break;
    d1-=w*d1/length(d1);
    if (length(d1)<1e-8) break;
    matrix<double,2,1> t = x-inv(d2)*d1;
    while (true) {
      double sum2 = 0;
      for (auto it=points.begin();it!=points.end();++it) {
        const matrix<double,2,1> y = get<1>(*it);
        sum2 += get<0>(*it)*length(t-y);
      }
      if (sum2<=sum) {
        x = t;
        break;
      } else {
        t = 0.5*(x+t);
      }
      if (length(x-t)<1e-8) {
        x = t;
        break;
      }
    }
    //cout << trans(x) << endl;
  }
  return x;
}

// w * ( x - m ) / || x - m ||₂

int main(int argc,char** argv) {
  vector<tuple<double,matrix<double,2,1>>> points;
  for (string line;getline(cin,line);) {
    double w,x0,x1;
    stringstream ss(line);
    ss >> w >> x0 >> x1;
    if (!ss) break;
    matrix<double,2,1> x{x0,x1};
    points.emplace_back(w,x);
  }
  auto median = geometric_median(points);
  cout << trans(median);
  /*
  for (auto it=data.begin();it!=pivot_end;++it) {
    cout << get<0>(*it) << " " << trans(get<1>(*it));
  }
  */
}
