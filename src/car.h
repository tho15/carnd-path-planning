#ifndef CAR_H
#define CAR_H

#include<math.h>

using namespace std;

class Car {
public:
	Car() { };
	Car(int id, double x, double y, double vx, double vy, double s, double d):
		car_id(id), car_x(x), car_y(y), car_vx(vx), car_vy(vy), car_s(s), car_d(d) { };
	
	~Car() { };
	
	int get_car_lane() {
		int i = 0;
		if(car_d > 8.0) i = 2;
		else if (car_d > 4.0) i = 1;
		
		return i;
	}
	
	double get_car_land_d() {
		double l = 2;
		if(car_d > 8.0) l = 10;
		else if(car_d > 4.0) l = 6;
		
		return l;
	}
	
	double get_car_vel() const {
		return sqrt(car_vx*car_vx + car_vy*car_vy)*0.02;
	}
	
	void print() const {
		cout << "car: " << car_id << " " << car_x << " " << car_y << " " << car_s << " " << " " << car_d << endl;
	}
	
	int		car_id;
	double	car_x;
	double	car_y;
	double	car_vx;
	double	car_vy;
	double	car_s;
	double	car_d;
};





#endif // CAR_H
