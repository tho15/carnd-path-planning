#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/LU"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#include "car.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

vector<double> map_waypoints_x;
vector<double> map_waypoints_y;
vector<double> map_waypoints_s;
vector<double> map_waypoints_dx;
vector<double> map_waypoints_dy;

vector<double> map_waypoints_fit_s;
vector<double> map_waypoints_fit_x;
vector<double> map_waypoints_fit_y;

enum {
	_CS_KL,
	_CS_KL_FL,
	_CS_LCL,
	_CS_LCR
};

vector<int>  CAR_STATE;  // { state, target lane }

double max_s = 6945.554;
double max_speed = 0.445;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

vector<double> JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */
    MatrixXd A(3, 3);
    A << pow(T, 3), pow(T, 4), pow(T, 5),
         3*pow(T, 2), 4*pow(T, 3), 5*pow(T, 4),
         6*T, 12*pow(T, 2), 20*pow(T, 3);
    
    VectorXd B(3);
    B[0] = end[0] - (start[0]+start[1]*T+start[2]*pow(T, 2)/2.0);
    B[1] = end[1] - (start[1]+start[2]*T);
    B[2] = end[2] - start[2];
    
    VectorXd x(3);
    x = A.inverse()*B;
         
    return {start[0],start[1],start[2]/2.0,x[0],x[1],x[2]};
}

double quintic_eval_s(vector<double> p, double t)
{
	return p[0]+p[1]*t+p[2]*pow(t,2)+p[3]*pow(t,3)+p[4]*pow(t,4)+p[5]*pow(t,5);
}

double quintic_eval_d(vector<double> p, double t)
{
	return p[1]+2*p[2]*t+3*p[3]*pow(t,2)+4*p[4]*pow(t,3)+5*p[5]*pow(t,4);
}

double quintic_eval_a(vector<double> p, double t)
{
	return 2*p[2]+6*p[3]*t+12*p[4]*pow(t, 2)+20*p[5]*pow(t, 3);
}

void max_jmt_acc_jerk(vector<double> r, double dt, double &max_acc, double &max_jerk)
{
	double ts = 1.0/dt;
	
	double t_am1 = (-24.0*r[4]+sqrt(pow(24*r[4], 2) - 4*6*r[3]*60*r[5]))/(2*60*r[5]);
	double t_am2 = (-24.0*r[4]-sqrt(pow(24*r[4], 2) - 4*6*r[3]*60*r[5]))/(2*60*r[5]);
	//cout << "max acc at: " << t_am1 << " " << t_am2 << endl;
	
	double max_a1 = 2*r[2]+6*r[3]*t_am1+12*r[4]*pow(t_am1, 2)+20*r[5]*pow(t_am1, 3);
	max_a1 = max_a1*ts*ts;	
	double max_a2 = 2*r[2]+6*r[3]*t_am2+12*r[4]*pow(t_am2, 2)+20*r[5]*pow(t_am2, 3);
	max_a2 = max_a2*ts*ts;
	cout << "calculated max acc: " << max_a1 << " " << max_a2 << endl;
	max_acc = fabs(max_a1);
	if(fabs(max_a2) > fabs(max_a1)) {
		max_acc = fabs(max_a2);
	}
	
	double t_jm = -0.2*r[4]/r[5]; 
	//cout << "max jerk at: " << t_jm << endl;
	max_jerk = 6*r[3]+24*r[4]*t_jm+60*r[5]*t_jm*t_jm;
	max_jerk = fabs(max_jerk*pow(ts, 3));
	cout << "calculated max jert: " << max_jerk << endl;
}

void dump_vector(const vector<double> &v)
{
	for(int i = 0; i < v.size(); ++i) {
		cout << v[i] << " ";
	}
	cout << endl;
}

void dump_vector(const vector<int> &v)
{
	for(int i = 0; i < v.size(); ++i) {
		cout << v[i] << " ";
	}
	cout << endl;
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int get_lane_num(double d, double &gd)
{
	int i = 0;
	
	gd = 2.0;
	if(d > 8.0) {
		gd = 10;
		i  = 2;
	}
	else if(d > 4.0) {
		gd = 6;
		i  = 1;
	}
	
	return i;
}

// we know s2 is ahead of s1
double get_s_distance(const double s1, const double s2)
{
	double d = s2 - s1;
	
	if(d < 0) {
		d = s2 + max_s - s1;
	}
	return d;
}

// get the car ahead within check_dist range, return -1 if no found, otherwise return the index 
// of the car in car vector
int get_closing_car_ahead(const vector<Car> & cars, double car_d, double car_s, double check_dist)
{
	int nc = cars.size();
	double td;
	
	int cid = -1;
	double min_ds = check_dist;
	
	//cout << "gcca: my car: " << car_s << " " << car_d << " " << check_dist << " " << cars.size() << endl;
	for(int i = 0; i < nc; ++i) {
		if(get_lane_num(car_d, td) != get_lane_num(cars[i].car_d, td)) continue;
		
		//cout << "gcca: checking car " << cars[i].car_id << " ";
		//cars[i].print();
		double d = cars[i].car_s - car_s;
		if(d < 0 && (car_s + check_dist) > max_s) {
			double ds = car_s + check_dist - max_s;
			if(cars[i].car_s < ds && (cars[i].car_s+max_s-car_s) < min_ds) {
				min_ds = cars[i].car_s + max_s - car_s;
				cid = i;
			}
			//cout << "gcca d<0: " << d << " " << ds << " " << min_ds << " " << cid << endl;
		} else if( d > 0 && d < min_ds ) {
			min_ds = d;
			cid = i;
			//cout << "gcca d > 0 " << d << " " << cid << " " << min_ds << endl;
		}
	}
	return cid;
}

// function return waypoints cross car_s for spline interpolation
vector<double> get_map_waypoints_inter(
	double car_s,
	const  vector<double>& map_s,
	double &local_car_s,
	vector<int>& map_idx,
	int num_prev,
	int num_next
)
{
	if(car_s > map_s[map_s.size()-1]) {
		car_s = map_s[map_s.size()-1];
	}

	auto it = lower_bound(map_s.begin(), map_s.end(), car_s);
	int idx = it - map_s.begin();
	
	for(int i = num_prev; i > 0; --i) {
		map_idx.push_back((map_s.size()+idx-i) % map_s.size());
	}
	
	for(int i = 0; i < num_next; ++i) {
		map_idx.push_back((idx+i) % map_s.size());
	}
	
	vector<double> local_map_waypoints_s;
	
	double s0 = map_s[map_idx[0]];
	local_map_waypoints_s.push_back(0.0);
	
	for(int i = 1; i < map_idx.size(); ++i) {
		double ds = map_s[map_idx[i]] - s0;
		
		if(ds < 0) {  // cross the end point, need special care
			ds = map_s[map_idx[i]] + max_s - s0;
		}
		local_map_waypoints_s.push_back(ds);
	}
	
	if(car_s < s0) {
		local_car_s = car_s + max_s - s0;
	} else {
		local_car_s = car_s - s0;
	}
	
	//cout << "map index to fit: ";
	//dump_vector(map_idx);
	//cout << "local s: " << local_car_s << endl;
	
	return local_map_waypoints_s;
}

void waypoints_spline_fit(
	double car_s,
	double& local_car_s,
	double& s0,
	tk::spline  &spline_fit_sx,
	tk::spline	&spline_fit_sy,
	tk::spline	&spline_fit_sdx,
	tk::spline	&spline_fit_sdy
)
{
	vector<int>	map_s_idx;
	
	auto map_fit_s = get_map_waypoints_inter(car_s, map_waypoints_s, local_car_s, map_s_idx, 10, 15);
	s0 = map_waypoints_s[map_s_idx[0]];

	vector<double> map_fit_x;
	vector<double> map_fit_y;
	vector<double> map_fit_dx;
	vector<double> map_fit_dy;
	
	for(int i = 0; i < map_s_idx.size(); ++i) {
		int idx = map_s_idx[i];
		map_fit_x.push_back(map_waypoints_x[idx]);
		map_fit_y.push_back(map_waypoints_y[idx]);
		map_fit_dx.push_back(map_waypoints_dx[idx]);
		map_fit_dy.push_back(map_waypoints_dy[idx]);
	}
	
	spline_fit_sx.set_points(map_fit_s, map_fit_x);
	spline_fit_sy.set_points(map_fit_s, map_fit_y);
	spline_fit_sdx.set_points(map_fit_s, map_fit_dx);
	spline_fit_sdy.set_points(map_fit_s, map_fit_dy);
}

vector<double> getXY_splines_fit(
	double s, 
	double d, 
	const tk::spline &spline_fit_sx, 
	const tk::spline &spline_fit_sy,
	const tk::spline &spline_fit_sdx,
	const tk::spline &spline_fit_sdy
)
{
	double x = spline_fit_sx(s);
	double y = spline_fit_sy(s);
	double dx = spline_fit_sdx(s);
	double dy = spline_fit_sdy(s);
	
	x += dx*d;
	y += dy*d;
	
	return {x, y};
}

vector<Car> get_cars_in_lane_segment(const vector<Car> &cars, double car_s, double car_d, double dist)
{
	double s1, s2;  // segment from s1 to s2
	vector<Car> cars_found;
	
	s1 = car_s -dist;
	if(s1 < 0) {
		s1 += max_s;
	}
	s2 = car_s +dist;
	if(s2 > max_s) {
		s2 -= max_s;
	}
	
	double td;
	for(int i = 0; i < cars.size(); ++i) {
		if(get_lane_num(car_d, td) != get_lane_num(cars[i].car_d, td)) continue;
			
		if(s2 > s1 && cars[i].car_s > s1 && cars[i].car_s < s2) {
			cars_found.push_back(cars[i]);
		}
		else if(s1 > s2 && (cars[i].car_s < s1 || cars[i].car_s > s2)) {
			cars_found.push_back(cars[i]);
		}
	}
	
	return cars_found;
}

// check if it is safe to change to left lane, car_s and car_d is
// current car's coord
bool safe_to_left_lane(const vector<Car> &cars, double car_d, double car_s, double car_speed)
{
	double td;
	
	if (0 == get_lane_num(car_d, td) )
		return false;  // can't go left
	
	// if there is no car in 10 meters ahead and following, saft to change lane
	auto left_cars = get_cars_in_lane_segment(cars, car_s, car_d-4.0, 30);
	cout << "left cars found: " << left_cars.size() << endl;
	if(left_cars.empty()) return true;
	
	/* vector<double> c_idd = get_closing_car_behind(cars, car_d-4.0, car_s, 20);
	if(c_idd[1] >= 10 && c_idd[1] < 20 && cars[int(c_idd[0])].get_car_vel() < car_speed) {
		return true;
	} */
	return false;
}

bool safe_to_right_lane(const vector<Car> &cars, double car_d, double car_s, double car_speed)
{
	double td;
	
	if (2 == get_lane_num(car_d, td) )
		return false;  // can't go left
	
	// if there is no car in 10 meters ahead and following, saft to change lane
	auto left_cars = get_cars_in_lane_segment(cars, car_s, car_d+4.0, 30);
	cout << "right cars found: " << left_cars.size() << endl;
	if(left_cars.empty()) return true;
	
	/* vector<double> c_idd = get_closing_car_behind(cars, car_d+4.0, car_s, 20);
	if(c_idd[1] >= 10 && c_idd[1] < 20 && cars[int(c_idd[0])].get_car_vel() < car_speed) {
		return true;
	} */
	return false;
}

// car_start is vector of {car_s, car_speed, car_a, car_d, car_d_speed, car_d_a}
int jmt_trajectory(
	const vector<double> &car_state,
	const vector<Car> &cars,
	const int closing_car,
	int num_steps,
	vector<double> &next_x,
	vector<double> &next_y,
	vector<vector<double> > &next_sva,
	int pidx
)
{
	double cost = 0.0;
	unsigned int next_st = CAR_STATE[0];
	
	double local_car_s, s0;
	vector<double> map_fit_waypoints_s;
	vector<double> map_fit_waypoints_x;
	vector<double> map_fit_waypoints_y;
	
	tk::spline  spline_fit_sx;
	tk::spline	spline_fit_sy;
	tk::spline	spline_fit_sdx;
	tk::spline	spline_fit_sdy;
	
	cout << "start state for trajectory calculation: " << endl;
	dump_vector(car_state);
	
	waypoints_spline_fit(car_state[0], local_car_s, s0, spline_fit_sx, spline_fit_sy, spline_fit_sdx, spline_fit_sdy);

	cout <<"local car_s is: " << local_car_s << endl;

	vector<double> s_start = {local_car_s, car_state[1], car_state[2]};
	vector<double> d_start = {car_state[3], car_state[4], car_state[5]};
	double max_dist_per_step = max_speed; //0.445;  // each step is 20 ms, about 20 m/s which is ~50MPH
	
	double goal_s = s_start[0] + max_dist_per_step*num_steps;
	
	if(car_state[1] < max_dist_per_step/2.5) {
		max_dist_per_step = 0.4;
		goal_s = s_start[0] + max_dist_per_step*num_steps/2.0;
		max_dist_per_step /= 2.0;
	}
	vector<double> s_goal = {goal_s, max_dist_per_step, 0};
	
	double goal_d = 2.0;
	if(car_state[3] > 8.0) {
		goal_d = 9.85;
	}
	else if(car_state[3] > 4.0) {
		goal_d = 6;
	}
	
	vector<double> d_goal = {goal_d, 0, 0};
	
	// TODO: based on cost function
	//int closing_car = get_closing_car_ahead(cars, car_state[3], car_state[0], 25);
	//cout << "jmt trajectory generation find closing car: " << closing_car << endl;
	
	if(CAR_STATE[0] == _CS_LCL) {
		double td; 
		if(get_lane_num(car_state[3], td) == CAR_STATE[1]) {
			cout << "change LCL to KL state" << endl;
			next_st = _CS_KL;
			//CAR_STATE[0] = _CS_KL;
		} else {
			cout << "keeping in lane change state!" << endl;
			d_goal[0] = CAR_STATE[1]*4.0 + 2.0;
		}
	} else if(CAR_STATE[0] == _CS_LCR) {
		double td; 
		if(get_lane_num(car_state[3], td) == CAR_STATE[1]) {
			cout << "change LCR to KL state" << endl;
			//CAR_STATE[0] = _CS_KL;
			next_st = _CS_KL;
		} else {
			cout << "keeping in right lane change state!" << endl;
			d_goal[0] = CAR_STATE[1]*4.0 + 2.0;
		}
	}
    else if(closing_car >= 0 && cars[closing_car].get_car_vel() < max_speed) {
   		double td;
    	if(safe_to_left_lane(cars, car_state[3], car_state[0], car_state[2])) {
    		cout << "now changing to left lane!" << endl;
    		d_goal[0] -= 4.0;
    		/* if(CAR_STATE[0] == _CS_KL_FL) {
    			s_goal[1] = cars[closing_car].get_car_vel();
    			s_goal[0] = s_start[0] + cars[closing_car].get_car_vel()*num_steps;
    		} */
    		//if(d_goal[0] == 2.0) {
    		//	d_goal[0] -= 0.1;
    		//}
    		//CAR_STATE[0] = _CS_LCL;
    		next_st = _CS_LCL;
    		CAR_STATE[1] = get_lane_num(d_goal[0], td);
    	}
    	else if(safe_to_right_lane(cars, car_state[3], car_state[0], car_state[2])) {
    		cout << "now changing to right lane!" << endl;
    		d_goal[0] += 4.0;
    		/* if(CAR_STATE[0] == _CS_KL_FL) {
    			s_goal[1] = cars[closing_car].get_car_vel();
    			s_goal[0] = s_start[0] + cars[closing_car].get_car_vel()*num_steps;
    		} */
    		//if(d_goal[0] == 10.0) {
    		//	d_goal[0] -= 0.1;
    		//}
    		//CAR_STATE[0] = _CS_LCR;
    		next_st = _CS_LCR;
    		CAR_STATE[1] = get_lane_num(d_goal[0], td);
    	} else {
     		cout << "following car ahead: " << closing_car << endl;
    		cars[closing_car].print();
    		next_st = _CS_KL_FL;
    		//klfl_get_best_jmt(car_state[0], s_start, cars[closing_car]);
    		//s_goal[0] = s_start[0] + cars[closing_car].get_car_vel()*num_steps;
    		double ds = get_s_distance(car_state[0], cars[closing_car].car_s);
    		//num_steps = 150;
    		s_goal[0] = s_start[0] + ds + cars[closing_car].get_car_vel()*num_steps - 30;
    		s_goal[1] = cars[closing_car].get_car_vel();
    		//pidx = -1;
    		//next_sva.clear();    		
		}   	
    }
    else {
    	cout << "nothing changed, keep going!" << endl;
    	next_st = _CS_KL;
    }
	// end of TODO
	
	cout << "s start: " << s_start[0] << " " << s_start[1] << " " << s_start[2] << " " << endl;
	cout << "s goal: " << s_goal[0] << " " << s_goal[1] << " " << s_goal[2] << " " << endl;	
	cout << "d start: " << d_start[0] << " " << d_start[1] << " " << d_start[2] << " " << endl;
	cout << "d goal: " << d_goal[0] << " " << d_goal[1] << " " << d_goal[2] << " " << endl;	
	
	auto jmt_s_params = JMT(s_start, s_goal, num_steps);
	auto jmt_d_params = JMT(d_start, d_goal, num_steps);

	double max_acc, max_jerk;
	max_jmt_acc_jerk(jmt_s_params, 0.02, max_acc, max_jerk);
		
	CAR_STATE[0] = next_st;
	
	double previous_s = 0;
	int num_reuse = 0;
	if(pidx > 0) num_reuse = 15;
	double px = 0.0;
	double py = 0.0;
	double pv = 0.0;
	for(int t = 0; t < num_steps/2; ++t) {	
		double s   = quintic_eval_s(jmt_s_params, t);
		double s_v = quintic_eval_d(jmt_s_params, t);
		double s_a = quintic_eval_a(jmt_s_params, t);
		
		if(t < num_reuse) {  // reuse previous path
			cout << "spline s is: " << s << " " << next_sva[pidx+t][0] << " " << s0;
			s   = next_sva[pidx+t][0] - s0;
			if(s < 0) s += max_s;
			cout << " reuse s is: " << s << endl;
			s_v = next_sva[pidx+t][1];
			s_a = next_sva[pidx+t][2];
		}
		
		double d   = quintic_eval_s(jmt_d_params, t);
		double d_v = quintic_eval_d(jmt_d_params, t);
		double d_a = quintic_eval_a(jmt_d_params, t);
		//cout << d << " ";
		
		if(t > 0 && (s - previous_s) > max_speed) {
			s = previous_s + max_speed;
		}
		
		if(t > 1) {
			double ds = s - previous_s;
			double da = ds - pv;
			double acc = fabs(da) * 50 * 50;
			if(acc > 9.5) {
				cout << "acc exceeded!!!!! " << s << " " << previous_s << " " << acc << endl;
			
				while(acc > 9.5) {
					s = s - 0.001*(fabs(da)/da);
					ds = s - previous_s;
					da = ds - pv;
					acc = fabs(da) * 50 * 50;
				}
				
				cout << "reduce acc to: " << acc << " s=" << s << endl;
			}
		}
		
		auto xy = getXY_splines_fit(s, d, spline_fit_sx, spline_fit_sy, spline_fit_sdx, spline_fit_sdy);
		if(t > 0) { // check result speed
			double dxy = distance(px, py, xy[0], xy[1]);
			if(dxy > max_speed) {
				cout << "fitting dist too high!, reduce s.." << dxy << endl;
				while(dxy > max_speed) {
					s -= 0.01;
					xy = getXY_splines_fit(s, d, spline_fit_sx, spline_fit_sy, spline_fit_sdx, spline_fit_sdy);
					dxy = distance(px, py, xy[0], xy[1]);
				}
			}
		}
		
		if( t > 0 ) pv = s - previous_s;
		
		previous_s = s;
		
		next_x.push_back(xy[0]);
		next_y.push_back(xy[1]);
		px = xy[0];
		py = xy[1];
		
		s += s0;
		if(s > max_s) s -= max_s;
		if(pidx < 0) {
			next_sva.push_back({s, s_v, s_a, d, d_v, d_a, xy[0], xy[1]});
		} else {
			//cout << "t is " << t << " " << next_sva.size() << endl;
			next_sva[t] = {s, s_v, s_a, d, d_v, d_a, xy[0], xy[1]};
		}
	}
	cout << "complete trajectory calculation!" << endl;
	
	return cost;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  /* vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy; */
  vector<vector<double> > prev_sdva;
  int num_traj_steps = 200;
  vector<Car> cars;


  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  //h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
  h.onMessage([&prev_sdva, &num_traj_steps, &cars](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
          	
          	cars.clear();
          	//cout << "cars on road: " << endl;
          	for(int i = 0; i < sensor_fusion.size(); ++i) {
          		Car car(sensor_fusion[i][0], sensor_fusion[i][1], sensor_fusion[i][2], sensor_fusion[i][3], sensor_fusion[i][4],
          			    sensor_fusion[i][5], sensor_fusion[i][6]);
          		cars.push_back(car);
          		//car.print();
          	}
          	//cout << endl;
          	int closing_car = get_closing_car_ahead(cars, car_d, car_s, 30);
          	if(closing_car >= 0) {
          		cout << "closing car: ";
          		cars[closing_car].print();
          	}

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
          	
			cout << "car s=" << car_s << " d=" << car_d << " x=" << car_x << " y=" << car_y << " speed=" << car_speed  << " yaw=" << car_yaw << endl;
          	vector<double> car_st = {car_x, car_y, car_s, car_d, car_yaw, car_speed};

			if(prev_sdva.size() == 0) { // car just start!
				double td;
				CAR_STATE = { _CS_KL, get_lane_num(car_d, td) };
				vector<double> car_state = {car_s, 0, 0, car_d, 0, 0};

				jmt_trajectory(car_state, cars, closing_car, num_traj_steps, next_x_vals, next_y_vals, prev_sdva, -1);
			} else {	
				int idx = prev_sdva.size() - previous_path_x.size();
				
				// here estimate car s and d verlocity and acceleration
				double s_v0 = prev_sdva[idx][0];
				double s_v1 = prev_sdva[idx+1][0];
				double s_v2 = prev_sdva[idx+2][0];
				double d_v0 = prev_sdva[idx][3];
				double d_v1 = prev_sdva[idx+1][3];
				double d_v2 = prev_sdva[idx+3][3];
				
				double est_sv = s_v1 - s_v0;
				double est_sa = (s_v2 - s_v1) - (s_v1 - s_v0);
				double est_dv = d_v1 - d_v0;
				double est_da = (d_v2 - d_v1) - (d_v1 - d_v0);
				cout << "idx is: " << idx << endl;
				if(idx > 50 || (closing_car >= 0 && CAR_STATE[0] == _CS_KL)) {																	
					vector<double> car_state = {prev_sdva[idx][0], prev_sdva[idx][1], prev_sdva[idx][2], prev_sdva[idx][3], prev_sdva[idx][4], prev_sdva[idx][5]};
					
					//prev_sdva.clear();
					jmt_trajectory(car_state, cars, closing_car, num_traj_steps, next_x_vals, next_y_vals, prev_sdva, idx);							
				} else {
					cout << "continue the trajectory" << endl;
					for(int i = 0; i < previous_path_x.size(); i++) {
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}
				}
          	}
          	
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































