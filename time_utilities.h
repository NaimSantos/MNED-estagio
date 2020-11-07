#ifndef CUSTOMTIMER_H
#define CUSTOMTIMER_H

#include <chrono>
#include <ctime>

class CustomTimer
{
	public:
		CustomTimer(){
			startpoint = std::chrono::steady_clock::now();
		}
		~CustomTimer(){
			TimeDestructor();
		}

		void TimeDestructor(){
			auto endpoint = std::chrono::steady_clock::now();

			auto start = std::chrono::time_point_cast<std::chrono::miliseconds>(startpoint).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<std::chrono::miliseconds>(endpoint).time_since_epoch().count();

			auto total = end - start;
			double s = total*0.001;
			
			std::cout << "\nTempo decorrido:" << s << " segundos (" << total << " milisegundos)" << std::endl;
		}
	
		private:
			std::chrono::time_point<std::chrono::steady_clock> startpoint;
};

std::chrono::time_point<std::chrono::steady_clock> capturetime(){
	auto res = std::chrono::steady_clock::now();
	return res;
}

double get_elapsed_time(std::chrono::time_point<std::chrono::steady_clock> ti, std::chrono::time_point<std::chrono::steady_clock> tf){
	double elapsed_time_ms = std::chrono::duration_cast <std::chrono::miliseconds> (tf - ti).count();
	return elapsed_time_ms;
}

void print_current_time(){
	std::time_t now = std::time(nullptr);

	std::cout << std::asctime(std::localtime(&now));
}

#endif //CUSTOMTIMER_H
