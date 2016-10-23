#include <string>
#include <ctime>
using namespace std;

struct Time{
	time_t firsttime;
	time_t secondtime;
	string cas_string(time_t casdoc);
	time_t taketime();
	double delta(time_t t1a, time_t t2a);
	void CPUtime();
	string cas_ASCII();
};

