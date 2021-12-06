#include <iostream>
#include <fstream>
#include "symbol_regression.h"
#include <cmath>

//reads points from file
void read_points(std::vector<double>*, std::vector<double>*);
//change notation from polish to normal
std::string polish_to_normal(std::vector<std::string>);
//creat points and write to file
void creat_points(int start, int endd);

int main()
{
	//creat_points(-5, 5);
	std::vector<double> X;
	std::vector<double> Y;
	read_points(&X, &Y);

	while (true)
	{
		int n;
		std::cout << "1 - Regression\n2 - exit()\n>>>";
		std::cin >> n;
		if (n == 1)
		{
			int gsize, psize;
			float pmut;
			std::cout << "Input population size: ";
			std::cin >> psize;
			std::cout << "Input generation size: ";
			std::cin >> gsize;
			std::cout << "Input mutation size (0, 1): ";
			std::cin >> pmut;

			std::vector<std::string> v = symbolregression(X, Y, psize, gsize, pmut);
			std::cout <<"Function: "<< polish_to_normal(v) << std::endl;
			
			std::ifstream summ("information.txt");
			std::string s_summ, sss;
			while (summ >> sss) s_summ = sss;
			std::cout << "Sum: " << s_summ << std::endl;
			summ.close();
		}
		else
		if (n == 2) return 0;	
	}
	
}

void read_points(std::vector<double>* X, std::vector<double>* Y)
{
	std::ifstream points("points.txt");
	std::string x_point, y_point;
	while (points >> x_point && points >> y_point)
	{
		X->push_back(atof(x_point.c_str()));
		Y->push_back(atof(y_point.c_str()));
	}
	points.close();
}

std::string polish_to_normal(std::vector<std::string> polish)
{
	std::stack<std::string> s;

	for (int i = polish.size() - 1; i >= 0; i--)
	{
		int t = type(polish[i]);
		if (t == 0 || t == 1)
		{
			s.push(polish[i]);
		}
		else
		if (t == 3)
		{
			std::string op1 = s.top();
			s.pop();
			s.push(polish[i] + "(" + op1 + ")");
		}
		else
		{
			std::string op1 = s.top();
			s.pop();
			std::string op2 = s.top();
			s.pop();
			s.push("(" + op1 + polish[i] + op2 + ")");
		}
	}
	return s.top();
}

void creat_points(int start, int endd)
{
	std::ofstream rrr("points.txt");
	
	for (int i = start; i < endd; i++)
	{
		rrr << i << " " << i*i + 3*i + 5 << "\n";
	}
	
	rrr.close();
}