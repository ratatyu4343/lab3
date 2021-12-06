#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <cmath>
#include <ga/ga.h> //GALib

//fuctions that need 2 arguments
const int FUNCTIONS_BINAR_SIZE = 3;
const std::string FUNCTIONS_BINAR[] = { "*", "+", "-" };

//fuctions that need 1 arguments
const int FUNCTIONS_UNAR_SIZE = 3;
const std::string FUNCTIONS_UNAR[] = { "sin", "cos", "abs"};
	
//max depth of notation
const int MAX_DEPTH = 2;

std::vector<double> X; // X
std::vector<double> Y; // F(X)
	
//creat genome with random polish notation
void initialization(GAGenome& g);
	
//calculate sum of delta^2
//where delta = x_real - x_that_genome_returns
//function return 1/sum, (maximum when sum->0)
float fitnes(GAGenome& g);
	
//does random activity
//like change function/var
//can build new function (p = 0.1)
int mutation(GAGenome& g, float pmut);
	
//make 2 random child by 2 parents 
int mycrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* g3, GAGenome* g4);

//return random polish nonation
std::vector<std::string> random_vector();
	
//cout << genome
void print(GA1DArrayGenome<std::string>& g);
	
//return type of object string s contains
//0 - x, 1 - number, 2 - binar function, 3 - unar function
int type(std::string s);
	
//return function(x)
double calculate(GA1DArrayGenome<std::string>& g, double x);

//coefitients algotythm
GA1DArrayGenome<std::string>* GLOBAL_GENOME = NULL;
void initial(GAGenome& g)
{
	GA1DArrayGenome<std::string>& gg = (GA1DArrayGenome<std::string>&)g;
	gg.resize(GLOBAL_GENOME->size());

	for (int i = 0; i < GLOBAL_GENOME->size(); i++)
	{
		if (type(GLOBAL_GENOME->gene(i)) == 1)
		{
			gg.gene(i, std::to_string(GARandomDouble(-10, 10)));
		}
		else
		{
			gg.gene(i, GLOBAL_GENOME->gene(i));
		}
	}

	g = gg;
}

int mutt(GAGenome& g, float pmut)
{
	if (pmut >= GARandomFloat(0, 1))
	{
		GA1DArrayGenome<std::string>& gg = (GA1DArrayGenome<std::string>&)g;
		for (int i = 0; i < gg.size(); i++)
		{
			if (type(gg.gene(i)) == 1)
			{
				int coin = GARandomInt(0, 1);
				if(coin)
					gg.gene(i, std::to_string(GARandomDouble(-10, 10)));
			}
		}
		g = gg;
		return 1;
	}
	
	return 0;
}

int mycros(const GAGenome& g1, const GAGenome& g2, GAGenome* g3, GAGenome* g4)
{
	int n = 0;
	GA1DArrayGenome<std::string>& gen1 = (GA1DArrayGenome<std::string> &)(g1);
	GA1DArrayGenome<std::string>& gen2 = (GA1DArrayGenome<std::string> &)(g2);
	if (g3)
	{
		GA1DArrayGenome<std::string>& gen3 = (GA1DArrayGenome<std::string> &)(*g3);
		gen3.resize(gen1.size());

		for (int i = 0; i < gen1.size(); i++)
		{
			int coin = GARandomInt(0, 1);
			if (coin)
			{
				gen3.gene(i, gen1.gene(i));
			}
			else
			{
				gen3.gene(i, gen2.gene(i));
			}
		}
		*g3 = gen3;
		n++;
	}
	if (g4)
	{
		GA1DArrayGenome<std::string>& gen4 = (GA1DArrayGenome<std::string> &)(*g4);
		gen4.resize(gen1.size());

		for (int i = 0; i < gen1.size(); i++)
		{
			int coin = GARandomInt(0, 1);
			if (coin)
				gen4.gene(i, gen1.gene(i));
			else
				gen4.gene(i, gen2.gene(i));
		}
		*g4 = gen4;
		n++;
	}

	return n;
}

std::vector<std::string> best_coef(GA1DArrayGenome<std::string>& g)
{
	GLOBAL_GENOME = new GA1DArrayGenome<std::string>(0, fitnes);
	GLOBAL_GENOME->resize(g.size());
	for (int i = 0; i < g.size(); i++)
	{
		GLOBAL_GENOME->gene(i, g.gene(i));
	}
	GA1DArrayGenome<std::string> new_gen(0, fitnes);

	new_gen.initializer(initial);
	new_gen.mutator(mutt);
	new_gen.crossover(mycros);
	GASimpleGA al (new_gen);
	
	al.populationSize(5);
	al.nGenerations(15);
	al.pMutation(0.3);

	al.evolve();

	delete GLOBAL_GENOME;
	
	std::vector<std::string> v;
	for (int i = 0; i < ((GA1DArrayGenome<std::string>&)al.statistics().bestIndividual()).size(); i++)
		v.push_back(((GA1DArrayGenome<std::string>&)al.statistics().bestIndividual()).gene(i));
	return v;
}

std::vector<std::string> symbolregression(std::vector<double> x,std::vector<double> y, int population, int count, double pmut)
{
	std::vector<std::string> v;

	X = x;
	Y = y;

	GA1DArrayGenome<std::string> g(0, fitnes);
	g.initializer(initialization);
	g.mutator(mutation);
	g.crossover(mycrossover);
	
	GASimpleGA alg(g);
	alg.populationSize(population);
	alg.nGenerations(count);
	alg.pMutation(pmut);
	alg.scoreFilename("information.txt");
	alg.selectScores(GAStatistics::Maximum);
	alg.flushFrequency(1);
	alg.evolve();

	GA1DArrayGenome<std::string>& best = (GA1DArrayGenome<std::string>&)alg.statistics().bestIndividual();
	for (int i = 0; i < best.size(); i++)
	{
		v.push_back(best.gene(i));
	}

	return v;
}

std::vector<std::string> random_vector()
{
	std::vector<std::string> v;
	int n = 1;
	while (n > 0)
	{
		if (n < MAX_DEPTH)
		{
			int r = GARandomInt(0, 2);
			if (r == 0)
			{
				int rr = GARandomInt(0, 1);
				if (rr == 0)
				{
					v.push_back(std::to_string(GARandomDouble(0, 10)));
				}
				else
				{
					v.push_back("x");
				}
				n--;
			}
			else
				if (r == 1)
				{
					v.push_back(FUNCTIONS_BINAR[GARandomInt(0, FUNCTIONS_BINAR_SIZE - 1)]);
					n++;
				}
				else
					if (r == 2)
					{
						v.push_back(FUNCTIONS_UNAR[GARandomInt(0, FUNCTIONS_UNAR_SIZE - 1)]);
					}
		}
		else
		{
			int rr = GARandomInt(0, 1);
			if (rr == 0)
			{
				v.push_back(std::to_string(GARandomDouble(0, 10)));
			}
			else
			{
				v.push_back("x");
			}
			n--;
		}
	}
	return v;
}

void initialization(GAGenome& g)
{
	std::vector<std::string> v = random_vector();

	auto s = (GA1DArrayGenome<std::string> &) g;
	s.resize(v.size());

	for (int i = 0; i < v.size(); i++)
	{
		s.gene(i, v[i]);
	}

	g = s;
}

void print(GA1DArrayGenome<std::string>& g)
{
	std::cout << "->";
	for (int i = 0; i < g.size(); i++)
	{
		std::cout << g.gene(i) << "->";
	}
	std::cout << "\n";
}

int type(std::string s)
{
	if (s == "x")
		return 0;
	for (int i = 0; i < FUNCTIONS_BINAR_SIZE; i++)
		if (s == FUNCTIONS_BINAR[i])
		{
			return 2;
		}
	for (int i = 0; i < FUNCTIONS_UNAR_SIZE; i++)
		if (s == FUNCTIONS_UNAR[i])
		{
			return 3;
		}
	return 1;
}

double calculate(GA1DArrayGenome<std::string>& g, double x)
{
	int len = g.size();

	std::stack<double> s;
	for (int i = len - 1; i >= 0; i--)
	{
		int typ = type(g.gene(i));
		if (typ == 0)
		{
			s.push(x);
		}
		else
		if (typ == 1)
		{
			s.push(atof(g.gene(i).c_str()));
		}
		else
		if (typ == 2)
		{
			double a = s.top();
			s.pop();
			double b = s.top();
			s.pop();
			double c = 0;
			if (g.gene(i) == "+")
			{
				c = a + b;
			}
			else
			if (g.gene(i) == "-")
			{
				c = a - b;
			}
			else
			if (g.gene(i) == "*")
			{
				c = a * b;
			}
			s.push(c);
		}
		else
		{
			double a = s.top();
			s.pop();
			double c = 0;
			if (g.gene(i) == "sin")
			{
				c = sin(a);
			}
			else
			if (g.gene(i) == "cos")
			{
				c = cos(a);
			}
			else
			if (g.gene(i) == "abs")
			{
				c = abs(a);
			}
			else
			if (g.gene(i) == "sgn")
			{
				if (a > 0)
					c = 1;
				else 
				if (a < 0)
					c = -1;
				else
					c = 0;
			}
			s.push(c);
		}
	}
	
	return s.top();
}

int mutation(GAGenome& g, float pmut)
{
	GA1DArrayGenome<std::string>& s = (GA1DArrayGenome<std::string> &)(g);

	double p = GARandomDouble(0, 1);
	if (pmut >= p)
	{
		int r = GARandomInt(0, s.size() - 1);
		if (type(s.gene(r)) == 0)
		{
			int coin = GARandomInt(0, 1);
			if(coin == 1)
				s.gene(r, std::to_string(GARandomDouble(0, 10)));
		}
		else
		if (type(s.gene(r)) == 1)
		{
			for (int i = 0; i < s.size() - 1; i++)
			{
				if (type(s.gene(i)) == 1)
				{
					s.gene(i, std::to_string(GARandomDouble(-10, 10)));
				}
			}
		}
		else
		if (type(s.gene(r)) == 2)
		{
			s.gene(r, FUNCTIONS_BINAR[GARandomInt(0, FUNCTIONS_BINAR_SIZE - 1)]);
		}
		else
		if (type(s.gene(r)) == 3)
		{
			s.gene(r, FUNCTIONS_UNAR[GARandomInt(0, FUNCTIONS_UNAR_SIZE - 1)]);
		}
		
		int n = GARandomInt(0, 15);
		if (n == 1 || n == 11)
		{
			std::vector<std::string> v = random_vector();
			s.resize(v.size());
			for (int i = 0; i < v.size(); i++)
			{
				s.gene(i, v[i]);
			}
		}
		else
		if (n == 2)
		{
			int t = type(s.gene(0));
			if (t == 3)
			{
				for (int i = 0; i < s.size() - 1; i++)
				{
					s.gene(i, s.gene(i + 1));
				}
				s.resize(s.size() - 1);
			}
		}
		else
		if (n == 3)
		{
			s.resize(s.size() + 1);
			for (int i = s.size() - 1; i > 0; i--)
			{
				s.gene(i, s.gene(i - 1));
			}
			s.gene(0, FUNCTIONS_UNAR[GARandomInt(0, FUNCTIONS_UNAR_SIZE - 1)]);
		}
		else
		if (n == 4)
		{
			s.resize(s.size() + 2);
			for (int i = s.size() - 3; i >= 0; i--)
			{
				s.gene(i+2, s.gene(i));
			}
			s.gene(0, FUNCTIONS_BINAR[GARandomInt(0, FUNCTIONS_BINAR_SIZE - 1)]);
			s.gene(1, std::to_string(GARandomDouble(0, 10)));
		}
		else
		if (n == 5)
		{
			if ((s.size() > 1) && (type(s.gene(0)) == 2) && (type(s.gene(1)) <= 1))
			{
				for (int i = 0; i < s.size() - 2; i++)
				{
					s.gene(i, s.gene(i + 2));
				}
				s.resize(s.size() - 2);
			}
		}
		else
		if (n == 6)
		{
			for(int i = 0; i < s.size(); i++)
				if (type(s.gene(i)) == 1)
				{
					std::vector<std::string> v = best_coef(s);
					s.resize(v.size());
					for (int i = 0; i < s.size(); i++)
					{
						s.gene(i, v[i]);
					}
					break;
				}
		}

		g = s;
		return 1;
	}
	return 0;
}

int mycrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* g3, GAGenome* g4)
{
	int n = 0;

	GA1DArrayGenome<std::string>& gen1 = (GA1DArrayGenome<std::string> &)(g1);
	GA1DArrayGenome<std::string>& gen2 = (GA1DArrayGenome<std::string> &)(g2);

	int len = gen1.size() < gen2.size() ? gen1.size() : gen2.size();
	int max_len = gen1.size() > gen2.size() ? gen1.size() : gen2.size();
	int coin;

	if (g3)
	{
		GA1DArrayGenome<std::string>& gen3 = (GA1DArrayGenome<std::string> &)(*g3);
		gen3.resize(gen1.size());

		gen3.copy(gen1, 0, 0, gen1.size());
		bool flag = false;
		for (int i = 0; i < len; i++)
		{
			coin = GARandomInt(0, 2);
			if (type(gen1.gene(i)) == type(gen2.gene(i)))
			{
				int coin2 = GARandomInt(0, 1);
				if (coin2 == 0)
				{
					gen3.gene(i, gen2.gene(i));
				}
			}
		}
		*g3 = gen3;
		n++;
	}
	if (g4)
	{
		GA1DArrayGenome<std::string>& gen4 = (GA1DArrayGenome<std::string> &)(*g4);
		gen4.resize(gen2.size());

		gen4.copy(gen2, 0, 0, gen2.size());
		for (int i = 0; i < len; i++)
		{
			coin = GARandomInt(0, 2);
			if (type(gen1.gene(i)) == type(gen2.gene(i)))
			{
				int coin2 = GARandomInt(0, 1);
				if (coin2 == 0)
				{
					gen4.gene(i, gen1.gene(i));
				}
			}
		}
		n++;
		*g4 = gen4;
	}
	return n;
}

float fitnes(GAGenome& g)
{
	GA1DArrayGenome<std::string>& gen = (GA1DArrayGenome<std::string>&)g;

	double sum = 0;
	for (int i = 0; i < X.size(); i++)
	{
		sum += pow(Y[i] - calculate(gen, X[i]), 2);
	}

	return 1 / sum;
}


