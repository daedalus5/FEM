#include <Partio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Particles.h"	
#include "Spring.h"
#include "TestCase.h"

using T = float;
const int dim = 3;

int main(int argc, char* argv[])
{   
	std::string s;
	std::ifstream infile;
	std::ofstream outfile;
	infile.open("actual_input.txt");
	outfile.open("actual_output.txt");

	if(infile.is_open()){
		while(!infile.eof()){
			getline(infile, s);
			TestCase<T,dim> test;
			test.initialize(s);
			if(outfile.is_open()){
				test.print_info(outfile);
			}
		}
	}
	infile.close();
	outfile.close();
    return 0;
}