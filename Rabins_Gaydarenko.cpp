// Rabins_Gaydarenko.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


#include "Rabins_Gaydarenko.h"



int main(int argc, char** argv)
{
	
	if (*argv[3] == 'f')
	{
		RabinSecret MyRabin(atoi(argv[1]), atoi(argv[2]), argv[4], argv[5]);
		MyRabin.RabinSharing();

	}
	else if (*argv[3] == 'b')
	{
		RabinSecret MyRabin(atoi(argv[1]), atoi(argv[2]), argv[4], argv[5], 1);
		MyRabin.RabinRecovery();

	}
	else std::cout << "use correct list of parameter";

	return 0;
}


