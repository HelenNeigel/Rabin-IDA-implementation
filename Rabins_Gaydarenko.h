#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "stdafx.h"
#include <direct.h>
#include <stdlib.h>
#include <fstream>
#include <ostream> 
#include <vector>
#include <random>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"

#include "Galois.h"

struct coord
{
	int i = -1;
	int j = -1;
};

class RabinSecret{

	unsigned char inverses[256];

	const int n_parties;
	const int k_threshold;

	char* secret_file_name;
	std::vector<char*> parties_files_names;

	std::vector<unsigned char> SecretBuffer;
	std::vector<std::vector<unsigned char>> PartiesBuffers;

	std::vector<std::vector<unsigned char>> MatrixCoefficients;
	std::vector<std::vector<unsigned char>> MinorMatrix;

	void Sharing();
	void Recovery();
	
	void MatrixRndGenerator();
	void GetRecoveryMatrix();
	unsigned char FindDenom(std::vector<coord>);
	
public:
	RabinSecret(int parties, int threshold, char* secret, char* shares, int set = 0) :n_parties(parties), k_threshold(threshold)
	{

		secret_file_name = new char[strlen(secret) + 1];
		strcpy(secret_file_name, secret);

		boost::filesystem::path dir_temp(boost::filesystem::current_path());

		std::string temp_s(dir_temp.string());
		temp_s.append("\\");
		temp_s.append(shares);
		boost::filesystem::path dir(temp_s);

		if (set == 1){

			boost::filesystem::directory_iterator end_iter;
			int i = 0;
			for (boost::filesystem::directory_iterator dir_itr(dir); dir_itr != end_iter; ++dir_itr)
			{
				std::string temp_str(temp_s);
				temp_str.append("\\");
				temp_str.append(dir_itr->path().filename().string());


				char* temp = new char[temp_str.size()];
				strcpy(temp, temp_str.c_str());

				parties_files_names.push_back(temp);

				i++;
			}

		}
		else
		{
			char buffer[3];
	
			boost::filesystem::path dir(boost::filesystem::current_path());

			for (int i = 0; i < n_parties; ++i)
			{
				std::string temp_str(dir.string());
				temp_str.append("\\");
				temp_str.append(shares);
				temp_str.append("\\");
				temp_str.append(secret);
				int pos = temp_str.find('.');
				temp_str.insert(pos, itoa(i + 1, buffer, 10));
				char* temp = new char[temp_str.size()];
				strcpy(temp, temp_str.c_str());

				parties_files_names.push_back(temp);
				
			}

		}
	}

	void RabinSharing();
	void RabinRecovery();

};



void RabinSecret::MatrixRndGenerator()
{

	std::random_device rd;
	for (int i = 0; i < n_parties; ++i)
	{
		MatrixCoefficients.resize(i + 1, std::vector<unsigned char>());
		MatrixCoefficients[i].push_back(1);
		for (int j = 1; j < k_threshold; ++j)
			MatrixCoefficients[i].push_back(rd() % 256);
	}

}

void RabinSecret::GetRecoveryMatrix()
{
	MinorMatrix.resize(k_threshold, std::vector<unsigned char>());
	std::vector<coord> mem_ignore;
	unsigned char denom = FindDenom(mem_ignore);

	//Fill MinorMatrix
	for (int i = 0; i < MatrixCoefficients.size(); ++i)
	for (int j = 0; j < MatrixCoefficients.size(); ++j)
	{
		coord addition;
		addition.i = i;
		addition.j = j;
		std::vector<coord> mem_ignore_t;
		mem_ignore_t.push_back(addition);

		MinorMatrix[i].push_back(FindDenom(mem_ignore_t));
	}
	// Transpose MinorMatrix
	for (int i = 0; i < MinorMatrix.size(); ++i)
	for (int j = i + 1; j < MinorMatrix.size(); ++j)
	{
		unsigned char temp = MinorMatrix[i][j];
		MinorMatrix[i][j] = MinorMatrix[j][i];
		MinorMatrix[j][i] = temp;
	}

	// Get Inverse Matrix as Transpose*(1/denom)
	for (int i = 0; i < MinorMatrix.size(); ++i)
	for (int j = 0; j < MinorMatrix.size(); ++j)
	{
		MinorMatrix[i][j] = GFMul(MinorMatrix[i][j], inverses[denom]);
	}

}

unsigned char RabinSecret::FindDenom(std::vector<coord> mem_ignore)
{
	unsigned char result = 0;

	if (mem_ignore.size() == k_threshold - 2)
	{
		unsigned char temp[4];
		int t = 0;

		for (int row = 0; row < MatrixCoefficients.size(); ++row)
		{
			int is_row_there = 0;

			for (int i = 0; i < mem_ignore.size(); ++i)
			if (mem_ignore[i].i == row) is_row_there = 1;

			if (!is_row_there)
			{
				for (int col = 0; col < MatrixCoefficients[0].size(); ++col)
				{
					int is_col_there = 0;

					for (int i = 0; i < mem_ignore.size(); ++i)
					if (mem_ignore[i].j == col) is_col_there = 1;

					if (!is_col_there)
					{
						temp[t] = MatrixCoefficients[row][col];
						++t;
					}
				}
			}
		}
		result = GFAdd(GFMul(temp[0], temp[3]), GFMul(temp[1], temp[2]));
	}

	else
	{
		unsigned char temp = 0;

		for (int row = 0; row < MatrixCoefficients.size(); ++row)
		{
			int is_row_there = 0;

			for (int i = 0; i < mem_ignore.size(); ++i)
			if (mem_ignore[i].i == row) is_row_there = 1;

			if (!is_row_there)
			{
				for (int col = 0; col < MatrixCoefficients[0].size(); ++col)
				{
					int is_col_there = 0;

					for (int i = 0; i < mem_ignore.size(); ++i)
					if (mem_ignore[i].j == col) is_col_there = 1;

					if (!is_col_there)
					{
						coord addition;
						addition.i = row;
						addition.j = col;


						mem_ignore.push_back(addition);

						temp = GFAdd(temp, GFMul(MatrixCoefficients[row][col], FindDenom(mem_ignore)));
						mem_ignore.pop_back();
					}
					
				}
					
				
				row = MatrixCoefficients.size();
				result = temp;
			}
		}
	}
 	return result;
}


void RabinSecret::Sharing()
{
	PartiesBuffers.resize(0, std::vector<unsigned char>());
	PartiesBuffers.resize(n_parties, std::vector<unsigned char>());
	int num_itr = (SecretBuffer.size() / k_threshold);

	for (int i = 0; i < num_itr; ++i)
	{
		for (int j = 0; j < n_parties; ++j)
		{
			unsigned char code_num = 0;
			for (int k = 0; k < k_threshold; ++k)
			{
				unsigned char temp = GFMul(SecretBuffer[i*k_threshold + k], MatrixCoefficients[j][k]);
				code_num = GFAdd(temp, code_num);
			}

			PartiesBuffers[j].push_back(code_num);
		}
	}
}

void RabinSecret::RabinSharing()
{
	//Download Secret File
	MatrixRndGenerator();
	unsigned char gab= 0;
	int chunk = k_threshold * 1024;
	
	
	std::ifstream OpenedFile(secret_file_name, std::ios::binary);

	std::vector<std::ofstream> OpenedFiles_Pat;
	
	bool flag = true;

	for (int i = 0; i < n_parties; ++i)
	{
		OpenedFiles_Pat.push_back(std::ofstream(parties_files_names[i], std::ofstream::binary));
		
		flag=flag&&OpenedFiles_Pat[i].is_open();
	}

	if (OpenedFile.is_open()&&flag)
	{
	
		while (!OpenedFile.eof())
		{
//Read chunk of File into buffer
			SecretBuffer.resize(0);
			char* buf = new char[chunk];
			
			OpenedFile.read(buf,chunk);
			
			int s = OpenedFile.gcount();
			SecretBuffer.insert(SecretBuffer.begin(), buf, buf + s);
			delete[] buf;
		
		if (SecretBuffer.size() < chunk)
		{
			int add= gab=SecretBuffer.size() % k_threshold;
			while (add>0)
			{
				SecretBuffer.push_back(0);
				--add;
			}
		}  
//Sharing chank into n chanks		
			Sharing();
//Write chunks from buffers into Files
			int i = 0;
			for (std::vector<std::vector<unsigned char>>::iterator it = PartiesBuffers.begin(); it != PartiesBuffers.end(); ++it)
			{
				std::ostream_iterator<unsigned char> file_itr(OpenedFiles_Pat[i]);
				if (flag)
				std::copy(MatrixCoefficients[i].begin(), MatrixCoefficients[i].end(), file_itr);

				std::copy(it->begin(), it->end(), file_itr);

				++i;

			}
			
			flag = false;
		}
		for (int i = 0; i < n_parties; ++i)
		{
			OpenedFiles_Pat[i] << gab;
			OpenedFiles_Pat[i].close();
		}
	
		OpenedFile.close();
	}
	else
	{
		std::cout << "can't open file";

	}
	}
void RabinSecret::Recovery()
{
	SecretBuffer.resize(0);

	for (int j = 0; j < PartiesBuffers[0].size(); ++j)
	{
		for (int str = 0; str < k_threshold; ++str)
		{
			unsigned char secret_num = 0;
			for (int i = 0; i < k_threshold; ++i)
			{
				unsigned char mul = GFMul(PartiesBuffers[i][j], MinorMatrix[str][i]);
				secret_num = GFAdd(mul, secret_num);
			}

			SecretBuffer.push_back(secret_num);
		}
	}

}

void RabinSecret::RabinRecovery()
{
	MatrixCoefficients.resize(k_threshold,std::vector<unsigned char>());
	GFInitInverse(inverses);

	int gab = 0;
 	int chunk = k_threshold * 1024;
	
	std::ofstream OpenedFile(secret_file_name, std::ofstream::binary);
	std::ostream_iterator<unsigned char> file_itr(OpenedFile);

	std::vector<std::ifstream> OpenedFiles_Pat;
	bool flag = true;

	for (int i = 0; i < k_threshold; ++i)
	{
		OpenedFiles_Pat.push_back(std::ifstream(parties_files_names[i], std::ios::binary));
		flag = flag&&OpenedFiles_Pat[i].is_open();
	}

	if (OpenedFile.is_open() && flag)
	{
		int iteration = 0;
		while (!OpenedFiles_Pat[0].eof())
		{
			PartiesBuffers.resize(0, std::vector<unsigned char>());
			PartiesBuffers.resize(k_threshold, std::vector<unsigned char>());
			//Read chunks of Files into buffers
			for (int f = 0; f < k_threshold; ++f)
			{
				char* buf = new char[chunk];

				OpenedFiles_Pat[f].read(buf, chunk);

				int s = OpenedFiles_Pat[f].gcount();
				if (iteration==0)
				{
					if (OpenedFiles_Pat[f].eof())
					{ 
						gab = buf[s-1];
						PartiesBuffers[f].insert(PartiesBuffers[f].begin(), buf + k_threshold, buf + s - 1);
					}
					else
					PartiesBuffers[f].insert(PartiesBuffers[f].begin(), buf+k_threshold, buf + s);
					MatrixCoefficients[f].insert(MatrixCoefficients[f].begin(),buf,buf+k_threshold);
				}
				else if (OpenedFiles_Pat[f].eof())
				{
					PartiesBuffers[f].insert(PartiesBuffers[f].begin(), buf, buf + s-1);
					gab = buf[s-1];
				}
				else
				PartiesBuffers[f].insert(PartiesBuffers[f].begin(), buf, buf + s);
				delete[] buf;
			}
			if (iteration == 0)GetRecoveryMatrix();
			
			//Recovery chank from k chanks		
			Recovery();
			//Write chunk from buffer into File
				
				std::copy(SecretBuffer.begin(), SecretBuffer.end(), file_itr);

		 ++iteration;
		}
		
		for (int i = 0; i < k_threshold; ++i)
		{
			OpenedFiles_Pat[i].close();
		}
		OpenedFile.close();
	}
	else
	{
		std::cout << "can't open file";
	}

}