#include<iostream>
#include <fstream>
#include<cmath>
#include<omp.h>
#include<iomanip>
#include<cmath>
#define Count_of_Threads  8
void Insert_In_File(double* new_array, size_t count, std::fstream& stream)
{
	if (stream.is_open())
	{
		for (size_t i = 0; i < count; i++)
			stream << new_array[i] << std::setprecision(16) << " ";
		stream << "\n";
		stream << "\n";
	}
	else
		throw "The file isn't exist";
}

inline double Get_A_Part_Of_Approximation(double prev_point, double curr_point, double next_point, double deltaX, double sigma, double K_i_prev, double K_i_curr,double K_i_next, double TimeCoef)
{
	
	return(  (sigma/(deltaX*deltaX)) * (prev_point+(K_i_prev*TimeCoef) - 2*(curr_point)-2*(K_i_curr*TimeCoef)+next_point+(K_i_next*TimeCoef) )   );
}

int main()
{
	double X_max = 0.0;
	double deltaX = 0.0;
	double T_max = 0.0;
	double deltaT = 0.0;
	double sigma = 0.0;
	size_t count = 0;
	double time = 0.0;
	double tmp = 0.0;
	const double PI = 3.141592653589793;
	int i = 1;


	std::fstream f;
	f.open("output.txt", std::fstream::in | std::fstream::out);

	std::ifstream f_in;
	f_in.open("const_initial.txt");
	if (f_in.is_open())
	{
		f_in >> X_max;
		f_in >> deltaX;
		f_in >> T_max;
		f_in >> deltaT;
		f_in >> sigma;
	}
	else
		throw "file wasn't opened";
	f_in.close();

	if (X_max == 0.0 || deltaX == 0.0)
		throw "count = 0 or count is infinity";

	count = static_cast<size_t>((int)(X_max / deltaX)) + 1;

	double* temp_array = nullptr;
	temp_array = new double[count];

	double* curr_array = nullptr;
	curr_array = new double[count];

	double* tmp_buffer = nullptr;
	tmp_buffer = new double[count];

	double* K1 = new double[count];
	double* K2 = new double[count];
	double* K3 = new double[count];
	double* K4 = new double[count];

	double t1 = 0, t2 = 0, dt = 0;


	for (size_t q = 0; q < count; q++)
		temp_array[q] = sin(PI * deltaX * q);

	temp_array[count - 1] = 0.0;


	Insert_In_File(temp_array, count, f);

	t1 = omp_get_wtime();
	int a = 1;

	while (time <= T_max)
	{
		time = deltaT * i;
		curr_array[0] = 0.0;
		curr_array[count - 1] = 0.0;
	
		K1[0] = 0;
		K1[count - 1] = 0;
		K2[0] = 0;
		K2[count - 1] = 0;
		K3[0] = 0;
		K3[count - 1] = 0;
		K4[0] = 0;
		K4[count - 1] = 0;
	
#pragma omp parallel for shared(K1,temp_array)  schedule (static,Count_of_Threads )
		for (int a = 1; a < count - 1; a++)
			K1[a] = Get_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, 0, 0, 0, 0);

#pragma omp parallel for shared(K1,K2,temp_array)  schedule (static,Count_of_Threads )
		for (int a = 1; a < count - 1; a++)
			K2[a] = Get_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1,K1[a - 1], K1[a], K1[a + 1], (deltaT / 2));

#pragma omp parallel for shared(K2,K3,temp_array) schedule (static,Count_of_Threads )
		for (int a = 1; a < count - 1; a++)
			K3[a] = Get_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, K2[a - 1], K2[a], K2[a + 1], (deltaT / 2));

#pragma omp parallel for shared(K3,K4,temp_array)  schedule (static,Count_of_Threads )
		for (int a = 1; a < count - 1; a++)
			K4[a] = Get_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, K3[a - 1], K3[a], K3[a + 1], deltaT );

#pragma omp parallel for shared(temp_array,curr_array,K1,K2,K3,K4) schedule (static,Count_of_Threads )
		for (int a = 1; a < count - 1; a++)
			curr_array[a] = temp_array[a] + (deltaT / 6) * (K1[a] + (2 * K2[a]) + (2 * K3[a]) + K4[a]);
		

			
		
		tmp_buffer = temp_array;
		temp_array = curr_array;
		curr_array = tmp_buffer;

		
		i++;
	}
	
	t2 = omp_get_wtime();
	dt = t2 - t1;
	std::cout << dt << std::endl;
	
	Insert_In_File(curr_array, count, f);
	f.close();
	delete[] temp_array;
	
	delete[] curr_array;	
	delete[] K1;
	delete[] K2;
	delete[] K3;
	delete[] K4;

	return 0;
}

