#include<iostream>
#include <fstream>
#include<cmath>
void Insert_In_File(double* new_array, size_t count, std::fstream& stream)
{
	if (stream.is_open())
	{
		for (size_t i = 0; i < count; i++)
			stream << new_array[i] << " ";
		stream << "\n";
		stream << "\n";
	}
	else
		throw "The file isn't exist";
}

inline double Get_A_Part_Of_Approximation(double prev_point, double curr_point, double next_point, double deltaX, double sigma, double K_i, double TimeCoef)
{
	double tmp_coef_with_K = K_i * TimeCoef;
	return(  (sigma/(deltaX*deltaX)) * (prev_point+tmp_coef_with_K - 2*(curr_point+tmp_coef_with_K)+next_point+tmp_coef_with_K  ));
}
inline double Get_Another_A_Part_Of_Approximation(double prev_point, double curr_point, double next_point, double deltaX, double sigma, double K_i, double TimeCoef)
{
	double tmp_coef_with_K = K_i * TimeCoef;
	return( tmp_coef_with_K +  (sigma / (deltaX * deltaX)) * (prev_point - 2*(curr_point)+next_point));
}

inline double Get_Another_Another_A_Part_Of_Approximation(double prev_point, double curr_point, double next_point, double deltaX, double sigma, double K_i, double TimeCoef)
{
	double tmp_coef_with_K = K_i * TimeCoef;
	return( (sigma / (deltaX * deltaX)) * (prev_point - 2 * (curr_point)+next_point+tmp_coef_with_K));
}

inline double Get_New_Point_In_Time(double prev_point, double curr_point, double next_point, double deltaT, double deltaX, double sigma)
{

	double tmp = 0.0;
	tmp = prev_point - 2 * curr_point + next_point;

	return((sigma * deltaT) / (deltaX * deltaX) * tmp + curr_point);
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

	


	for (size_t q = 0; q < count; q++)
		temp_array[q] = sin(PI * deltaX * q);

	temp_array[count - 1] = 0.0;


	Insert_In_File(temp_array, count, f);

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

		for (int a = 1; a < count - 1; a++)
		{
			K1[a] = Get_Another_Another_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, 0, 0);
			K2[a] = Get_Another_Another_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, K1[a],  (deltaT / 2));
			K3[a] = Get_Another_Another_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, K2[a], (deltaT / 2));
			K4[a] = Get_Another_Another_A_Part_Of_Approximation(temp_array[a - 1], temp_array[a], temp_array[a + 1], deltaX, 1, K3[a], deltaT);
			curr_array[a] = temp_array[a] + (deltaT / 6) * (K1[a] + (2 * K2[a]) + (2 * K3[a]) + K4[a]);
			
		}
		if (i % 100000 == 0)
		{
			//Insert_In_File(curr_array, count, f);
			std::cout << time << std::endl;
		}
		//for (size_t q = 0; q < count; q++)
		//	temp_array[q] = curr_array[q];

		tmp_buffer = temp_array;
		temp_array = curr_array;
		curr_array = tmp_buffer;

		//std::cout << time << std::endl;
		i++;
	}

	Insert_In_File(curr_array, count, f);
	f.close();
	delete[] temp_array;
	//delete[] tmp_buffer;
	delete[] curr_array;	
	delete[] K1;
	delete[] K2;
	delete[] K3;
	delete[] K4;

	return 0;
}

