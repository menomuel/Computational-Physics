#include <iostream>
#include <limits>
#include <iomanip>

template<typename T>
T machineEpsilon(bool printFlag = false)
{
	auto eps = T(1);
	auto count = 0;
	while (1 + eps / 2 > 1)
	{
		eps /= 2;
		++count;
	}
	if(printFlag)
		std::cout << "Had to divide by 2 " << count << " times\n";
	return eps;
}


template<typename T>
void bitOrder()
{
	int typeSize = sizeof(T) * 8;
	
	auto eps = T(1);
	auto countMantissa = 0;
	while (1 + eps / 2 > 1)
	{
		eps /= 2;
		++countMantissa;
	}

	int countPow = typeSize - 1 - countMantissa;


	std::cout << "Type size - " << typeSize << ", where sign " << 1
													<< ", pow " << countPow
													<< ", mantissa " << countMantissa << std::endl;	
}

template<typename T>
int maxPow()
{
	auto num = T(1);
	int pow = 0;

	while (num * 10 != num)
	{
		num *= 10;
		++pow;
	}

	return --pow; // cut 1 power for inf
}


//subnormal
template<typename T>
int minPow()
{
	auto num = T(1);
	int pow = 0;

	while (num / 10 != num)
	{
		num /= 10;
		--pow;
	}

	return pow;
}

/*
template <class T>
int maxMantissaPow() {
	T add = T(1);
	int count = 0;
	
	do {
		++count;
		add /= 10;
	} while (1 != 1 + add);
	return count;
}

auto num = T(1);
	int pow = 0, countPow = 0;
	while (num * 2 != num)
	{
		num *= 2;
		++pow;
	}

	while(pow != 1)
	{
		pow /= 2;
		++countPow;
	}
*/



int main()
{
	using test_type = float;
	bitOrder<test_type>();
	
	std::cout << "std: " << std::numeric_limits<test_type>::epsilon() << std::endl;
	std::cout << "calc: " << machineEpsilon<test_type>() << std::endl;
	std::cout << "min/max pow: " << maxPow<test_type>() << std::endl;
	std::cout << "denorm min pow: " << minPow<test_type>() << std::endl;

	std::cout << std::setprecision(maxPow<test_type>());
	std::cout << test_type(1) << std::endl;
	std::cout << 1 + machineEpsilon<test_type>() / 2 << std::endl;
	std::cout << 1 + machineEpsilon<test_type>() << std::endl;
	std::cout << 1 + machineEpsilon<test_type>() + machineEpsilon<test_type>() / test_type(2) << std::endl;

	
	return 0;
}