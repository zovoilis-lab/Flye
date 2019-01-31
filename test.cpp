#include "src/common/bfcontainer.h"
#include <algorithm>
#include <iostream>
#include <chrono>
#include <deque>

void dequeTest(size_t VEC_SIZE)
{
	std::cout << "deque test" << std::endl;

	auto timeStart = std::chrono::system_clock::now();
	std::deque<int> vec;
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		vec.push_back(rand() % VEC_SIZE);
	}

	std::cout << "Initialized " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();

	size_t next = 0;
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		next = vec[next];
	}

	std::cout << next << std::endl;
	std::cout << "Access " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();

	std::sort(vec.begin(), vec.end());
	std::cout << "Sorted " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
}


void vectorTest(size_t VEC_SIZE)
{
	std::cout << "vector test" << std::endl;

	auto timeStart = std::chrono::system_clock::now();
	std::vector<int> vec;
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		vec.push_back(rand() % VEC_SIZE);
	}

	std::cout << "Initialized " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();

	size_t next = 0;
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		next = vec[next];
	}

	std::cout << next << std::endl;
	std::cout << "Access " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();

	std::sort(vec.begin(), vec.end());
	std::cout << "Sorted " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
}

void bfTest(size_t VEC_SIZE)
{
	std::cout << "bfTest" << std::endl;
	typedef int TestType;
	ChunkPool<TestType> pool;
	BFContainer<TestType> container(pool);

	auto timeStart = std::chrono::system_clock::now();
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		container.push_back(rand() % VEC_SIZE);
	}

	std::cout << "Initialized " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();

	size_t next = 0;
	for (size_t i = 0; i < VEC_SIZE; ++i)
	{
		next = container[next];
	}

	std::cout << next << std::endl;
	std::cout << "Access " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
	timeStart = std::chrono::system_clock::now();


	std::sort(container.begin(), container.end());
	for (size_t i = 0; i < container.size() - 1; ++i)
	{
		if (container[i] > container[i + 1]) throw std::runtime_error("failure!");
	}

	std::cout << "Sorted " << 
		std::chrono::duration_cast<std::chrono::duration<float>>
			(std::chrono::system_clock::now() - timeStart).count() << std::endl;
}

int main()
{
	size_t SIZE = 10000000;
	vectorTest(SIZE);
	dequeTest(SIZE);
	bfTest(SIZE);
}
