//(c) 2019 by Authors
//This file is a part of the Flye program.
//Released under the BSD license (see LICENSE file)

//Big Functional Container - a substitute for vector
//for very large arrays (up to several billion elements).
//Basically, it is a vector of vectors - so we have 
//a list of (relatively big) chunks instead of one
//gigantic memory chunk (several Gb) - allocating which might be
//difficult for OS, cause memory fragmentation and slowdown.
//Includes memory pool of chunks that could be shared by multuple threads

#pragma once

#include <vector>
#include <mutex>

template <class T, int ChunkSize = 32 * 1024 * 1024>
class ChunkPool
{
public:
	ChunkPool(){}
	ChunkPool(const ChunkPool& other) = delete;
	ChunkPool(ChunkPool&& other) = delete;
	ChunkPool& operator=(ChunkPool& other) = delete;

	~ChunkPool()

	{
		for (T* chunk : _freeChunks) delete[] chunk;
		_freeChunks.clear();
	}

	T* getChunk()
	{
		std::lock_guard<std::mutex> lock(_chunkMutex);
		if (!_freeChunks.empty())
		{
			T* ref = _freeChunks.back();
			_freeChunks.pop_back();
			return ref;
		}
		else
		{
			return new T[chunkCapacity()];
		}
	}

	void returnChunk(T* chunk)
	{
		std::lock_guard<std::mutex> lock(_chunkMutex);
		_freeChunks.push_back(chunk);
	}

	size_t chunkCapacity()
	{
		return ChunkSize / sizeof(T);
	}

	size_t numberChunks() {return _freeChunks.size();}
	
private:
	std::mutex 		_chunkMutex;
	std::vector<T*> _freeChunks;
};


template <class T>
class BFContainer
{
public:
	BFContainer(ChunkPool<T>& chunkPool): 
		_pool(chunkPool), _chunkCapacity(chunkPool.chunkCapacity()),
		_size(0), _lastChunkOffset(_chunkCapacity) {}

	BFContainer(const BFContainer& other) = delete;
	BFContainer(BFContainer&& other) = delete;
	BFContainer& operator=(const BFContainer& other) = delete;

	~BFContainer()
	{
		this->clear();
	}

	void push_back(T&& elem)
	{
		if (_lastChunkOffset == _chunkCapacity)
		{
			_chunks.push_back(_pool.getChunk());
			_lastChunkOffset = 0;
		}
		_chunks.back()[_lastChunkOffset++] = elem;
		++_size;
	}

	template <typename ... Args>
	void emplace_back(Args&& ... args)
	{
		if (_lastChunkOffset == _chunkCapacity)
		{
			_chunks.push_back(_pool.getChunk());
			_lastChunkOffset = 0;
		}
		new (_chunks.back() + _lastChunkOffset) T(std::forward<Args>(args)...);
		++_lastChunkOffset;
		++_size;
	}

	size_t size() {return _size;}

	T& operator[](size_t index)
	{
		return _chunks[index / _chunkCapacity][index % _chunkCapacity];
	}

	void clear()
	{
		for (T* chunk : _chunks)
		{
			_pool.returnChunk(chunk);
			_chunks.clear();
			_lastChunkOffset = _chunkCapacity;
		}
	}

	class BFIterator: public std::iterator<std::random_access_iterator_tag, T>
	{
	private:
		BFContainer* _container;
		size_t _index;

    public:
        BFIterator(BFContainer* cont=nullptr, size_t idx=0): 
			_container(cont), _index(idx) {}

        BFIterator& operator=(const BFIterator &rhs) 
			{_container = rhs._container; _index = rhs._index; return *this;}
        BFIterator& operator+=(size_t rhs) {_index += rhs; return *this;}
        BFIterator& operator-=(size_t rhs) {_index -= rhs; return *this;}
		T& operator*() {return (*_container)[_index];}
        //T* operator->() {return &(*_container)[_index];}
		T& operator[](size_t rhs) {return _container[_index + rhs];}

    	bool operator==(const BFIterator& rhs) {return _index == rhs._index;}
        bool operator!=(const BFIterator& rhs) {return !(*this == rhs);}
        bool operator> (const BFIterator& rhs) {return _index > rhs._index;}
        bool operator< (const BFIterator& rhs) {return _index < rhs._index;}
        bool operator>=(const BFIterator& rhs) {return _index >= rhs._index;}
        bool operator<=(const BFIterator& rhs) {return _index <= rhs._index;}
		
		BFIterator& operator++() {++_index; return *this;}
        BFIterator& operator--() {--_index; return *this;}
        BFIterator& operator++(int) {BFIterator tmp(*this); ++_index; return tmp;}
        BFIterator& operator--(int) {BFIterator tmp(*this); --_index; return tmp;}

        BFIterator operator+(size_t rhs) 
			{return BFIterator(_container, _index + rhs);}
		friend BFIterator operator+(size_t lhs, const BFIterator& rhs) 
			{return rhs + lhs;}
		
        BFIterator operator-(size_t rhs) 
			{return BFIterator(_container, _index - rhs);}
        size_t operator-(const BFIterator& rhs) 
			{return _index - rhs._index;}

	};

	BFIterator begin() {return BFIterator(this, 0);}
	BFIterator end() {return BFIterator(this, _size);}

private:

	ChunkPool<T>& 	_pool;
	std::vector<T*> _chunks;
	size_t 			_chunkCapacity;
	size_t 			_size;
	size_t 			_lastChunkOffset;
};

