#pragma once

#include "../../lib/minimap2/minimap.h"

class MinimapBuffer
{
public:
	MinimapBuffer();
	~MinimapBuffer();

	MinimapBuffer(const MinimapBuffer&) = delete;
	void operator=(const MinimapBuffer&) = delete;

	mm_tbuf_t* get() const;

private:
	mm_tbuf_t* _buffer;
};