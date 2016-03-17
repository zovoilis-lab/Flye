
#include <string>

class Kmer
{
public:
	Kmer(const std::string& kmerString);

	Kmer reverseComplement() const;
	Kmer rotateRight(char dnaSymbol) const;
	Kmer rotateLeft(char dnaSymbol) const;
	std::string dnaRepresentation() const;

private:
	int64_t _representation;
};
