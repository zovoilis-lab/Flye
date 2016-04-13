#include "Worker.h"
#include <chrono>
#include <sstream>

namespace
{

std::vector<std::string> 
splitString(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}
}

Worker::Worker(const std::string& scoreMatPath):
	_scoreMat(5, 5),
	_filePos(0)
{
	std::ofstream file;
	file.open("results.txt");
	std::chrono::time_point<std::chrono::system_clock> now;
	now = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(now);
	file << "File was produced at: " << std::ctime(&time);
	file << "\n";
	file.close();

	//Parse scoring matrix
	_scoreMat.loadMatrix(scoreMatPath);
}

void Worker::run(const std::string& dataPath, const std::string& format) 
{
	//Parse Fasta file (DNA)
	this->readFasta(_reads, dataPath);

	//Find the candidate
	std::string candidate = _reads[0];
	Record rec;
	std::string prev_candidate = "";
	this->outputSeparator();
	
	while(candidate.compare(prev_candidate)) {
		prev_candidate = candidate;
		std::transform(prev_candidate.begin(), prev_candidate.end(), 
					   prev_candidate.begin(), ::tolower);
		this->runOneToAll(candidate, rec);
		candidate = rec.read;
		if (format == "verbose")
			this->outputRecord(rec);
	}

	//Record the rec
	if (format == "short")
		this->outputRecord(rec);

	this->outputSeparator();
}

void Worker::run(size_t start, size_t stop, const std::string& dataPath, 
				 const std::string& format) 
{
	size_t initial = start;
	double interval = 5;
	double prev = 0;

	while (start < stop) {
		this->progressUpdate(start, stop, initial, interval, prev);
		
		//Parse Fasta file (DNA)
		//Find the candidate
		std::string candidate = this->readFastaSpecial(_reads, dataPath);
		if (_reads.size() == 0) 
		{
			return;
		}

		Record rec;
		std::string prev_candidate = "";

		outputSeparator(); //---------

		while (candidate.compare(prev_candidate)) 
		{
			prev_candidate = candidate;
			this->runOneToAll(candidate, rec);
			candidate = rec.read;
			if (format == "verbose") 
				outputRecord(rec);
		}

		//Record the rec
		if (format == "short")
			outputRecord(rec);

		outputSeparator(); //---------
		start++;
	}
}

void Worker::runOneToAll(const std::string& candidate, Record& rec) 
{
	double score = 0;
	Alignment align(_reads.size());
	//Global
	//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
	for (size_t i = 0; i < _reads.size(); ++i) 
	{
		score += align.globalAlignment(candidate, _reads[i], &_scoreMat, i);
	}

	rec.methodUsed = "global";
	rec.score = score;
	rec.read = candidate;
	std::transform(rec.read.begin(), rec.read.end(), 
				   rec.read.begin(), ::tolower);

	//Test
	//std::cout << "Global: " << candidate << " score: " << score << std::endl;

	//Deletion
	for (size_t del_index = 0; del_index < candidate.size(); del_index++) 
	{
		score = 0;

		//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
		for (size_t i = 0; i < _reads.size(); i++) {
			score += align.addDeletion(i, del_index + 1);
		}

		//Test-------------------------------------------------
		//std::string strI = candidate;
		//std::cout << "Deletion: " << strI.erase(del_index, 1) 
		//			<< " score: " << score << std::endl;
		//Test-------------------------------------------------

		//Record if less
		if (score > rec.score) {
			std::string str = candidate;
			rec.methodUsed = "deletion";
			rec.score = score;
			rec.read = str.erase(del_index, 1);
			rec.del_index = del_index;
		}
	}

	//Substitution
	char alphabet[4] = {'A', 'C', 'G', 'T'};
	for (size_t sub_index = 0; sub_index < candidate.size(); sub_index++) 
	{
		//for each (char letter in alphabet) {
		for (char letter : alphabet)
		{
			if (letter == toupper(candidate[sub_index]))
				continue;
			score = 0;

			//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
			for (size_t i = 0; i < _reads.size(); i++) {
				score += align.addSubstitution(i, sub_index + 1, letter, 
											   _reads[i], &_scoreMat);
			}

			//Test-------------------------------------------------
			//std::string strI = candidate;
			//strI.erase(sub_index, 1);
			//std::cout << "Substitution: " << strI.insert(sub_index, 1, letter) 
			//			<< " score: " << score << std::endl;
			//Test-------------------------------------------------

			//Record if less
			if (score > rec.score) 
			{
				std::string str = candidate;
				rec.methodUsed = "substitution";
				rec.score = score;
				str.erase(sub_index, 1);
				rec.read = str.insert(sub_index, 1, letter);
				rec.sub_index = sub_index;
				rec.sub_letter = letter;
			}
		}
	}

	//Insertion
	for (size_t ins_index = 0; ins_index < candidate.size()+1; ins_index++) 
	{
		for (char letter : alphabet)
		{
			score = 0;

			//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
			for (size_t i = 0; i < _reads.size(); i++) 
			{
				score += align.addInsertion(i, ins_index + 1, letter, 
											_reads[i], &_scoreMat);		
			}

			//Test-------------------------------------------------
			//std::string strI = candidate;
			//std::cout << "Insertion: " << strI.insert(ins_index, 1, letter) 
			//			<< " score: " << score << std::endl;
			//Test-------------------------------------------------


			//Record if less
			if (score > rec.score) 
			{
				std::string str = candidate;
				rec.methodUsed = "insertion";
				rec.score = score;
				rec.read = str.insert(ins_index, 1, letter);
				rec.ins_index = ins_index;
				rec.ins_letter = letter;
			}
		}
	}	
}

void Worker::outputRecord(const Record& rec) {
	std::ofstream file;
	file.open("results.txt", std::ios::app);

	file << std::fixed
		 << std::setw(22) << std::left << "Consensus: " 
		 << std::right << rec.read << "\n"
		 << std::setw(22) << std::left << "Score: " << std::right 
		 << std::setprecision(2) << rec.score << "\n"
		 << std::setw(22) << std::left << "Last method applied: " 
		 << std::right << rec.methodUsed << "\n";
	

	if (rec.methodUsed == "deletion") {
		file << "Char at index: " << rec.del_index << " was deleted. \n";
	}
	else if (rec.methodUsed == "substitution") {
		file << "Char at index " << rec.sub_index << " was substituted with " 
			<< "'" << rec.sub_letter << "'" << ".\n";
	}
	else if (rec.methodUsed == "insertion") {
		file << "'"<< rec.ins_letter << "'" 
			 << " was inserted at index " << rec.ins_index << ".\n";
	}
	file << std::endl;
	file.close();
}

void Worker::outputSeparator() 
{
	std::ofstream file;
	file.open("results.txt", std::ios::app);
	file << "------------------------------------------ \n";
	file.close();
}

void Worker::readFasta(std::vector<std::string>& reads, 
					   const std::string& path) 
{
	std::string line;
	std::ifstream file(path);

	if (file.is_open()) {
		while (getline(file, line)) {
			if (line[0] == '>' || line == "") {
				continue;
			}
			reads.push_back(line);
		}
	}
}

std::string Worker::readFastaSpecial(std::vector<std::string>& reads, 
									 const std::string& path)
{
	std::string line;
	std::ifstream file(path, std::ios::binary);
	std::string candidate;
	reads.clear();

	if (file.is_open()) {
		//Get to the current bubble
		file.seekg(_filePos, file.beg);

		getline(file, line);
		if (line == "") {
			std::cerr << "End of the file reached. 'end' value might exceed number of bubbles. \n";
			return "";
		}
		std::vector<std::string> elems = splitString(line, ' ');
		
		std::string marker = elems[0];
		//int index = stoi(elems[1]);
		int numOfReads = stoi(elems[2]);
		//int sizeOfCand = stoi(elems[3]);

		if (marker.compare(">current") == 0) {
			getline(file, candidate);
		}

		//Get all appropriate string for the current bubble
		int count = 0;
		while (getline(file, line) && count < numOfReads) {
			if (line[0] == '>' || line == "") {
				continue;
			}
			reads.push_back(line);
			count++;
			_filePos = file.tellg();
		}
	}
	file.close();
	return(candidate);
}

//std::vector<std::string>& 
//Worker::split(const std::string &s, char delim, std::vector<std::string> &elems) {
//	
//	return elems;
//}

void Worker::progressUpdate(int start, int stop, int initial, 
							int interval, double& prev) 
{
	double done = start - initial;
	double total = stop - initial;
	double progress = done / total * 100;
	if (progress >= prev) {
		std::cout << "Completed: " << std::setw(2) 
				  << progress << "%" << std::endl;
		prev += interval;
	}
}
