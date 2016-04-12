#include "Worker.h"
#include <chrono>


Worker::Worker(string scoreMatPath){
	ofstream file;
	file.open("results.txt");
	std::chrono::time_point<std::chrono::system_clock> now;
	now = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(now);
	file << "File was produced at: " << std::ctime(&time);
	file << "\n";
	file.close();

	//Parse scoring matrix
	scoreMat = new ScoringMatrix(5, 5);
	scoreMat->loadMatrix(scoreMatPath);
	filePos = 0;
}

void Worker::run(string dataPath, string format) {
	//Parse Fasta file (DNA)
	readFasta(reads, dataPath);

	//Find the candidate
	std::string candidate = reads[0];
	record rec;
	std::string prev_candidate = "";
	outputSeparator();
	
	while(candidate.compare(prev_candidate)) {
		prev_candidate = candidate;
		std::transform(prev_candidate.begin(), prev_candidate.end(), 
					   prev_candidate.begin(), ::tolower);
		runOneToAll(candidate, rec);
		candidate = rec.read;
		if (format == "verbose") {
			outputRecord(rec);
		}
	}

	//Record the rec
	if (format == "short") {
		outputRecord(rec);
	}
	outputSeparator();
}

void Worker::run(size_t start, size_t stop, string dataPath, string format) {
	size_t initial = start;
	double interval = 5;
	double prev = 0;

	while (start < stop) {
		progressUpdate(start, stop, initial, interval, prev);
		
		//Parse Fasta file (DNA)
		//Find the candidate
		std::string candidate = readFastaSpecial(reads, dataPath);
		if (reads.size() == 0) {
			return;
		}

		record rec;
		std::string prev_candidate = "";

		outputSeparator(); //---------

		while (candidate.compare(prev_candidate)) {
			prev_candidate = candidate;
			runOneToAll(candidate, rec);
			candidate = rec.read;
			if (format == "verbose") {
				outputRecord(rec);
			}
		}

		//Record the rec
		if (format == "short") {
			outputRecord(rec);
		}

		outputSeparator(); //---------
		start++;
	}
}

void Worker::runOneToAll(std::string candidate, record& rec) {
	double score = 0;
	Alignment align(reads.size());
	//Global
	//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
	for (int i = 0; i < reads.size(); ++i) {
		score += align.globalAlignment(candidate, reads[i], scoreMat, i);
	}

	rec.methodUsed = "global";
	rec.score = score;
	rec.read = candidate;
	std::transform(rec.read.begin(), rec.read.end(), rec.read.begin(), ::tolower);

	//Test
	//std::cout << "Global: " << candidate << " score: " << score << std::endl;

	//Deletion
	for (int del_index = 0; del_index < candidate.size(); del_index++) {
		score = 0;

		//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
		for (int i = 0; i < reads.size(); i++) {
			score += align.addDeletion(i, del_index + 1);
		}

		//Test-------------------------------------------------
		//std::string strI = candidate;
		//std::cout << "Deletion: " << strI.erase(del_index, 1) << " score: " << score << std::endl;
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
	char alphabet[4] = { 'A', 'C', 'G', 'T' };
	for (int sub_index = 0; sub_index < candidate.size(); sub_index++) {
		//for each (char letter in alphabet) {
		for (char letter : alphabet)
		{
			if (letter == toupper(candidate[sub_index])) {
				continue;
			}
			score = 0;

			//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
			for (int i = 0; i < reads.size(); i++) {
				score += align.addSubstitution(i, sub_index + 1, letter, reads[i], scoreMat);
			}

			//Test-------------------------------------------------
			//std::string strI = candidate;
			//strI.erase(sub_index, 1);
			//std::cout << "Substitution: " << strI.insert(sub_index, 1, letter) << " score: " << score << std::endl;
			//Test-------------------------------------------------

			//Record if less
			if (score > rec.score) {
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
	for (int ins_index = 0; ins_index < candidate.size()+1; ins_index++) {
		for (char letter : alphabet)
		{
			score = 0;

			//#pragma omp parallel for schedule(dynamic) reduction(+ : score)
			for (int i = 0; i < reads.size(); i++) {
				score += align.addInsertion(i, ins_index + 1, letter, reads[i], scoreMat);		
			}

			//Test-------------------------------------------------
			//std::string strI = candidate;
			//std::cout << "Insertion: " << strI.insert(ins_index, 1, letter) << " score: " << score << std::endl;
			//Test-------------------------------------------------


			//Record if less
			if (score > rec.score) {
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

void Worker::outputRecord(record rec) {
	ofstream file;
	file.open("results.txt", ios::app);

	file << fixed;
	file << setw(22) << left << "Consensus: " << right << rec.read << "\n";
	file << setw(22) << left << "Score: " <<right << std::setprecision(2) << rec.score << "\n";
	file << setw(22) << left << "Last method applied: " << right << rec.methodUsed << "\n";
	

	if (rec.methodUsed == "deletion") {
		file << "Char at index: " << rec.del_index << " was deleted. \n";
	}
	else if (rec.methodUsed == "substitution") {
		file << "Char at index " << rec.sub_index << " was substituted with " 
			<<"'"<< rec.sub_letter <<"'"<< ".\n";
	}
	else if (rec.methodUsed == "insertion") {
		file << "'"<<rec.ins_letter <<"'"<< " was inserted at index " << rec.ins_index << ".\n";
	}
	file << "\n";
	//file << "\n";
	file.close();
}

void Worker::outputSeparator() {
	ofstream file;
	file.open("results.txt", ios::app);
	file << "------------------------------------------ \n";
	file.close();
}

void Worker::readFasta(vector<string>& reads, string path) {
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

string Worker::readFastaSpecial(vector<string>& reads, string path) {
	std::string line;
	std::ifstream file(path, std::ios::binary);
	string candidate;
	reads.clear();

	if (file.is_open()) {
		//Get to the current bubble
		file.seekg(filePos, file.beg);

		getline(file, line);
		if (line == "") {
			std::cerr << "End of the file reached. 'end' value might exceed number of bubbles. \n";
			return "";
		}
		vector<string> elems = split(line, ' ');
		
		string marker = elems[0];
		int index = stoi(elems[1]);
		int numOfReads = stoi(elems[2]);
		int sizeOfCand = stoi(elems[3]);

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
			filePos = file.tellg();
		}
	}
	file.close();
	return(candidate);
}

std::vector<std::string>& Worker::split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

void Worker::progressUpdate(int start, int stop, int initial, int interval, double& prev) {
	double done = start - initial;
	double total = stop - initial;
	double progress = done / total * 100;
	if (progress >= prev) {
		std::cout << "Completed: " << setw(2) << progress << "%" << std::endl;
		prev += interval;
	}
}

std::vector<std::string> Worker::split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

Worker::~Worker() {
	delete scoreMat;
}
